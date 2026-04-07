"""
Example Usage
----------------
The following filters will return cis-bidentate N-O donors with up to 50 atoms that do not contain any CH2 groups, only contain C, H, N, O atoms in total, and have been observed to coordinate to Pt+2, Pt+4, Pd, or Ni metal centers. We speed up the filtering by only loading 1000 ligands from the MetaLig database.

.. code-block:: python

    from DARTassembler import LigandFilters
    filters = LigandFilters(db='metalig', n=1000)
    db = filters.run(
            filters=[
                {'filter': 'property', 'name': 'n_atoms', 'range': [1, 50]},
                {'filter': 'property', 'name': 'archetype', 'values': ['2-cis']},
                {'filter': 'smarts', 'smarts': '[C&H2]', 'should_contain': False, 'include_metal': True},
                {'filter': 'composition', 'elements': 'CHNO', 'instruction': 'must_at_least_contain', 'only_donors': False},
                {'filter': 'composition', 'elements': 'NO', 'instruction': 'must_contain_and_only_contain', 'only_donors': True},
                {'filter': 'parents', 'metal_centers': ['Pt+2', 'Pt+4', 'Pd', 'Ni']}
            ],
            outpath='filtered_ligand_db.jsonlines',
            dbinfo=True,
            metal=True,
    )
"""
from copy import deepcopy
from pathlib import Path
from typing import Union, Generator, Optional
import pandas as pd
from DARTassembler.src.metalig.db import LigandDB
from DARTassembler.src.metalig.utils_molecule import get_standardized_stoichiometry_from_atoms_list, \
    stoichiometry2atomslist
from DARTassembler.src.misc.io import get_correct_ligand_db_path_from_input, read_yaml
from DARTassembler.src.modules.modules import BaseModule
from DARTassembler.src.constants.paths import default_ligandfilters_yml_path

class LigandFilters(BaseModule):
    """
    This module applies user-defined filters to a ligand database to obtain a subset of ligands with desired properties.
    """

    def __init__(self, db: Union[str, Path, None], n: Union[int, None] = None):
        """
        Initialize the DART LigandFilters module. The options set here applied to all batches.

         .. tip:: All the parameters below are available as well via the ligandfilters .yml file as global options.

        :param str | None db: .jsonlines ligand db filepath or None to use the entire MetaLig database.
        :param int | None n: Maximum number of ligands to load from the database. If None, load all ligands.
        :return: None
        :rtype: None
        """
        super().__init__()
        self.n_max_ligands = n
        self.input_ligand_db_path = get_correct_ligand_db_path_from_input(db)
        self.db = LigandDB.from_json(path=self.input_ligand_db_path, n_max=n)

    def apply_filters(self, filters: list[dict]) -> LigandDB:
        """
        Apply a sequence of filters to the loaded ligand database and return the filtered database.

        Each filter in the list must be a dictionary with a key 'filter' specifying the filter type. Supported filters:

        - 'property' : filter by a named global property (e.g. 'charge', 'archetype').
        - 'composition' : filter by element composition or stoichiometry (elements may be a string that is converted to an atom list).
        - 'parents' : filter by parent metal centers (e.g. ['Fe', 'Co']).
        - 'smarts' : filter by a SMARTS pattern; note that a bond-order validity property filter is prepended automatically.

        :param filters: Ordered list of filter specification dictionaries. Each dictionary must contain at least the key 'filter' and other filter-specific keys.
        :type filters: list[dict]
        :raises ValueError: If a filter type specified in any filter dict is not recognized.
        :return: LigandDB object containing only the ligands that passed all filters.
        :rtype: LigandDB
        """
        self.filter_tracking = []
        filters = deepcopy(filters)

        # Prepend a bond-order filter if a SMARTS filter is used to make sure the SMARTS filter can be applied.
        filters = self._add_filter_for_valid_bond_orders_if_smarts(filters)

        unames = list(self.db.db.keys())
        self.n_ligands_before = len(unames)
        self.df_all_ligands = self._get_ligand_df()
        self.df_all_ligands['filter'] = None    # initialize column for filter tracking

        for idx, filter in enumerate(filters):
            n_ligands_before = len(unames)
            try:
                filtername = filter.pop('filter')
            except KeyError:
                raise ValueError(f'Filter does not have a "filter" key specifying the filter type: {filter}')
            if filtername == 'property':
                unames = [uname for uname in unames if self.db.db[uname].property_filter(**filter)]
                name_appendix = filter['name']
            elif filtername == 'composition':
                if isinstance(filter['elements'], str):
                    filter['elements'] = stoichiometry2atomslist(filter['elements'])
                unames = [uname for uname in unames if self.db.db[uname].composition_filter(**filter)]
                name_appendix = get_standardized_stoichiometry_from_atoms_list(filter['elements'])
            elif filtername == 'parents':
                unames = [uname for uname in unames if self.db.db[uname].parents_filter(**filter)]
                name_appendix = ', '.join(filter['metal_centers'])
            elif filtername == 'smarts':
                unames = [uname for uname in unames if self.db.db[uname].smarts_filter(**filter)]
                name_appendix = filter['smarts']
            else:
                raise ValueError(f'Filtername "{filtername}" is not valid! Supported filter names are "property", "composition", "parents", and "smarts".')

            n_ligands_after = len(unames)
            ligand_was_filtered = ~self.df_all_ligands.index.isin(unames) & (self.df_all_ligands['filter'].isna())
            unique_filtername = f"Filter {idx + 1:02d}: {filtername}: {name_appendix}"
            self.df_all_ligands.loc[ligand_was_filtered, 'filter'] = unique_filtername

            self.filter_tracking.append({
                "filter": filtername,
                "unique_filtername": unique_filtername,
                "n_ligands_before": n_ligands_before,
                "n_ligands_after": n_ligands_after,
                "n_ligands_removed": n_ligands_before - n_ligands_after,
                "full_filter_options": {name: option for name, option in filter.items()}
            })
        self.n_ligands_after = len(unames)

        self.df_all_ligands.fillna({'filter': 'Passed'}, inplace=True)      # fill in 'Passed' for ligands that were not filtered out
        self.df_all_ligands.set_index('unique_name', inplace=True)    # set index to ligand ID, making sure that the column in the csv is named 'Ligand ID'
        columns = ['filter'] + [col for col in self.df_all_ligands.columns if col != 'filter']
        self.df_all_ligands = self.df_all_ligands[columns]                # move 'filter' column to the front
        self.df_all_ligands = self.df_all_ligands.sort_values(by='filter')# sort by filter name

        filtered_db = self.db.get_sub_db(ligand_names=unames)

        return filtered_db

    def _add_filter_for_valid_bond_orders_if_smarts(self, filters: list[dict]):
        """
        Prepend a property filter requiring valid bond orders when any SMARTS filter is present.

        SMARTS matching requires correct bond-order information. If any filter in the provided list
        has 'filter' == 'smarts', a property filter demanding ligands with valid bond orders is
        inserted at the front of the filter list.

        :param filters: List of filter specification dictionaries to inspect and possibly modify.
        :type filters: list[dict]
        :return: Updated list of filters with a bond-order property filter prepended when necessary.
        :rtype: list[dict]
        """
        uses_smarts_filter = any(filter['filter'] == 'smarts' for filter in filters)

        if uses_smarts_filter:
            bo_filter = {
                'filter': 'property',
                'name': 'has_all_bond_orders_valid',
                'values': [True]
            }
            filters = [bo_filter] + filters

        return filters

    def _get_filter_tracking_string(self) -> str:
        """
        Produce a human-readable summary of filter actions and final statistics.

        The returned string contains:

        - a summary of filter options used,
        - a compact table of ligands passed/removed per filter,
        - overall counts before/after filtering,
        - counts of passed archetypes and sample stoichiometries,
        - filename hints for saved outputs.

        :return: Multi-line formatted summary string intended for user-facing info files or stdout.
        :rtype: str
        """
        df_filters = pd.DataFrame(self.filter_tracking)
        df_filters = df_filters[['unique_filtername', 'n_ligands_removed', 'n_ligands_after', 'full_filter_options']]
        df_filters = df_filters.rename(columns={'n_ligands_removed': 'Ligands removed', 'n_ligands_after': 'Ligands passed', 'unique_filtername': 'Filters', 'full_filter_options': 'Filter options'})
        df_filters = df_filters.set_index('Filters')

        output = f"{'  Filter Options  ':=^80}\n"
        max_colwidth = 45
        for filter, filter_options in df_filters['Filter options'].items():
            if len(filter) > max_colwidth:
                filter = filter[:max_colwidth-3] + '...'
            filter_options = ', '.join(f'{option}: {value}' for option, value in filter_options.items())
            output += f"{filter: <{max_colwidth+2}}{filter_options}\n"

        output += f"{'  Filter Results  ':=^80}\n"
        output += df_filters[['Ligands passed', 'Ligands removed']].to_string(justify='center', index_names=False, max_colwidth=max_colwidth) + '\n'

        # Count denticities of all passed ligands
        archetype_count = pd.Series([lig.archetype for lig in self._all_ligands_left()]).value_counts().to_dict()
        archetype_output = ', '.join(sorted([f'{geom} ({count})' for geom, count in archetype_count.items()]))

        n_ligands_before = len(self.df_all_ligands)
        n_ligands_after = len(self.unames)

        output += f"{'  Total summary of DART Ligand Filters run  ':=^80}\n"
        output += f"Before filtering:           {n_ligands_before} ligands\n"
        output += f"Filtered out:               {n_ligands_before - n_ligands_after} ligands\n"
        output += f"Passed:                     {n_ligands_after} ligands\n"
        output += f"Passed archetypes:          {archetype_output}\n"

        # Print the stoichiometries of the first five passed ligands.
        stoichiometries = ', '.join([lig.stoichiometry for idx, lig in enumerate(self._all_ligands_left()) if idx < 5])
        ellipsis = ', ...' if n_ligands_after > 5 else ''
        output += f'Passed stoichiometries:     {stoichiometries}{ellipsis}\n'
        if self.output_ligand_db_path is not None:
            output += f"Filtered ligand database with {n_ligands_after} ligands was saved to `{self.output_ligand_db_path.name}`.\n"
            if self.output_info:
                output += f"Info on filtered ligands saved to directory `{self.outdir.name}`.\n"
        output += "Done! All ligands filtered. Exiting DART Ligand Filters Module."

        return output

    def save_filtered_ligands_output(self) -> None:
        """
        Write informational outputs for filtered ligands: a summary text, a CSV overview, and concatenated XYZ files.

        The method creates an 'info' directory next to the specified output ligand DB path and writes:

        - filters.txt : human readable filter summary produced by _get_filter_tracking_string,
        - ligands_overview.csv : table with ligand metadata and filter assignment,
        - concat_*.xyz : concatenated xyz files for 'Passed' ligands and for each filter-specific removal group.

        :return: None
        :rtype: None
        """
        # Create directory structure
        self.xyz_outdir = Path(self.outdir, 'concat_xyz')  # directory for concatenated xyz files
        self.xyz_outdir.mkdir(parents=True, exist_ok=True)

        # Save stdout output of filtering to info directory
        with open(Path(self.outdir, "filters.txt"), 'w') as f:
            f.write(self._get_filter_tracking_string())

        # Save a csv with an overview of all ligands to info directory
        self.df_all_ligands.to_csv(Path(self.outdir, "ligands_overview.csv"), index=True)

        # Save concatenated xyz files
        modes = ['Passed'] + [filter['unique_filtername'] for filter in self.filter_tracking]
        for mode in modes:
            # Get ligand IDs that were filtered out with this filter or passed
            filtered_ligand_ids = self.df_all_ligands.index[self.df_all_ligands['filter'] == mode]

            # Remove spaces from mode name so that a file has never a space in its name. Also remove ":" because it is not allowed in filenames on Windows.
            xyz_filename = f"concat_{mode.replace(' ', '_').replace(':', '')}.xyz"  #
            xyz_filepath = Path(self.xyz_outdir, xyz_filename)

            # Write concatenated xyz file
            if len(filtered_ligand_ids) > 0:
                with open(xyz_filepath, 'w') as f:
                    xyz_string = '\n'.join([self.db.db[ligand_id].get_xyz_string(comment=None, with_metal=self.output_metal) for ligand_id in filtered_ligand_ids])
                    f.write(xyz_string)

        return

    def _all_ligands_left(self) -> Generator[object, None, None]:
        """
        Yield ligand objects that passed the filters.

        Iterates over the list self.unames (unique ligand identifiers that remain after filtering)
        and yields the corresponding ligand objects from the in-memory ligand database.

        :return: An iterator over ligand objects (in the same order as self.unames).
        :rtype: Iterator[object]
        """
        for uname in self.unames:
            yield self.db.db[uname]

    def _get_ligand_df(self):
        """
        Construct a pandas DataFrame with CSV-relevant metadata for the current ligand subset.

        The method queries each ligand via get_csv_info(max_entries=5) and assembles a DataFrame
        with ligand metadata suitable for writing to CSV or summary tables.

        :return: DataFrame where each row corresponds to a ligand (index are unique names until set elsewhere).
        :rtype: pd.DataFrame
        """
        ligands = {uname: self.db.db[uname].get_csv_info(max_entries=5) for uname in self.db.db.keys()}
        df = pd.DataFrame.from_dict(ligands, orient='index')
        return df

    def _save_ligand_info_csv(self):
        """
        Save a CSV file named 'ligand_info.csv' next to the output ligand database containing ligand metadata.

        The CSV is written to the parent directory of self.output_ligand_db_path with filename 'ligand_info.csv'.
        The DataFrame includes the same columns produced by _get_ligand_df.

        :return: None
        :rtype: None
        """
        self.df_ligand_info = self._get_ligand_df()
        outpath = Path(self.output_ligand_db_path.parent, "ligand_info.csv")
        self.df_ligand_info.to_csv(outpath, index=False)

        return

    @classmethod
    def run_from_yaml(cls, input: Union[str,Path,None], n: Optional[int]=None) -> "LigandFilters":
        """
        Create and run a LigandFilters instance from a YAML specification file.

        If input is None, a default template ligandfilters.yml is used. The YAML file
        must contain top-level keys 'db' (optional), 'n' (optional) and the filter list and options
        required by the run(...) method.

        :param input: Path to the filter input file (.yml) or None to use a default template.
        :type input: Union[str, Path, None]
        :param n: Number of ligand objects to include in the output. Takes precedence over the 'n' value in the YAML file if provided.
        :type n: Union[int, None]
        :return: A LigandFilters instance after executing the configured filters.
        :rtype: LigandFilters
        """
        if input is None:
            input = default_ligandfilters_yml_path

        input_dict = read_yaml(input)
        input_db_file = input_dict.pop('db', None)
        n_yaml = input_dict.pop('n', None)
        if n is None:
            n = n_yaml

        filter = LigandFilters(db=input_db_file, n=n)
        filter.run(**input_dict)

        return filter

    @classmethod
    def run_from_cli(cls, input: Union[str, Path, None] = None, n: Optional[int] = None) -> "LigandFilters":
        """
        Run ligand filtering using command-line style setup helpers and a YAML input.

        This method wraps run_from_yaml with BaseModule CLI pre/post hooks to provide
        standardized CLI logging and argument printing.

        :param input: Path to the filter input file (.yml) or None to use the default template.
        :type input: Union[str, Path, None]
        :param n: Number of ligand objects to include in the output. Takes precedence over the 'n' value in the YAML file if provided.
        :type n: Union[int, None]
        :return: A LigandFilters instance after executing the configured filters.
        :rtype: LigandFilters
        """
        super()._before_run_from_cli()
        super()._print_cli_input(input=input)
        ligandfilters = cls.run_from_yaml(input=input, n=n)
        super()._after_run_from_cli()

        return ligandfilters

    def run(self,
            filters: list[dict],
            outpath: Union[str, Path] = 'filtered_ligand_db.jsonlines',
            dbinfo: bool = True,
            metal: bool = True,
            ) -> LigandDB:
        """
        Apply provided filters, save the filtered ligand database, and optionally save auxiliary info.

        .. tip::

            All the parameters below are available as well via the assembler .yml file as batch options (i.e. indented in the ``batches:`` list).

        The method executes the filtering pipeline, writes the filtered LigandDB to outpath (if provided),
        and optionally writes human-readable information (CSV, XYZ concatenations) controlled by dbinfo.
        It returns the list of unique ligand identifiers that passed all filters.

        :param list[dict] filters: List of filter specification dictionaries to apply in sequence.
        :param [None|str] outpath: Path to the output ligand database file. If None, no ligand DB file is written.
        :param bool dbinfo: If True, write additional info files (CSV, concatenated XYZ) to an info directory.
        :param bool metal: If True, in the concatenated XYZ files, include a pseudo metal center in the ligand structure for visualization.
        :return: LigandDB object containing only the ligands that passed all filters.
        :rtype: LigandDB
        """
        self.output_info = dbinfo
        self.output_metal = metal
        self.output_ligand_db_path = Path(outpath)
        outdirname = f'info_{self.output_ligand_db_path.with_suffix("").name}'
        self.outdir = Path(self.output_ligand_db_path.parent, outdirname)  # directory for full output info

        print(f"Starting DART Ligand Filters Module.")
        print(f"Input ligand db file: `{self.input_ligand_db_path.name}`")
        print(f"Output ligand db file: `{self.output_ligand_db_path.name}`")
        filtered_db = self.apply_filters(filters=filters)
        self.unames = list(filtered_db.db.keys())

        # Make a directory for the output if specified
        self.output_ligand_db_path.parent.mkdir(parents=True, exist_ok=True)
        filtered_db._to_json(self.output_ligand_db_path, json_lines=True,
                             desc=f'Save ligand db to `{self.output_ligand_db_path.name}`')

        if dbinfo:
            self.save_filtered_ligands_output()

        self.output = self._get_filter_tracking_string()
        print(self.output)

        return filtered_db



if __name__ == '__main__':
    filters = LigandFilters(db='metalig')
    #%% Filters from FeRedox
    db = filters.run(
        filters=[
            {'filter': 'property', 'name': 'n_denticities', 'values': [2, 3, 4]},
            {'filter': 'composition', 'elements': 'NO', 'instruction': 'must_only_contain_in_any_amount', 'only_donors': True},
            {'filter': 'composition', 'elements': 'CHNOPS', 'instruction': 'must_only_contain_in_any_amount', 'only_donors': False},
            # {'filter': 'smarts', 'smarts': 'S(~O)(~O)~O', 'should_contain': True, 'include_metal': False},
            # {'filter': 'smarts', 'smarts': 'S(=O)(=O)O', 'should_contain': True, 'include_metal': False},
            # {'filter': 'smarts', 'smarts': 'P(~O)(~O)~O', 'should_contain': True, 'include_metal': False},
            # {'filter': 'smarts', 'smarts': 'P(=O)(O)O', 'should_contain': True, 'include_metal': False},
            {'filter': 'smarts', 'smarts': '[C;D3](~[O;D1;H0,H1])(~[O;D1;H0,H1])~[!#1]', 'should_contain': True, 'include_metal': False},
            {'filter': 'property', 'name': 'n_haptic_groups', 'values': [0]},
            {'filter': 'parents', 'metal_centers': ['Fe']},
            {'filter': 'smarts', 'smarts': '[$([CX3H1](=[OX1])[OX2H1,O-]),$([CX3]([#6])(=[OX1])[OX2H1,O-])]',
             'should_contain': True, 'include_metal': False},
            {'filter': 'smarts', 'smarts': '[$([CX3H1](=[OX1])[OX2H1,O-]),$([CX3]([#6])(=[OX1])[OX2H1,O-])]',
             'should_contain': True, 'include_metal': True},

        ],
        outpath='/Users/timosommer/Downloads/test_DART/test_modules/test2/feredox_ligands.jsonlines',
    )