import shutil
from copy import deepcopy
from pathlib import Path
from typing import Union
import pandas as pd
from DARTassembler.src.ligand_extraction.DataBase import LigandDB
from DARTassembler.src.ligand_extraction.utilities_Molecule import get_standardized_stoichiometry_from_atoms_list
from DARTassembler.src.metalig.metalig_utils import get_correct_ligand_db_path_from_input

filters = [
    {
        'filter': 'property',
        'name': 'n_atoms',
        'range': (1,80)
    },
    {
        'filter': 'property',
        'name': 'charge',
        'values': [-2, -1, 0]
    },
    {
        'filter': 'parent_complexes',
        'metal_centers': ['Fe', 'Mn', 'Cr', 'Ru', 'Co', 'Ni', 'Cu', 'Ir', 'Mo'],
    },
    {
        'filter': 'composition',
        'elements': ['O'],
        'instruction': 'must_exclude',
        'only_donors': True,
    },
    {
        'filter': 'property',
        'name': 'molecular_weight',
        'range': (30.1, 300.2)
    },
    {
        'filter': 'property',
        'name': 'planarity',
        'range': (0.2, 1.0)
    },
    {
        'filter': 'property',
        'name': 'n_ligand_instances',
        'range': (5, 999999999)
    },
    {
        'filter': 'smarts',
        'smarts': 'N=N',
        'should_contain': False,
        'include_metal': True
    },
    {
        'filter': 'smarts',
        'smarts': '[C&H2]',
        'should_contain': False,
        'include_metal': True
    },
    {
        'filter': 'property',
        'name': 'n_beta_hydrogens',
        'range': (0, 0)
    }
]



oer_ligand_filters = [
    {
        'filter': 'property',
        'name': 'n_donors',
        'range': (2, 3)
    },
    {
        'filter': 'property',
        'name': 'n_haptic_atoms',
        'range': (0, 0)
    },
    {
        'filter': 'property',
        'name': 'n_beta_hydrogens',
        'range': (0, 0)
    },
    {
        'filter': 'parent_complexes',
        'metal_centers': ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Mo', 'W', 'Tc', 'Re', 'Ru', 'Os', 'Rh', 'Ir']
    },
    {
        'filter': 'composition',
        'elements': ['C', 'H', 'N', 'O', 'F'],
        'instruction': 'must_only_contain_in_any_amount',
        'only_donors': False
    },
    {
        'filter': 'composition',
        'elements': ['C', 'H'],
        'instruction': 'must_at_least_contain',
        'only_donors': False
    },
    {
        'filter': 'composition',
        'elements': ['N', 'O'],
        'instruction': 'must_only_contain_in_any_amount',
        'only_donors': True
    },
    {
        'filter': 'smarts',
        'smarts': '[C&H2]',
        'should_contain': False,
        'include_metal': True
    },
    {
        'filter': 'smarts',
        'smarts': '[N]~[N!R]',
        'should_contain': False,
        'include_metal': False
    },
]
# OER ligand filters in old format
# filters:
#   - filter: denticities
#     denticities: [2, 3]
#
#   - filter: remove_ligands_with_adjacent_coordinating_atoms
#     remove_ligands_with_adjacent_coordinating_atoms: True
#
#   - filter: remove_ligands_with_beta_hydrogens
#     remove_ligands_with_beta_hydrogens: True
#
#   - filter: metal_ligand_binding_history
#     metal_ligand_binding_history: [Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, Mo, W, Tc, Re, Ru, Os, Rh, Ir]  #[Cr, Mn, Fe, Co, Ru]
#     apply_to_denticities:
#
#   - filter: ligand_composition
#     elements: [C, H, N, O, F]
#     instruction: must_only_contain_in_any_amount
#     apply_to_denticities:
#
#   - filter: ligand_composition
#     elements: [C, H]
#     instruction: must_at_least_contain
#     apply_to_denticities:
#
#   - filter: coordinating_atoms_composition
#     elements: [N, O]
#     instruction: must_only_contain_in_any_amount
#     apply_to_denticities:
#
#   - filter: remove_ligands_with_missing_bond_orders
#     remove_ligands_with_missing_bond_orders: True
#     apply_to_denticities:
#
#   - filter: smarts  # remove CH2 units only, leave CH3
#     smarts: '[C&H2]'
#     should_contain: False
#     include_metal: True
#     apply_to_denticities:
#
#   - filter: smarts  # remove N~N bonds not in rings with arbitrary bond orders
#     smarts: '[N]~[N!R]'
#     should_contain: False
#     include_metal: False
#     apply_to_denticities:



# New filters, everything should be setup for these filters
possible_filters = {
    'property': ['name', 'range', 'values'],
    'composition': ['elements', 'instruction', 'only_donors'],
    'parent_complexes': ['metal_centers'],
    'smarts': ['smarts', 'should_contain', 'include_metal'],
}

# Refactor the LigandFilters class
class NewLigandFilters(object):

    def __init__(self, input_db_file: Union[str, Path], n_max_ligands: Union[int, None] = None):
        self.n_max_ligands = n_max_ligands
        self.input_ligand_db_path = get_correct_ligand_db_path_from_input(input_db_file)
        self.db = LigandDB.load_from_json(path=input_db_file, n_max=n_max_ligands)

    def _apply_filters(self, filters: list[dict]) -> list[str]:
        self.filter_tracking = []
        filters = deepcopy(filters)

        # Prepend a bond-order filter if a SMARTS filter is used to make sure the SMARTS filter can be applied.
        filters = self._add_filter_for_valid_bond_orders_if_smarts(filters)

        self.unames = list(self.db.db.keys())
        self.n_ligands_before = len(self.unames)
        self.df_all_ligands = self._get_ligand_df()
        self.df_all_ligands['filter'] = None    # initialize column for filter tracking

        for idx, filter in enumerate(filters):
            n_ligands_before = len(self.unames)
            filtername = filter.pop('filter')
            if filtername == 'property':
                self.unames = [uname for uname in self.unames if self.db.db[uname].has_global_property_in_range(**filter)]
                name_appendix = filter['name']
            elif filtername == 'composition':
                self.unames = [uname for uname in self.unames if self.db.db[uname].has_specified_stoichiometry(**filter)]
                name_appendix = get_standardized_stoichiometry_from_atoms_list(filter['elements'])
            elif filtername == 'parent_complexes':
                self.unames = [uname for uname in self.unames if self.db.db[uname].has_specified_metal_centers(**filter)]
                name_appendix = ', '.join(filter['metal_centers'])
            elif filtername == 'smarts':
                self.unames = [uname for uname in self.unames if self.db.db[uname].has_specified_smarts(**filter)]
                name_appendix = filter['smarts']
            else:
                raise ValueError(f'Filter "{filtername}" not recognized!')

            n_ligands_after = len(self.unames)
            ligand_was_filtered = ~self.df_all_ligands.index.isin(self.unames) & (self.df_all_ligands['filter'].isna())
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
        self.n_ligands_after = len(self.unames)

        self.df_all_ligands.fillna({'filter': 'Passed'}, inplace=True)      # fill in 'Passed' for ligands that were not filtered out
        self.df_all_ligands.set_index('unique_name', inplace=True)    # set index to ligand ID, making sure that the column in the csv is named 'Ligand ID'
        columns = ['filter'] + [col for col in self.df_all_ligands.columns if col != 'filter']
        self.df_all_ligands = self.df_all_ligands[columns]                # move 'filter' column to the front
        self.df_all_ligands = self.df_all_ligands.sort_values(by='filter')# sort by filter name

        return self.unames

    def _add_filter_for_valid_bond_orders_if_smarts(self, filters: list[dict]):
        """
        In case of a SMARTS filter, ligands with missing bond orders are removed. This filter adds a filter that removes ligands with missing bond orders.
        :return: updated list of filters
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
        geometry_count = pd.Series([lig.geometry for lig in self._all_ligands_left()]).value_counts().to_dict()
        geometry_output = ', '.join(sorted([f'{geom} ({count})' for geom, count in geometry_count.items()]))

        n_ligands_before = len(self.df_all_ligands)
        n_ligands_after = len(self.unames)

        output += f"{'  Total summary of DART Ligand Filters run  ':=^80}\n"
        output += f"Before filtering:           {n_ligands_before} ligands\n"
        output += f"Filtered out:               {n_ligands_before - n_ligands_after} ligands\n"
        output += f"Passed:                     {n_ligands_after} ligands\n"
        output += f"Passed geometries:          {geometry_output}\n"

        # Print the stoichiometries of the first five passed ligands.
        stoichiometries = ', '.join([lig.stoichiometry for idx, lig in enumerate(self._all_ligands_left()) if idx < 5])
        ellipsis = ', ...' if n_ligands_after > 5 else ''
        output += f'Passed stoichiometries:     {stoichiometries}{ellipsis}\n'

        output += f"Filtered ligand database with {n_ligands_after} ligands was saved to `{self.output_ligand_db_path.name}`.\n"
        if self.output_info:
            output += f"Info on filtered ligands saved to directory `{self.outdir.name}`.\n"
        output += "Done! All ligands filtered. Exiting DART Ligand Filters Module."

        return output

    def save_filtered_ligands_output(self) -> None:
        """
        Saves a directory with an overview of all ligands that were filtered out and passed. This overview contains both a csv file with all ligands and one concatenated xyz file for each filter plus the passed ligands.
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
                    for ligand_id in filtered_ligand_ids:
                        xyz_string = self.db.db[ligand_id].get_xyz_file_format_string(comment=None, with_metal=self.output_metal)
                        f.write(xyz_string)

        return

    def _all_ligands_left(self):
        for uname in self.unames:
            yield self.db.db[uname]

    def _get_ligand_df(self):
        ligands = {uname: self.db.db[uname].get_ligand_output_info(max_entries=5) for uname in self.unames}
        df = pd.DataFrame.from_dict(ligands, orient='index')
        return df

    def _save_ligand_info_csv(self):
        self.df_ligand_info = self._get_ligand_df()
        outpath = Path(self.output_ligand_db_path.parent, "ligand_info.csv")
        self.df_ligand_info.to_csv(outpath, index=False)

        return

    def get_filtered_db(self,
                        filters: list[dict],
                        output_db_file: Union[None, str, Path] = 'filtered_ligand_db.jsonlines',
                        output_ligands_info: bool = True,
                        output_metal: bool = True,
                        pre_delete: bool = False
                        ):
        """
        Apply filters to the ligand database and save the filtered ligands to a new ligand database file.
        :param filters: list of dictionaries with filter options.
        :param output_db_file: path to the output ligand database file. If None, no output files are saved.
        :param output_ligands_info: if a directory with info on the filtered ligands should be saved additionally to the ligand database.
        :param output_metal: if a pseudo metal center should be added to the ligands in the concatenated xyz files.
        :param pre_delete: if the output ligand db file and the info directory should be deleted before the new files are saved.
        :return: list of unique names of the ligands that passed all filters.
        """
        self.output_ligand_db_path = Path(output_db_file)
        self.output_info = output_ligands_info
        self.output_metal = output_metal
        outdirname = f'info_{self.output_ligand_db_path.with_suffix("").name}'
        self.outdir = Path(self.output_ligand_db_path.parent, outdirname)  # directory for full output info

        if pre_delete:
            self.output_ligand_db_path.unlink(missing_ok=True)
            shutil.rmtree(self.outdir, ignore_errors=True)

        print(f"Starting DART Ligand Filters Module.")
        print(f"Input ligand db file: `{self.input_ligand_db_path.name}`")
        print(f"Output ligand db file: `{self.output_ligand_db_path.name}`")
        self._apply_filters(filters=filters)

        # Make a directory for the output if specified
        if self.output_ligand_db_path is not None:
            filtered_db = LigandDB({uname: self.db.db[uname] for uname in self.unames})
            self.output_ligand_db_path.parent.mkdir(parents=True, exist_ok=True)
            filtered_db.to_json(self.output_ligand_db_path, json_lines=True,
                                desc=f'Save ligand db to `{self.output_ligand_db_path.name}`')

            if output_ligands_info:
                self.save_filtered_ligands_output()

        self.output = self._get_filter_tracking_string()
        print(self.output)

        return self.unames



if __name__ == '__main__':
    n_max = None
    outpath = 'data/oer_ligandfilters/filtered_ligands.jsonlines'
    output_info = True

    # #%% ==============    Refactor ligand filters    ==================
    filter = NewLigandFilters(input_db_file='metalig', n_max=n_max)
    ligands = filter.get_filtered_db(oer_ligand_filters, output_db_file=outpath, output_ligands_info=output_info)

    # #%% ==============    Doublecheck refactoring    ==================
    # from dev.test.Integration_Test import IntegrationTest
    # old_dir = Path('/Users/timosommer/PhD/projects/DARTassembler/dev/DART_refactoring_to_v1_1_0/data/benchmark_ligandfilters')
    # if old_dir.exists():
    #     test = IntegrationTest(new_dir=Path(outpath).parent, old_dir=old_dir)
    #     test.compare_all()
    #     print('Test for ligand filters passed!')
    # else:
    #     print(f'ATTENTION: could not find benchmark folder "{old_dir}"!')