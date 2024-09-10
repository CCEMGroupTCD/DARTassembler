import shutil
import sys

from tqdm import tqdm
from DARTassembler.src.ligand_extraction.DataBase import LigandDB
from DARTassembler.src.constants.Paths import project_path
from DARTassembler.src.ligand_filters.FilteringStage import FilterStage
from pathlib import Path
import pandas as pd
pd.options.mode.chained_assignment = None   # silence pandas SettingWithCopyWarning
from typing import Union
from DARTassembler.src.assembly.Assembly_Input import LigandFilterInput, _mw, _filter, _ligand_charges, _ligcomp, _coords, \
    _metals_of_interest, _denticities_of_interest, _remove_ligands_with_neighboring_coordinating_atoms, \
    _remove_ligands_with_beta_hydrogens, _strict_box_filter, _acount, _acount_min, _acount_max, _denticities, \
    _ligcomp_atoms_of_interest, _ligcomp_instruction, _mw_min, _mw_max, _graph_hash_wm, _stoichiometry, _min, _max, \
    _interatomic_distances, _occurrences, _planarity, _remove_missing_bond_orders, _atm_neighbors, \
    _atom, _neighbors, _smarts_filter, _smarts, _should_be_present, _include_metal





class LigandFilters(object):

    def __init__(self, filepath: Union[str, Path], max_number: Union[int, None] = None, output_ligand_db_path: Union[None, str, Path] = None, delete_output_dir: bool = False):
        self.filepath = filepath
        self.max_number = max_number
        self.input = LigandFilterInput(path=self.filepath)

        self.ligand_db_path = self.input.ligand_db_path
        self.output_ligand_db_path = output_ligand_db_path or self.input.output_ligand_db_path  # Overwrite the output path if specified
        self.output_info = self.input.output_filtered_ligands
        self.filters = self.input.filters

        if delete_output_dir:
            shutil.rmtree(self.output_ligand_db_path.parent, ignore_errors=True)

        # Make a directory for the output if specified
        if self.output_info:
            outdirname = f'info_{self.output_ligand_db_path.with_suffix("").name}'
            self.outdir = Path(self.output_ligand_db_path.parent, outdirname)   # directory for full output
            self.outdir.mkdir(parents=True, exist_ok=True)
            self.xyz_outdir = Path(self.outdir, 'concat_xyz')                   # directory for concatenated xyz files
            if self.xyz_outdir.exists():    # Delete the existing directory so that there are no old files left that go into the new directory
                shutil.rmtree(self.xyz_outdir)
            self.xyz_outdir.mkdir(parents=True, exist_ok=True)

        self.filter_tracking = []

    def get_filtered_db(self) -> LigandDB:
        print(f"Starting DART Ligand Filters Module.")
        print(f"Input ligand db file: `{self.ligand_db_path.name}`")
        print(f"Output ligand db file: `{self.output_ligand_db_path.name}`")

        db = LigandDB.load_from_json(
            self.ligand_db_path,
            n_max=self.max_number,
        )

        self.Filter = FilterStage(db)
        self.n_ligands_before = len(self.Filter.database.db)
        self.df_all_ligands = self.get_ligand_df()
        self.df_all_ligands['Filter'] = None    # initialize column for filter tracking
        if self.output_info:
            self.all_xyz_strings = {ligand_id: ligand.get_xyz_file_format_string(comment=None, with_metal=True) for ligand_id, ligand in self.Filter.database.db.items()}   # save xyz strings for outputting later and keep a placeholder comment for replacing later

        # mandatory filters
        self.Filter.filter_charge_confidence(filter_for="confident")
        self.Filter.filter_unconnected_ligands()

        for idx, filter in tqdm(enumerate(self.filters), desc="Applying filters", unit=" filters", total=len(self.filters), file=sys.stdout):
            filtername = filter[_filter]
            unique_filtername = f"Filter {idx+1:02d}: {filtername}"    # name for printing filters for the user
            n_ligands_before = len(self.Filter.database.db)

            if filtername == _denticities_of_interest:
                self.Filter.denticity_of_interest_filter(denticity_of_interest=filter[_denticities_of_interest])

            elif filtername == _graph_hash_wm:
                self.Filter.graph_hash_with_metal_filter(graph_hashes_with_metal=filter[_graph_hash_wm])

            elif filtername == _remove_ligands_with_neighboring_coordinating_atoms:
                if filter[_remove_ligands_with_neighboring_coordinating_atoms]:
                    self.Filter.filter_neighbouring_coordinating_atoms(denticities=filter[_denticities])

            elif filtername == _remove_ligands_with_beta_hydrogens:
                if filter[_remove_ligands_with_beta_hydrogens]:
                    self.Filter.filter_betaHs(denticities=filter[_denticities])

            # ====== Denticity dependent filters ======
            elif filtername == _acount:
                self.Filter.filter_atom_count(min=filter[_acount_min], max=filter[_acount_max], denticities=filter[_denticities])

            elif filtername == _ligcomp:
                self.Filter.filter_ligand_atoms(
                    denticity=filter[_denticities],
                    atoms_of_interest=filter[_ligcomp_atoms_of_interest],
                    instruction=filter[_ligcomp_instruction])

            elif filtername == _metals_of_interest:
                self.Filter.metals_of_interest_filter(
                    denticity=filter[_denticities],
                    metals_of_interest=filter[_metals_of_interest])

            elif filtername == _ligand_charges:
                self.Filter.filter_ligand_charges(
                    denticity=filter[_denticities],
                    charge=filter[_ligand_charges])

            elif filtername == _coords:
                self.Filter.filter_coordinating_group_atoms(
                    denticity=filter[_denticities],
                    atoms_of_interest=filter[_ligcomp_atoms_of_interest],
                    instruction=filter[_ligcomp_instruction])

            elif filtername == _mw:
                self.Filter.filter_molecular_weight(
                    min=filter[_mw_min],
                    max=filter[_mw_max],
                    denticities=filter[_denticities]
                )
            elif filtername == _occurrences:
                self.Filter.filter_occurrences(
                    min=filter[_min],
                    max=filter[_max],
                    denticities=filter[_denticities]
                )
            elif filtername == _interatomic_distances:
                self.Filter.filter_interatomic_distances(
                    min=filter[_min],
                    max=filter[_max],
                    denticities=filter[_denticities]
                )
            elif filtername == _planarity:
                self.Filter.filter_planarity(
                    min=filter[_min],
                    max=filter[_max],
                    denticities=filter[_denticities]
                )
            elif filtername == _stoichiometry:  # deprecated
                self.Filter.stoichiometry_filter(
                    stoichiometry=filter[_stoichiometry],
                    denticities=filter[_denticities]
                )
            elif filtername == _remove_missing_bond_orders:
                self.Filter.filter_missing_bond_orders(
                    denticities=filter[_denticities]
                )
            elif filtername == _atm_neighbors:
                self.Filter.filter_atomic_neighbors(
                    atom=filter[_atom],
                    neighbors=filter[_neighbors],
                    denticities=filter[_denticities]
                )
            elif filtername == _smarts_filter:
                self.Filter.filter_smarts_substructure_search(
                    smarts=filter[_smarts],
                    should_be_present=filter[_should_be_present],
                    include_metal=filter[_include_metal],
                    denticities=filter[_denticities]
                )
            else:
                raise ValueError(f"Unknown filter: {filtername}")

            # To the dataframe with all ligands, add a column specifying which filter was applied to filter this ligand. This is important for outputting a csv with all filtered out ligands later.
            ligand_was_filtered = ~self.df_all_ligands.index.isin(self.Filter.database.db.keys()) & (self.df_all_ligands['Filter'].isna())
            self.df_all_ligands.loc[ligand_was_filtered, 'Filter'] = unique_filtername

            n_ligands_after = len(self.Filter.database.db)
            self.filter_tracking.append({
                "filter": filtername,
                "unique_filtername": unique_filtername,
                "n_ligands_before": n_ligands_before,
                "n_ligands_after": n_ligands_after,
                "n_ligands_removed": n_ligands_before - n_ligands_after,
                "full_filter_options": {name: option for name, option in filter.items() if name != _filter}
            })

        self.n_ligands_after = len(self.Filter.database.db)

        # Clean up the ligand df
        self.df_all_ligands.fillna({'Filter': 'Passed'}, inplace=True)      # fill in 'Passed' for ligands that were not filtered out
        self.df_all_ligands.set_index('Ligand ID', inplace=True)    # set index to ligand ID, making sure that the column in the csv is named 'Ligand ID'
        columns = ['Filter'] + [col for col in self.df_all_ligands.columns if col != 'Filter']
        self.df_all_ligands = self.df_all_ligands[columns]                # move 'Filter' column to the front
        self.df_all_ligands = self.df_all_ligands.sort_values(by='Filter')# sort by filter name

        return self.Filter.database

    def get_filter_tracking_string(self) -> str:
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
        denticity_count = pd.Series(
            [lig.denticity for lig in self.Filter.database.db.values()]).value_counts().to_dict()
        dent_output = ', '.join(sorted([f'{dent}: {count}' for dent, count in denticity_count.items()]))

        output += f"{'  Total summary of DART Ligand Filters run  ':=^80}\n"
        output += f"Before filtering:  {self.n_ligands_before} ligands\n"
        output += f"Filtered out:      {self.n_ligands_before - self.n_ligands_after} ligands\n"
        output += f"Passed:            {self.n_ligands_after} ligands\n"
        output += f"Denticities:       {dent_output}\n"

        # If the number of ligands after filtering is small, print them explicitly
        if self.n_ligands_after <= 10:
            stoichiometries = ','.join([ligand.stoichiometry for ligand in self.Filter.database.db.values()])
            output += f'Passed ligands:    {stoichiometries}\n'

        output += f"Filtered ligand database with {self.n_ligands_after} ligands was saved to `{self.output_ligand_db_path.name}`.\n"
        if self.output_info:
            output += f"Info on filtered ligands saved to directory `{self.outdir.name}`.\n"
        output += "Done! All ligands filtered. Exiting DART Ligand Filters Module."

        return output

    def save_filtered_ligand_db(self):
        filtered_db = self.get_filtered_db()

        if not self.output_ligand_db_path.parent.exists():
            self.output_ligand_db_path.parent.mkdir(parents=True)
        filtered_db.to_json(self.output_ligand_db_path, json_lines=True, desc=f'Save ligand db to `{self.output_ligand_db_path.name}`')

        self.output = self.get_filter_tracking_string()

        # Optionally output filtered ligands info
        if self.output_info:
            self.save_filtered_ligands_output()

        print(self.output)

        return

    def save_filtered_ligands_output(self) -> None:
        """
        Saves a directory with an overview of all ligands that were filtered out and passed. This overview contains both a csv file with all ligands and one concatenated xyz file for each filter plus the passed ligands.
        """

        # Save stdout output of filtering to info directory
        with open(Path(self.outdir, "filters.txt"), 'w') as f:
            f.write(self.get_filter_tracking_string())

        # Save a csv with an overview of all ligands to info directory
        self.df_all_ligands.to_csv(Path(self.outdir, "ligands_overview.csv"), index=True)

        # Save concatenated xyz files
        modes = ['Passed'] + [filter['unique_filtername'] for filter in self.filter_tracking]
        for mode in modes:
            # Get ligand IDs that were filtered out with this filter or passed
            filtered_ligand_ids = self.df_all_ligands.index[self.df_all_ligands['Filter'] == mode]

            # Remove spaces from mode name so that a file has never a space in its name
            xyz_filename = f"concat_{mode.replace(' ', '')}.xyz"
            xyz_filepath = Path(self.xyz_outdir, xyz_filename)

            # Write concatenated xyz file
            with open(xyz_filepath, 'w') as f:
                for ligand_id in filtered_ligand_ids:
                    xyz_string = self.all_xyz_strings[ligand_id]
                    f.write(xyz_string)

        return

    def get_ligand_df(self):
        db = self.Filter.database
        ligands = {uname: ligand.get_ligand_output_info(max_entries=5) for uname, ligand in db.db.items()}
        return pd.DataFrame.from_dict(ligands, orient='index')

    def save_ligand_info_csv(self):
        self.df_ligand_info = self.get_ligand_df()
        outpath = Path(self.output_ligand_db_path.parent, "ligand_info.csv")
        self.df_ligand_info.to_csv(outpath, index=False)

        return



if __name__ == "__main__":
    ligand_filter_path = project_path().extend(*'testing/integration_tests/ligand_filters/data_input/ligandfilters.yml'.split('/'))
    max_number = 1000


    filter = LigandFilters(filepath=ligand_filter_path, max_number=max_number)
    filter.save_filtered_ligand_db()

