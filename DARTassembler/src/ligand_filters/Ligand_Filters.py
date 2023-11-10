from DARTassembler.src.ligand_extraction.DataBase import LigandDB
from DARTassembler.src.constants.Paths import project_path
from DARTassembler.src.ligand_filters.FilteringStage import FilterStage
from dev.test.Integration_Test import IntegrationTest
from pathlib import Path
import pandas as pd
from typing import Union
from DARTassembler.src.assembly.Assembly_Input import LigandFilterInput, _mw, _filter, _ligand_charges, _ligcomp, _coords, \
    _metals_of_interest, _denticities_of_interest, _remove_ligands_with_neighboring_coordinating_atoms, \
    _remove_ligands_with_beta_hydrogens, _strict_box_filter, _acount, _acount_min, _acount_max, _denticities, \
    _ligcomp_atoms_of_interest, _ligcomp_instruction, _mw_min, _mw_max, _graph_hash_wm, _stoichiometry, _min, _max, \
    _md_bond_length, _interatomic_distances, _occurrences, _planarity





class LigandFilters(object):

    def __init__(self, filepath: Union[str, Path], max_number: Union[int, None] = None):
        self.filepath = filepath
        self.max_number = max_number
        self.input = LigandFilterInput(path=self.filepath)

        self.ligand_db_path = self.input.ligand_db_path
        self.output_ligand_db_path = self.input.output_ligand_db_path
        self.filters = self.input.filters

        self.filter_tracking = []

    def get_filtered_db(self) -> LigandDB:
        print(f"Filtering ligand database: {self.ligand_db_path} --> {self.output_ligand_db_path}...")
        db = LigandDB.load_from_json(
            self.ligand_db_path,
            n_max=self.max_number,
        )

        self.Filter = FilterStage(db)
        self.n_ligands_before = len(self.Filter.database.db)

        # mandatory filters
        self.Filter.filter_charge_confidence(filter_for="confident")
        self.Filter.filter_unconnected_ligands()

        for filter in self.filters:
            filtername = filter[_filter]
            n_ligands_before = len(self.Filter.database.db)

            if filtername == _stoichiometry:
                self.Filter.stoichiometry_filter(stoichiometry=filter[_stoichiometry], denticities=filter[_denticities])

            if filtername == _denticities_of_interest:
                self.Filter.denticity_of_interest_filter(denticity_of_interest=filter[_denticities_of_interest])

            elif filtername == _graph_hash_wm:
                self.Filter.graph_hash_with_metal_filter(graph_hashes_with_metal=filter[_graph_hash_wm])

            elif filtername == _remove_ligands_with_neighboring_coordinating_atoms:
                if filter[_remove_ligands_with_neighboring_coordinating_atoms]:
                    self.Filter.filter_neighbouring_coordinating_atoms()

            elif filtername == _remove_ligands_with_beta_hydrogens:
                if filter[_remove_ligands_with_beta_hydrogens]:
                    self.Filter.filter_betaHs()

            elif filtername == _strict_box_filter:
                if filter[_strict_box_filter]:
                    self.Filter.box_excluder_filter()

            # Denticity dependent filters
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
            elif filtername == _md_bond_length:
                self.Filter.filter_metal_donor_bond_lengths(
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
            elif filtername == _stoichiometry:
                self.Filter.stoichiometry_filter(
                    stoichiometry=filter[_stoichiometry],
                    denticities=filter[_denticities]
                )
            else:
                raise ValueError(f"Unknown filter: {filtername}")


            n_ligands_after = len(self.Filter.database.db)
            self.filter_tracking.append({
                "filter": filtername,
                "n_ligands_before": n_ligands_before,
                "n_ligands_after": n_ligands_after,
                "n_ligands_removed": n_ligands_before - n_ligands_after,
                "full_filter_options": {name: option for name, option in filter.items() if name != _filter}
            })

            # todo: subgraph matching filter

        self.n_ligands_after = len(self.Filter.database.db)

        return self.Filter.database


    def get_filter_tracking_string(self) -> str:
        output= "===================== FILTER TRACKING =====================\n"
        output += f"--> Filtered ligand database was saved as '{self.output_ligand_db_path}'.\nRemoved ligands:\n"


        for filter in self.filter_tracking:
            output += f"  - {filter['filter']}: {filter['n_ligands_removed']}\n"
            options = ', '.join([f"{name}: {option}" for name, option in filter['full_filter_options'].items()])
            output += f"        options --> {options}\n"
        output += f"\nNumber of ligands before filtering: {self.n_ligands_before}\n"
        output += f"Number of ligands filtered out: {self.n_ligands_before - self.n_ligands_after}\n"
        output += f"Number of ligands after filtering: {self.n_ligands_after}\n"

        denticity_count = pd.Series([lig.denticity for lig in self.Filter.database.db.values()]).value_counts().to_dict()
        dent_output = ', '.join(sorted([f'{dent}: {count}' for dent, count in denticity_count.items()]))
        output += f"Number of ligands per denticity: {dent_output}\n"

        if self.n_ligands_after <= 10:
            stoichiometries = ','.join([ligand.stoichiometry for ligand in self.Filter.database.db.values()])
            output += f'  --> Remaining ligands: {stoichiometries}\n'

        return output

    def print_filter_tracking(self):
        print(self.get_filter_tracking_string())
        return

    def save_filter_tracking(self):
        outpath = Path(self.output_ligand_db_path.parent, "filter_tracking.txt")
        with open(outpath, 'w') as f:
            f.write(self.get_filter_tracking_string())

        return

    def save_filtered_ligand_db(self):
        filtered_db = self.get_filtered_db()

        if not self.output_ligand_db_path.parent.exists():
            self.output_ligand_db_path.parent.mkdir(parents=True)
        filtered_db.to_json(self.output_ligand_db_path, json_lines=True)

        return

    def save_ligand_info_csv(self):
        db = self.Filter.database
        ligands = {uname: ligand.get_ligand_output_info(max_entries=5) for uname, ligand in db.db.items()}
        self.df_ligand_info = pd.DataFrame.from_dict(ligands, orient='index')
        outpath = Path(self.output_ligand_db_path.parent, "ligand_info.csv")
        self.df_ligand_info.to_csv(outpath, index=False)

        return



if __name__ == "__main__":
    ligand_filter_path = project_path().extend(*'dev/test/2-1-1-Cl-example/ligandfilters.yml'.split('/'))
    max_number = 100


    filter = LigandFilters(filepath=ligand_filter_path, max_number=max_number)
    filter.save_filtered_ligand_db()
    filter.save_filter_tracking()
    filter.print_filter_tracking()

