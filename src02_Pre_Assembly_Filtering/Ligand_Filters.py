from src01.DataBase import LigandDB
from constants.Paths import project_path, default_ligand_db_path
from src02_Pre_Assembly_Filtering.FilteringStage import FilterStage
from test.Integration_Test import IntegrationTest
from pathlib import Path
from typing import Union
from src05_Assembly_Refactor.Assembly_Input import LigandFilterInput, _mw, _filter, _ligand_charges, _ligcomp, _coords, \
    _metals_of_interest, _denticities_of_interest, _remove_ligands_with_neighboring_coordinating_atoms, \
    _remove_ligands_with_beta_hydrogens, _only_confident_charges, _strict_box_filter, _filter_even_odd_electron_count, \
    _acount, _acount_min, _acount_max, _denticities, _ligcomp_atoms_of_interest, _ligcomp_instruction, _mw_min, _mw_max
from src01.io_custom import load_unique_ligand_db

# todo: add mandatory filter which checks if a charge is present

class LigandFilters(object):

    def __init__(self, filepath: Union[str,Path], max_number: Union[int, None] = None):
        self.filepath = Path(filepath)
        self.max_number = max_number
        self.input = LigandFilterInput(path=self.filepath)

        self.ligand_db_path = self.input.ligand_db_path
        self.output_ligand_db_path = self.input.output_ligand_db_path
        self.filters = self.input.filters

        self.filter_tracking = []

    def get_filtered_db(self) -> LigandDB:
        db = LigandDB.from_json(
                                json_=self.ligand_db_path,
                                type_="Ligand",
                                max_number=self.max_number
                                )

        self.Filter = FilterStage(db)
        self.n_ligands_before = len(self.Filter.database.db)

        for filter in self.filters:
            filtername = filter[_filter]
            n_ligands_before = len(self.Filter.database.db)

            if filtername == _denticities_of_interest:
                self.Filter.denticity_of_interest_filter(denticity_of_interest=filter[_denticities_of_interest])

            elif filtername == _remove_ligands_with_neighboring_coordinating_atoms:
                if filter[_remove_ligands_with_neighboring_coordinating_atoms]:
                    self.Filter.filter_neighbouring_coordinating_atoms()

            elif filtername == _only_confident_charges:
                if filter[_only_confident_charges]:
                    self.Filter.filter_charge_confidence(filter_for="confident")

            elif filtername == _remove_ligands_with_beta_hydrogens:
                if filter[_remove_ligands_with_beta_hydrogens]:
                    self.Filter.filter_betaHs()

            elif filtername == _strict_box_filter:
                if filter[_strict_box_filter]:
                    self.Filter.box_excluder_filter()

            elif filtername == _filter_even_odd_electron_count:
                self.Filter.filter_even_odd_electron(filter_for=filter[_filter_even_odd_electron_count])

            # Denticity dependent filters
            elif filtername == _acount:
                min_atom_count = filter[_acount_min]
                max_atom_count = filter[_acount_max]
                denticities = filter[_denticities]
                if min_atom_count is not None:
                    self.Filter.filter_atom_count(denticity=denticities, number=min_atom_count, instruction="greater_than")
                if max_atom_count is not None:
                    self.Filter.filter_atom_count(denticity=denticities, number=max_atom_count, instruction="less_than")

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
                                                denticity=filter[_denticities],
                                                atomic_weight_min=filter[_mw_min],
                                                atomic_weight_max=filter[_mw_max]
                                                )

            n_ligands_after = len(self.Filter.database.db)
            self.filter_tracking.append({
                "filter": filtername,
                "n_ligands_before": n_ligands_before,
                "n_ligands_after": n_ligands_after,
                "n_ligands_removed": n_ligands_before - n_ligands_after,
                "full_filter_options": {name: option for name, option in filter.items() if name != _filter}
                })

            # todo: subgraph match
        self.n_ligands_after = len(self.Filter.database.db)

        return self.Filter.database

    def get_filter_tracking_string(self) -> str:
        output= "===================== FILTER TRACKING =====================\n"
        output += f"An overview of filters applied to '{self.ligand_db_path}'. The new ligand database is saved as '{self.output_ligand_db_path}'.\nRemoved ligands:\n"
        for filter in self.filter_tracking:
            output += f"  - {filter['filter']}: {filter['n_ligands_removed']}\n"
            options = ', '.join([f"{name}: {option}" for name, option in filter['full_filter_options'].items()])
            output += f"        options --> {options}\n"
        output += f"\nNumber of ligands before filtering: {self.n_ligands_before}\n"
        output += f"Number of ligands filtered out: {self.n_ligands_before - self.n_ligands_after}\n"
        output += f"Number of ligands after filtering: {self.n_ligands_after}\n"

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


if __name__ == "__main__":

    ligand_filter_path = project_path().extend('src05_Assembly_Refactor', 'ligandfilters.yml')
    max_number = None

    filter = LigandFilters(filepath=ligand_filter_path, max_number=max_number)
    filter.save_filtered_ligand_db()
    filter.save_filter_tracking()
    filter.print_filter_tracking()

    # Check if the new filtered db is the same as the old one
    benchmark_dir = project_path().extend("src14_Assembly_Unit_Test", 'ligandfilters_benchmark')
    if benchmark_dir.exists():
        test = IntegrationTest(new_dir=filter.output_ligand_db_path.parent, old_dir=benchmark_dir)
        test.compare_all()
    else:
        print("No benchmark directory found. Could not perform integration test.")