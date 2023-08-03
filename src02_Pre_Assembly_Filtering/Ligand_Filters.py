from src01.DataBase import LigandDB
from constants.Paths import project_path, default_ligand_db_path
from src02_Pre_Assembly_Filtering.FilteringStage import FilterStage
from test.Integration_Test import IntegrationTest
from pathlib import Path
from typing import Union
from src05_Assembly_Refactor.Assembly_Input import LigandFilterInput
from src01.io_custom import load_unique_ligand_db

class LigandFilters(object):

    def __init__(self, filepath: Union[str,Path], max_number: Union[int, None] = None):
        self.filepath = Path(filepath)
        self.max_number = max_number
        self.input = LigandFilterInput(path=self.filepath)

        self.ligand_db_path = self.input.ligand_db_path if self.input.ligand_db_path is not None else default_ligand_db_path
        self.output_ligand_db_path = self.input.output_ligand_db_path if self.input.output_ligand_db_path is not None else self.get_output_ligand_db_path()
        self.denticities_of_interest = self.input.denticities_of_interest
        self.remove_ligands_with_neighboring_coordinating_atoms = self.input.remove_ligands_with_neighboring_coordinating_atoms if self.input.remove_ligands_with_neighboring_coordinating_atoms is not None else False
        self.only_confident_charges = self.input.only_confident_charges if self.input.only_confident_charges is not None else False
        self.remove_ligands_with_beta_hydrogens = self.input.remove_ligands_with_beta_hydrogens if self.input.remove_ligands_with_beta_hydrogens is not None else False
        self.strict_box_filter = self.input.strict_box_filter if self.input.strict_box_filter is not None else False
        self.filter_even_odd_electron_count = self.input.filter_even_odd_electron_count
        self.denticity_dependent_filters = self.input.dentfilters


    def get_output_ligand_db_path(self):
        """
        Returns a path to the output ligand database. If the path already exists, it will be renamed to avoid overwriting.
        """
        path = Path('filtered_ligand_db.json').resolve()

        idx = 1
        while path.exists():
            path = Path(f'{path}_{idx}')
            idx += 1

        return path

    def get_filtered_db(self) -> LigandDB:
        db = LigandDB.from_json(
                                json_=self.ligand_db_path,
                                type_="Ligand",
                                max_number=self.max_number
                                )

        Filter = FilterStage(db)
        if self.denticities_of_interest is not None:
            Filter.denticity_of_interest_filter(denticity_of_interest=self.denticities_of_interest)
        if self.remove_ligands_with_neighboring_coordinating_atoms:
            Filter.filter_neighbouring_coordinating_atoms()
        if self.only_confident_charges:
            Filter.filter_charge_confidence(filter_for="confident")
        if self.remove_ligands_with_beta_hydrogens:
            Filter.filter_betaHs()
        if self.strict_box_filter:
            Filter.box_excluder_filter()
        if self.filter_even_odd_electron_count is not None:
            Filter.filter_even_odd_electron(filter_for=self.filter_even_odd_electron_count)

        # Denticity dependent filters
        for dent, settings in self.denticity_dependent_filters.items():

            # Get settings for this denticity
            acount_min, acount_max, ligcomp_atoms_of_interest, ligcomp_instruction, \
            ligand_charges, metals_of_interest, coords_atoms_of_interest, coords_instruction, \
            mw_min, mw_max = self.input.check_and_return_denticity_dependent_filter_settings(denticity=dent)

            if metals_of_interest is not None:
                Filter.metals_of_interest_filter(denticity=dent, metals_of_interest=metals_of_interest)
            if coords_atoms_of_interest is not None:
                Filter.filter_coordinating_group_atoms(denticity=dent, atoms_of_interest=coords_atoms_of_interest,
                                                       instruction=coords_instruction)
            if ligcomp_atoms_of_interest is not None:
                Filter.filter_ligand_atoms(denticity=dent, atoms_of_interest=ligcomp_atoms_of_interest,
                                           instruction=ligcomp_instruction)
            if mw_max is not None or mw_min is not None:
                Filter.filter_molecular_weight(denticity=dent, atomic_weight_min=mw_min, atomic_weight_max=mw_max)
            if ligand_charges is not None:
                Filter.filter_ligand_charges(denticity=dent, charge=ligand_charges)
            if acount_min is not None:
                Filter.filter_atom_count(denticity=dent, number=acount_min, instruction="greater_than")
            if acount_max is not None:
                Filter.filter_atom_count(denticity=dent, number=acount_max, instruction="less_than")

            # todo: subgraph match

        return Filter.database

    def save_filtered_ligand_db(self):
        filtered_db = self.get_filtered_db()


        if not self.output_ligand_db_path.parent.exists():
            self.output_ligand_db_path.parent.mkdir(parents=True)
        filtered_db.to_json(self.output_ligand_db_path, json_lines=True)

        return


if __name__ == "__main__":

    ligand_filter_path = project_path().extend('src05_Assembly_Refactor', 'ligandfilters.yml')
    save_db_path = project_path().extend("src14_Assembly_Unit_Test", "ligandfilters", "filtered_ligand_db_v1.7.json")
    max_number = 5000

    filter = LigandFilters(filepath=ligand_filter_path, max_number=max_number)
    filter.save_filtered_ligand_db()

    test = IntegrationTest(new_dir=save_db_path.parent, old_dir=Path(save_db_path.parent.parent, 'ligandfilters_benchmark'))
    test.compare_all()