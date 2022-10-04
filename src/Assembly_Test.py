# Todo: Move this file to src03; is here just for the moment due to pickle reasons

import pickle
from src03_Assembly.Assembly_Draft import *
from process import LigandDatabase
import random
from src03_Assembly.assembly_setup import list_of_metals, implemented_topologies


class RandomComplexAssembler:
    """
    is kind of the hub to handle the random assembly.
    Stores the ligand database and stores the configuration for the random assembly
    """

    def __init__(self, database_path: str, store_path: str):

        self.ligand_dict = self.dataload(database_path=database_path)

        self.possible_topologies = implemented_topologies

        try:
            self.assembly_functions = {f"{top}": globals()[get_top_string(top)] for top in self.possible_topologies}
        except Exception as e:
            print("More topologies requested than implemented")
            raise e

        self.list_of_metals = list_of_metals

        self.store_path = store_path

    @staticmethod
    def dataload(database_path):

        ligand_db = pickle.load(open(database_path, "rb"))

        if not isinstance(ligand_db, LigandDatabase):
            print("No correct database selected")
            return {}

        else:
            return ligand_db.ligand_dict

    def create_random_TMC(self, visualize_: bool, optimize_: bool = False):
        try:
            # random choices:
            (metal, charge) = random.choice(self.list_of_metals)
            #
            comp = random.choice(self.possible_topologies)
            #
            ligands = {i: random.choice(self.ligand_dict[index]) for i, index in enumerate(comp)}  # ligands

            #
            # build the center atom with 6 connections
            metal_bb = stk.BuildingBlock(smiles='[Hg+2]',
                                         functional_groups=(stk.SingleAtom(stk.Hg(0, charge=2)) for i in range(6)),
                                         position_matrix=np.ndarray([0, 0, 0])
                                         )

            # build the metal block with the new metal atom
            smiles_str = f"[{metal}{charge}]"
            stk_metal_func = getattr(stk_, metal)
            functional_groups = (stk.SingleAtom(stk_metal_func(0, charge=charge)) for i in range(6))
            final_metal_bb = stk.BuildingBlock(smiles=smiles_str,
                                               functional_groups=functional_groups,
                                               position_matrix=np.ndarray([0, 0, 0])
                                               )

            ligand_bb_dict = {}

            for k, lig in enumerate(ligands.values()):
                # fill ligand_bb_dict as  {nr : (ligand, ligand_bb)}
                lig_assembly_dict = lig.get_assembly_dict()

                ligand_to_mol(ligand=lig, xyz_path="../tmp/lig_xyz.xyz", target_path="../tmp/lig_mol.mol")

                ligand_bb_dict[k] = (lig, build_ligand(type_list=lig_assembly_dict["type"],
                                                       index_list=lig_assembly_dict["index"],
                                                       path_="../tmp/lig_mol.mol"))
                os.remove("../tmp/lig_mol.mol")

            complex_ = self.assembly_functions[str(comp)](metal_bb=metal_bb,
                                                          final_metal_bb=final_metal_bb,
                                                          ligand_bb_dict=ligand_bb_dict,
                                                          optimize_=optimize_,
                                                          planar_=planar_ceck(ligand_bb_dict)
                                                          )

            if complex_ is not None:
                remove_Hg(complex_,
                          name=f"{ligands[0].csd_code}",
                          visualize_=visualize_,
                          print_to_xyz=True,
                          path=self.store_path
                          )

            return complex_

        except rdkit.Chem.rdchem.AtomValenceException as ex:
            print(f"The standard Error: {ex} has occured. Solveable by turning sanitize in stk off")
            return None

        except Exception as e:
            print(f"Oh. A new error: {e}!!")
            return None


if __name__ == "__main__":

    RCA = RandomComplexAssembler(database_path="../data/ligand_db.pickle",
                                 store_path="../data/Assembled_Molecules")

    for i in range(10):
        random_complex = RCA.create_random_TMC(visualize_=True, optimize_=False)
        input("press Enter to continue")

    print("done")
