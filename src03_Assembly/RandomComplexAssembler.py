"""
Alt, alles neue wird in RandomComplexAssmebler ausprobiert
"""

import random
import yaml
from yaml import SafeLoader
import stk

stk_ = __import__("stk")

from src03_Assembly.Topology_Assembly_Methods import *
from src03_Assembly.utilities_stk import create_placeholder_Hg_bb, build_ligand, remove_Hg
from src03_Assembly.utilities_assembly import get_topology_string
from src03_Assembly.TransitionMetalComplex import TransitionMetalComplex as TMC

from src01.DataBase import LigandDB


class RandomComplexAssembler:
    """
    is kind of the hub to handle the random assembly.
    Stores the ligand database and stores the configuration for the random assembly
    """

    def __init__(self, database: LigandDB, store_path: str = "../data/Assembled_Molecules"):

        self.ligand_dict = database.get_lig_db_in_old_format()

        with open("assembly_setup.yml", "r") as handle:
            setup = yaml.load(handle, SafeLoader)

        self.possible_topologies = setup["implemented_topologies"]
        self.list_of_metals = setup["list_of_metals"]

        try:
            self.assembly_functions = {f"{top}": globals()[get_topology_string(top)] for top in
                                       self.possible_topologies}
        except Exception as e:
            print("More topologies requested than implemented")
            raise e

        self.store_path = store_path

        print("RCA init complete")

    def choose_random_parts(self, rs):
        if rs is not None:
            random.seed(rs)

        # random choices:
        metal, charge = random.choice(self.list_of_metals)
        #
        comp = random.choice(self.possible_topologies)
        #
        ligands = {i: random.choice(self.ligand_dict[index]) for i, index in enumerate(comp)}  # ligands

        return metal, charge, ligands, comp

    @staticmethod
    def create_metal_building_block(metal, charge) -> stk.BuildingBlock:

        # build the metal block with the new metal atom
        smiles_str = f"[{metal}{charge}]"
        stk_metal_func = getattr(stk_, metal)
        functional_groups = (stk.SingleAtom(stk_metal_func(0, charge=charge)) for i in range(6))
        final_metal_bb = stk.BuildingBlock(smiles=smiles_str,
                                           functional_groups=functional_groups,
                                           position_matrix=np.ndarray([0, 0, 0])
                                           )

        return final_metal_bb

    @staticmethod
    def create_ligand_building_blocks(ligands) -> dict:
        """
        # todo: GLaube an den LigandBuildingBlock dicts kann man noch arbeiten. Vielleicht ist das aber
        auch best moeglich so
        """
        ligand_bbs = {}

        for k, lig in enumerate(ligands.values()):
            # fill ligand_bb_dict as  {nr : (ligand, ligand_bb)}
            lig_assembly_dict = lig.get_assembly_dict()

            ligand_to_mol(ligand=lig, xyz_path="../tmp/lig_xyz.xyz", target_path="../tmp/lig_mol.mol")

            ligand_bbs[k] = (lig, build_ligand(type_list=lig_assembly_dict["type"],
                                               index_list=lig_assembly_dict["index"],
                                               path_="../tmp/lig_mol.mol"))
            os.remove("../tmp/lig_mol.mol")

        return ligand_bbs

    @staticmethod
    def planar_ceck(ligand_bb_dict):
        """
        If a tri or tetradentate ligand is contained in the topology, we need to evaluate whether this is
        a planar one or not.
        returns True if yes and false otherwise
        """
        for key, (lig, lig_bb) in ligand_bb_dict.items():
            if lig.denticity == 4 or lig.denticity == 3:
                return lig.planar_check()

        return False

    def create_random_TMC(self, visualize_: bool, optimize_: bool = False, random_seed: int = None):
        """
        Assembles a transition metal complex at random
        visualize
        """

        # First we choose elements at random
        metal, charge, ligands, topology = self.choose_random_parts(rs=random_seed)

        # With the so choosen metal we create the metal_building block based on the randomly choosen metal
        final_metal_bb = self.create_metal_building_block(metal, charge)

        # Now we need also to convert the ligands into building blocks themselves
        ligand_bb_dict = self.create_ligand_building_blocks(ligands)

        # Further, we need to check if the important ligand (tri or tetra) is planar or not
        planar = self.planar_ceck(ligand_bb_dict)

        # and thus we got everything to assemble the complex
        complex_ = self.assembly_functions[str(topology)](metal_bb=create_placeholder_Hg_bb(),
                                                          final_metal_bb=final_metal_bb,
                                                          ligand_bb_dict=ligand_bb_dict,
                                                          # optimize_=optimize_,            # todo: Talk this through with Cian
                                                          planar_=planar
                                                          )

        if complex_ is not None:
            # todo: Also hier muss ich auf jeden Fall auch nochmal ran, wie wir hier die Infos behalten
            #   das ist noch ein suboptimaler Datenflow
            # if we were able to create a complex it still remains to remove the placeholder Hg atoms
            # We return an RCA molecule
            ase_complex, xyz_str = remove_Hg(complex_, visualize_=visualize_)
            return TMC(mol=ase_complex, metal_symbol=metal, ligands=ligands, xyz_str=xyz_str)
        else:
            return complex_  # which will be None
