import random
import yaml
from yaml import SafeLoader
import logging
import stk
import numpy as np

from src03_Assembly.building_block_utility import rotate_tridentate_bb, rotate_tetradentate_bb, penta_as_tetra, \
    get_optimal_rotation_angle_tridentate, Bidentate_Rotator, nonplanar_tetra_solver
from src03_Assembly.stk_utils import create_placeholder_Hg_bb
from src03_Assembly.TransitionMetalComplex import TransitionMetalComplex as TMC
import src03_Assembly.stk_extension as stk_e

from src01.DataBase import LigandDB
from src01.Molecule import RCA_Ligand


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
        stk_metal_func = getattr(stk, metal)
        functional_groups = (stk.SingleAtom(stk_metal_func(0, charge=charge)) for i in range(6))
        final_metal_bb = stk.BuildingBlock(smiles=smiles_str,
                                           functional_groups=functional_groups,
                                           position_matrix=np.ndarray([0, 0, 0])
                                           )

        return final_metal_bb

    @staticmethod
    def planar_ceck(ligands):
        """
        If a tri or tetradentate ligand is contained in the topology, we need to evaluate whether this is
        a planar one or not.
        returns True if yes and false otherwise
        """
        for key, lig in ligands.items():
            if lig.denticity == 4 or lig.denticity == 3:
                return lig.planar_check()

        return True

    def convert_ligand_to_building_block_for_complex(self, ligands: dict[RCA_Ligand]) -> dict:

        # initial parameters
        topology_determining_ligand_planar = self.planar_ceck(ligands)

        # now the actual conversion
        ligand_buildingblocks = {}
        for i, ligand in enumerate(ligands.values()):

            if ligand.denticity == 0:
                # this is the reactant, which in our case will always be monodentate
                monodentate_topology = stk_e.Monodentate(metals=create_placeholder_Hg_bb(),
                                                         ligands=ligand.to_stk_bb())
                if topology_determining_ligand_planar is False:
                    monodentate_topology.set_ligand_coordinates(coordinates=[-1.2, 1.2, 0])

                bb_for_complex = stk.BuildingBlock.init_from_molecule(stk.ConstructedMolecule(
                    topology_graph=monodentate_topology),
                    functional_groups=[stk.SmartsFunctionalGroupFactory(smarts='[Hg+2]', bonders=(0,), deleters=())]
                )

            elif ligand.denticity == 1:
                # only chance, we are in 4-1-0
                monodentate_topology = stk_e.Monodentate(metals=create_placeholder_Hg_bb(),
                                                         ligands=ligand.to_stk_bb())
                if topology_determining_ligand_planar is False:
                    monodentate_topology.set_ligand_coordinates(coordinates=[-1.2, -1.2, 0])
                else:
                    monodentate_topology.set_ligand_coordinates(coordinates=[0, 0, -1.9])

                bb_for_complex = stk.BuildingBlock.init_from_molecule(stk.ConstructedMolecule(
                    topology_graph=monodentate_topology),
                    functional_groups=[stk.SmartsFunctionalGroupFactory(smarts='[Hg+2]', bonders=(0,), deleters=())]
                )

            elif ligand.denticity == 4:
                building_block = ligand.to_stk_bb()

                if topology_determining_ligand_planar is True:
                    # Then some rotation needs to be done
                    building_block = rotate_tetradentate_bb(building_block, ligand_=ligand)

                    tetra_topology_graph = stk.metal_complex.Porphyrin(metals=create_placeholder_Hg_bb(),
                                                                       ligands=building_block
                                                                       )

                    bb_for_complex = stk.BuildingBlock.init_from_molecule(
                        stk.ConstructedMolecule(topology_graph=tetra_topology_graph),
                        functional_groups=[
                            stk.SmartsFunctionalGroupFactory(smarts='[Hg+2]', bonders=(0,), deleters=(), )]
                    )

                else:
                    # todo: Das muss noch gemacht werden
                    new_bb_path = nonplanar_tetra_solver(bb=building_block,
                                                         lig=ligand)

                    bb_for_complex = stk.BuildingBlock.init_from_file(new_bb_path, functional_groups=[
                        stk.SmartsFunctionalGroupFactory(
                            smarts='[Hg]', bonders=(0,),
                            deleters=(), ), ], )

            elif ligand.denticity == 3:

                building_block = ligand.to_stk_bb()

                if topology_determining_ligand_planar is False:
                    print("Not implemented yet, we should not end up here")
                    raise NotImplementedError
                else:
                    building_block = rotate_tridentate_bb(tridentate_bb_=building_block, ligand_=ligand)

                    tridentate_toplogy = stk_e.Tridentate(metals=create_placeholder_Hg_bb(),
                                                          ligands=building_block
                                                          )

                    compl_constructed_mol = stk.ConstructedMolecule(topology_graph=tridentate_toplogy)

                    compl_constructed_mol = compl_constructed_mol.with_rotation_about_axis(
                        axis=np.array((0, 0, 1)),
                        angle=float(np.radians(
                            get_optimal_rotation_angle_tridentate(compl_constructed_mol, 10.0, 0.0, 0.0, ligand))),
                        origin=np.array((0, 0, 0))
                    )

                    bb_for_complex = stk.BuildingBlock.init_from_molecule(
                        compl_constructed_mol,
                        functional_groups=[stk.SmartsFunctionalGroupFactory(smarts='[Hg+2]', bonders=(0,), deleters=())]
                    )

            elif ligand.denticity == 2:

                bidentate_topology = stk_e.Bidentate(metals=create_placeholder_Hg_bb(),
                                                     ligands=ligand.to_stk_bb()
                                                     )

                # todo: Hier muss ich auch nochmal ran
                complex_bidentate = stk.ConstructedMolecule(topology_graph=bidentate_topology)

                new_mol_path = Bidentate_Rotator(ligand_bb=complex_bidentate, ligand=ligand)

                bb_for_complex = stk.BuildingBlock.init_from_file(new_mol_path, functional_groups=[
                    stk.SmartsFunctionalGroupFactory(smarts='[Hg+2]', bonders=(0,), deleters=(), ), ], )

            elif ligand.denticity == 5:
                """
                The whole trick is to identify a pentadentate as a planar tetradentate ligand and place it correctly.
                # We first aim to treat a penta dentate as a planar tetradentate and then follow the process above for a planar
                # tetra dentate
                """

                tetra_bb_for_penta, position_index = penta_as_tetra(ligand=ligand)

                tetra_bb_for_penta = rotate_tetradentate_bb(tetra_bb_for_penta, ligand)

                tip_position = list(tetra_bb_for_penta.get_atomic_positions(atom_ids=[int(position_index), ]))

                if float(tip_position[0][2]) > 0:
                    # an additional rotation is required
                    tetra_bb_for_penta = tetra_bb_for_penta.with_rotation_about_axis(angle=np.radians(180),
                                                                                     axis=np.array((1, 0, 0)),
                                                                                     origin=np.array((0, 0, 0))
                                                                                     )

                elif float(tip_position[0][2]) < 0:
                    # no rotation required
                    pass
                else:
                    print("Something went wrong")
                    raise ValueError

                penta_topology = stk.metal_complex.Porphyrin(metals=create_placeholder_Hg_bb(),
                                                             ligands=tetra_bb_for_penta
                                                             )

                bb_for_complex = stk.BuildingBlock.init_from_molecule(
                    stk.ConstructedMolecule(topology_graph=penta_topology),
                    functional_groups=[stk.SmartsFunctionalGroupFactory(smarts='[Hg+2]', bonders=(0,), deleters=(), ), ]
                )

            else:
                print("Unknown Ligand Denticity, something went wrong")
                raise ValueError

            ligand_buildingblocks[i] = bb_for_complex

        return ligand_buildingblocks

    def building_block_assembly(self,
                                ligands: dict,
                                metal: str,
                                charge: str,
                                ) -> stk.ConstructedMolecule:
        """
        We now generalize the whole assembly process
        """

        building_blocks = self.convert_ligand_to_building_block_for_complex(ligands)

        if len(building_blocks) == 3:
            complex_top = stk_e.complex_topology_three(metals=self.create_metal_building_block(metal, charge),
                                                       ligands={building_block: (i,) for i, building_block in
                                                                building_blocks.items()}
                                                       )
        else:
            complex_top = stk_e.complex_topology_two(metals=self.create_metal_building_block(metal, charge),
                                                     # kann einfach dann die klassenmethode fuer metall_bb werdedn
                                                     ligands={building_block: (i,) for i, building_block in
                                                              building_blocks.items()}
                                                     )

        complex_ = stk.ConstructedMolecule(topology_graph=complex_top)

        return complex_

    def create_random_TMC(self, random_seed: int = None):
        """
        Assembles a transition metal complex at random
        What this method does in brief:

        Choose random metals and ligands from the underlying RCA database
        convert them to stk.Building Blocks
        and assembles them according to the approiate topology
        """

        # First we choose elements at random
        metal, charge, ligands, topology = self.choose_random_parts(rs=random_seed)

        # Further, we need to check if the important ligand (tri or tetra) is planar or not
        planar = self.planar_ceck(ligands)  # is True fo 5-0!

        if topology == [3, 2, 0] and planar is False:
            complex_ = None
            logging.info("Not implemented yet")
            return complex_
        else:
            # and thus we got everything to assemble the complex
            complex_ = self.building_block_assembly(ligands, metal, charge)

            return TMC(compl=complex_, ligands=ligands, metal=metal, metal_charge=int(charge))
