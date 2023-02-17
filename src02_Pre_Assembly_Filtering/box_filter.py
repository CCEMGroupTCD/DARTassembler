from src01.Molecule import RCA_Molecule, RCA_Ligand
from src04_Assembly.building_block_utility import *
from src03_Pre_Assembly_Filter.constants_BoxExcluder import get_boxes
from ase import io
from stk import *

import logging
import stk
import numpy as np

from src04_Assembly.building_block_utility import rotate_tridentate_bb, rotate_tetradentate_bb, penta_as_tetra, \
    get_optimal_rotation_angle_tridentate, Bidentate_Rotator, nonplanar_tetra_solver, get_energy_stk
from src04_Assembly.stk_utils import create_placeholder_Hg_bb
import src04_Assembly.stk_extension as stk_e




def delete_Hg_improved(complex_):

    mol = complex_.to_rdkit_mol()
    Hg_index = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 80:
            Hg_index.append(int(atom.GetIdx()))
        else:
            pass
    mol_RW = Chem.EditableMol(mol)
    for index in Hg_index:
        mol_RW.RemoveAtom(index - Hg_index.index(index))
    non_editable_mol = mol_RW.GetMol()
    output_stk_bb = stk.BuildingBlock.init_from_rdkit_mol(non_editable_mol)
    return output_stk_bb


def box_excluder_descision(stk_Building_Block, denticity_, planar=True, threshhold=0.):
    boxes = get_boxes(denticity=denticity_, planar=planar, input_topology=None)
    atom_positions = list(stk_Building_Block.get_atomic_positions())

    score = 0

    for atom in atom_positions:  # loop through atoms
        point = [atom[j] for j in range(3)]

        for Box in boxes:
            if Box.point_in_box(point):
                score += Box.intensity / (1.0 + (
                        Box.sharpness * ((point[0]) - ((Box.x2 - Box.x1) / 2.0) + Box.x1) ** 2))
                score += Box.intensity / (1.0 + (
                        Box.sharpness * ((point[1]) - ((Box.y2 - Box.y1) / 2.0) + Box.y1) ** 2))
                score += Box.intensity / (1.0 + (
                        Box.sharpness * ((point[2]) - ((Box.z2 - Box.z1) / 2.0) + Box.z1) ** 2))
                break
    print("score: "+str(score))
    if score > threshhold:
        print("ligand has failed")
        # did not pass
        return False
    else:
        print("ligand has passed")
        return True


def box_filter(ligand: RCA_Ligand, box_default_descicion=False) -> bool:
    """
    Returns True if a ligand looks good and can pass
    and false if not
    """
    try:
        print("################################################################################################################################################")
        metal, charge = "Fe", "+2"

        # build the metal block with the new metal atom
        smiles_str = f"[{metal}{charge}]"
        stk_metal_func = globals()[f"{metal}"]
        # stk_metal_func = getattr(__import__("stk"), metal)
        functional_groups = (stk.SingleAtom(stk_metal_func(0, charge=charge)) for i in range(6))
        final_metal_bb = stk.BuildingBlock(smiles=smiles_str,
                                           functional_groups=functional_groups,
                                           position_matrix=np.ndarray([0, 0, 0])
                                           )

        ligand_bb = ligand.to_stk_bb()

        print("The ligand denticity is: " + str(ligand.denticity))

        # build the building blocks
        if len(ligand.coordinates) == 2 and ligand.denticity > 1:
            print("!!!Warning!!! -> A ligand with 2 atoms was detected -> skipping ligand")
            raise ValueError
        if ligand.denticity == 3 and ligand.planar_check() is True:
            print("I have encountered a tridentate ligand")
            building_block = rotate_tridentate_bb(tridentate_bb_=ligand_bb, ligand_=ligand)

            tridentate_topology = stk_e.Tridentate(metals=create_placeholder_Hg_bb(),
                                                   ligands=building_block
                                                   )

            compl_constructed_mol = stk.ConstructedMolecule(topology_graph=tridentate_topology)

            compl_constructed_mol = compl_constructed_mol.with_rotation_about_axis(
                axis=np.array((0, 0, 1)),
                angle=float(np.radians(
                    get_optimal_rotation_angle_tridentate(compl_constructed_mol, 10.0, 0.0, 0.0, ligand))),
                origin=np.array((0, 0, 0))
            )
            position_matrix = compl_constructed_mol.get_position_matrix()
            position_matrix[0] = [-0.6, 0, 0]
            compl_constructed_mol = compl_constructed_mol.with_position_matrix(position_matrix=position_matrix)
            bb_for_comp = stk.BuildingBlock.init_from_molecule(compl_constructed_mol, functional_groups=[stk.SmartsFunctionalGroupFactory(smarts='[Hg+2]', bonders=(0,), deleters=())])

        elif ligand.denticity == 3 and ligand.planar_check() is False:
            logging.info("Not implemented yet")
            print("Not implemented yet")
            # I am a little concerned by this, letting non_planar_tridentates through
            return True  # will always let it through


        elif ligand.denticity == 4 and ligand.planar_check() is True:
            print("I have encountered a tetradentate ligand")
            # so essentially what I need to do here is rotate the molecule to the postion it should be in
            ligand_bb = rotate_tetradentate_bb(ligand_bb, ligand_=ligand)

            tetra_topology_graph = stk.metal_complex.Porphyrin(metals=create_placeholder_Hg_bb(),
                                                               ligands=ligand_bb
                                                               )

            bb_for_comp = stk.BuildingBlock.init_from_molecule(stk.ConstructedMolecule(topology_graph=tetra_topology_graph),
                                                               functional_groups=[stk.SmartsFunctionalGroupFactory(smarts='[Hg+2]', bonders=(0,), deleters=(), )])

        elif ligand.denticity == 4 and ligand.planar_check() is False:
            print("Not implemented yet")
            return True  # will always let it through


        elif ligand.denticity == 5:
            print("I have encountered a pentadentate ligand")
            tetra_bb_for_penta, position_index = penta_as_tetra(ligand=ligand)

            tetra_bb_for_penta = rotate_tetradentate_bb(tetra_bb_for_penta, ligand)

            tip_position = list(tetra_bb_for_penta.get_atomic_positions(atom_ids=[int(position_index), ]))

            if float(tip_position[0][2]) > 0:
                # Additional rotation is required
                tetra_bb_for_penta = tetra_bb_for_penta.with_rotation_about_axis(angle=np.radians(180), axis=np.array((1, 0, 0)), origin=np.array((0, 0, 0)))
            elif float(tip_position[0][2]) < 0:
                # No rotation is required
                pass
            else:
                print("Tip position for penta = 0, unexpected")
                raise ValueError

            penta_topology = stk.metal_complex.Porphyrin(metals=create_placeholder_Hg_bb(),
                                                         ligands=tetra_bb_for_penta
                                                         )

            bb_for_comp = stk.BuildingBlock.init_from_molecule(
                stk.ConstructedMolecule(topology_graph=penta_topology),
                functional_groups=[stk.SmartsFunctionalGroupFactory(smarts='[Hg+2]', bonders=(0,), deleters=(), ), ]
            )

        elif ligand.denticity == 2:
            ligand_name = str(ligand.name)
            print("The ligand name is " + str(ligand_name))
            print("I have encountered a bidentate ligand")
            # Rember that the Bidentate_planar_Right_needs to be edited to accomidate a bidentate ligand in different positions
            bidentate_topology = stk_e.Bidentate_Planar_Right(metals=create_placeholder_Hg_bb(), ligands=ligand_bb)
            complex_bidentate = stk.ConstructedMolecule(topology_graph=bidentate_topology)

            print("The number of atoms are " + str(complex_bidentate.get_atomic_positions()))
            new_mol_path = Bidentate_Rotator(ligand_bb=complex_bidentate, ligand=ligand, top_list=[2, 2], bool_placed=False)
            bb_for_comp = stk.BuildingBlock.init_from_file(new_mol_path, functional_groups=[stk.SmartsFunctionalGroupFactory(smarts='[Hg+2]', bonders=(0,), deleters=(), ), ], )

        else:
            print("Something went wrong")
            raise ValueError

        # convert the building blocks to topologies
        if ligand.denticity == 5:
            complex_ = stk.ConstructedMolecule(topology_graph=complex_topology_two(metals=final_metal_bb, ligands=bb_for_comp))
            complex_ = delete_Hg_improved(complex_)
            print("Analysis completed_1")
            return box_excluder_descision(complex_, denticity_=ligand.denticity)


        else:
            complex_ = stk.ConstructedMolecule(topology_graph=complex_topology_three(metals=final_metal_bb, ligands=bb_for_comp))
            complex_ = delete_Hg_improved(complex_)
            print("Analysis completed_2")
            return box_excluder_descision(complex_, denticity_=ligand.denticity, planar=ligand.planar_check())

    except Exception as e:
        print(f"Box Excluder filter failed, will return default ({box_default_descicion}) (I.e. dont let it pass)"
              f"\n Reason for failing: {e}")
        return box_default_descicion
