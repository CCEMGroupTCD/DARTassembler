import stk
import os
import ast

from DARTassembler.src.assembler.utils import generate_pronounceable_word


def format_topologies(top_string):
    # This simply splits the topology and similarity(instruction) lists
    output_list = str(top_string).split("--")
    topology = output_list[0]
    instruction = output_list[1]
    topology = ast.literal_eval(topology)
    instruction = ast.literal_eval(instruction)
    return topology, instruction

def get_lig_db_in_old_format(db):
    """
    sort it by denticity
    """
    lig_dict_old_format = {}

    for lig in db.values():

        denticity = lig['denticity']
        if denticity not in lig_dict_old_format:
            lig_dict_old_format[denticity] = [lig]
        else:
            lig_dict_old_format[denticity].append(lig)

    return lig_dict_old_format

def visualize(input_complex):
    # This method allows to visualize in a blocking way during debug but is not essential at all
    print("initializing visualization")
    stk.MolWriter().write(input_complex, 'input_complex.mol')
    os.system('obabel .mol input_complex.mol .xyz -O  output_complex.xyz')
    os.system("ase gui output_complex.xyz")
    os.system("rm -f input_complex.mol")
    os.system("rm -f output_complex.xyz")
    print("visualization complete")



