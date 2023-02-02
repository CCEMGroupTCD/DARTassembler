from openbabel import openbabel as ob
from openbabel import pybel
import matplotlib.pyplot as plt
import os
from tqdm import tqdm

Func_indexes = [43, 46, 44, 41, 2, 1]
Func_indexes2 = [25, 24, 23, 2, 3, 21]


def reset_coords(file_read, file_write, indexes):
    with open(file_read, "r") as f1, open(file_write, "r+") as f2:
        i = 1
        data = []
        for line_read1, line_read2 in zip(f1.readlines(), f2.readlines()):
            if any(int(i - 5) == index for index in indexes):
                data.append(line_read1)
            else:
                data.append(line_read2)
            i = i + 1
    os.system('touch  tmp_uff_tmp.mol')
    with open('tmp_uff_tmp.mol', "w") as f3:
        for item in data:
            f3.write(item)
    f1.close()
    f2.close()
    f3.close()
    os.system('mv tmp_uff_tmp.mol tmp_uff.mol')
    return None


def movie(input_file):
    os.system('obabel .mol {} .xyz -O  movie_input.xyz'.format(input_file))
    os.system('touch  movie.xyz')
    os.system('cat movie_input.xyz >> movie.xyz ')


def get_energy(input_file):
    mol = next(pybel.readfile("mol", input_file))
    obmol = mol.OBMol
    ff = ob.OBForceField_FindType("uff")
    assert (ff.Setup(obmol))
    kj_to_kcal = 1.0 / 4.184
    ff.SetCoordinates(mol.OBMol)
    uffE = ff.Energy(False) * kj_to_kcal
    return uffE


def optimiser(input_file, output_file, n):
    mol = next(pybel.readfile("mol", "{}".format(input_file)))
    mol.localopt(forcefield='uff', steps=n)
    mol.write("mol", "{}".format(output_file), overwrite="True")


def optimise_complex_step1(input_file, output_file, setting):
    string = 'obabel .xyz {} .mol -O  tmp.mol'.format(input_file)
    os.system(string)
    if setting == "lock":
        i = 0
        for i in range(300):
            print(i)
            if i == 0:
                optimiser("tmp.mol", "tmp_uff.mol", n=1)                # tmp.mol is your original file with the original coordinates and tmp_uff.mol is the output_file
                reset_coords("tmp.mol", "tmp_uff.mol", Func_indexes2)
            else:
                optimiser("tmp_uff.mol", "tmp_uff.mol", n=1)
                reset_coords("tmp.mol", "tmp_uff.mol", Func_indexes2)
        optimiser("tmp_uff.mol", "tmp_uff_2.mol", n=1000)
        os.system("mv tmp_uff_2.mol {}".format(output_file))
        os.system("rm tmp_uff.mol")
    else:
        optimiser("tmp.mol", "tmp_uff.mol", n=1000)


if __name__ == "__main__":
    optimise_complex_step1("tmp2.xyz", "tmp_uff_step1.mol", "lock")

"""
optimise_complex("tmp2.xyz")
os.system("ase-gui movie.xyz")
os.system("rm -f movie_input.xyz")
#### REMBER TO DELETE THE MOVIE FILE AFTER EVERY RUN ####
# ['gaff', 'ghemical', 'mmff94', 'mmff94s', 'uff']
"""
plt.show()
