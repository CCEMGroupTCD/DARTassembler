"""
The script for extracting all mol2 files from the CSD and saving them batch wise (in orde to keep them below 100Mb each
-> Githubs limit)
"""
from ccdc import io
from tqdm import tqdm


def CSD_2_mol2():
    csd_reader = io.EntryReader("CSD")
    batch_nr = 0

    for i, el in tqdm(enumerate(csd_reader.molecules())):

        if i % 40000 == 0:
            batch_nr += 1

        try:
            if batch_nr < 10:
                with open(f"../data/Batch_0{batch_nr}.mol2", "a+") as file:
                    file.write("\n\n\n")
                    file.write(el.to_string("mol2"))
            else:
                with open(f"../data/Batch_{batch_nr}.mol2", "a+") as file:

                    file.write("\n\n\n")
                    file.write(el.to_string("mol2"))

        except Exception as e:
            print(f"Something went wrong {e}")


if __name__ == "__main__":
    CSD_2_mol2()
