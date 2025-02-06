# Here we will load in a ligin database and delete liands that we do not want
from DARTassembler.src.ligand_extraction.DataBase import LigandDB

if __name__ == "__main__":
    data = LigandDB.load_from_json("/Users/cianclarke/Documents/PhD/Complex_Assembly/DART/DART/dev/Assembler_revision_jan_2025/assembler/ligand_filtering/ligandfilters/data_output/hexagonal.json")

    for name, ligand in data.db.items():
        if ligand.eta != 0:
            print("done")
        else:
            pass