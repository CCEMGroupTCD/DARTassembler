# Here we will load in a ligin database and delete liands that we do not want
from DARTassembler.src.ligand_extraction.DataBase import LigandDB
from copy import deepcopy

if __name__ == "__main__":
    old_data = LigandDB.load_from_json("/Users/cianclarke/Documents/PhD/Complex_Assembly/DART/DART/dev/Assembler_revision_jan_2025/assembler/ligand_filtering/JSON/hexagonal.json")
    new_data = deepcopy(old_data)
    print(len(new_data.db))

    """
    # identify OH ligand
    for name, ligand in old_data.db.items():
        if (ligand.kappa == 1) and (ligand.pred_charge == -1) and (len(ligand.atomic_props["atoms"]) == 2) and ("O" in ligand.atomic_props["atoms"]) and ("H" in ligand.atomic_props["atoms"]):
            pass
        else:
            new_data.db.pop(name)"""




    # identify planar hexadenate ligands
    for name, ligand in old_data.db.items():
        print(ligand.geometry)
        if ligand.eta == 6 and ligand.kappa == 0 and ligand.geometry == "1_monodentate":
            pass
        else:
            new_data.db.pop(name)

    print(len(new_data.db))
    new_data.to_json("/Users/cianclarke/Documents/PhD/Complex_Assembly/DART/DART/dev/Assembler_revision_jan_2025/assembler/ligand_filtering/JSON/hexagonal_eta.json")