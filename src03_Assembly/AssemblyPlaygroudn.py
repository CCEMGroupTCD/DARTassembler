"""
My File to play arround with the assmebly
"""
from src01.DataBase import LigandDB


if __name__ == "__main__":

    F = LigandDB.from_json(json_="../data/Filtered_Jsons/filteredLigDB.json", type_="Ligand")
