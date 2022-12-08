from pathlib import Path
import json
import pandas as pd
import networkx
import ijson

class ComplexDB():

    def __init__(self, json_):
        self.json_path = Path(json_)

    def extract_ligands(self):

        with open(str(self.json_path), 'rb') as input_file:
            jsonobj = ijson.items(input_file, 'item')
            jsons = (o for o in jsonobj)
            for j in jsons:
                print(j)
                print('Hi')

        with open(str(self.json_path), 'rb') as f:
            for n in ijson.items(f, 'categoryAspects.item'):
                print(n)





if __name__ == '__main__':
    all_complexes_json = '../data/tmQMG_Jsons_test/tmQMG.json'  # Folder where we want to store the jsons

    cdb = ComplexDB(all_complexes_json)
    cdb.extract_ligands()
    print('Done!')