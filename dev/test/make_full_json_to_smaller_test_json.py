"""
This script just takes a json path and a number as input and outputs a smaller json with only the chosen number of items. Useful to get a json for debugging purposes from the whole big json file.
"""
import json

if __name__ == '__main__':

    input_json = '../data/tmQMG_Jsons_fixed_gbl_props_cutoffs_full/tmQMG.json'
    num_output_items = 1000
    output_json = '../data/tmQMG_Jsons_fixed_gbl_props_cutoffs_test/tmQMG.json'



    # Read in json.
    with open(input_json, 'r') as file:
        print('Read in json.')
        d = json.load(file)
    if len(d) < num_output_items:
        raise ValueError('Input dictionary is smaller than the chosen output size. Please choose a smaller num_output_items.')

    # Reduce to chosen number of items
    d = {key: value for i, (key, value) in enumerate(d.items()) if i < num_output_items}
    assert len(d) == num_output_items

    # Save smaller json.
    with open(output_json , 'w') as file:
        print(f'Save smaller json with only {num_output_items} items.')
        json.dump(d, file)

    print('Done!')

