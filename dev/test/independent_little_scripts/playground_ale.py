import pandas as pd
from DARTassembler.src.ligand_extraction.io_custom import load_json
from DARTassembler.src.ligand_extraction.utilities_graph import graph_from_graph_dict

if __name__ == '__main__':

    n_complexes = 100      # None for all complexes
    dataset = f'../../data/final_db_versions/complex_db_v1.7.json'  # all complexes

    print('Load complex database.')
    db_dict = load_json(dataset, n_max=n_complexes)
    df = pd.DataFrame.from_dict(db_dict, orient='index')
    df['graph'] = df['graph_dict'].apply(lambda graph_dict: graph_from_graph_dict(graph_dict))

    # Important transition metal complex properties:
    # - metal_oxi_state
    # - charge
    # - stoichiometry
    # - atomic_props
    # - graph_dict
    # - metal_atomic_number

    # features:
    # - 3D structure:
        #   - SOAP descriptors (in package dscribe)
    # - stoichiometry:
        #   - MAGPIE descriptors (in package chemml, maybe also in matminer)
    # target:
        # - metal_oxi_state

    # Important ML steps:
    # - chose a model
    # - read in data
    # - split into train and test set
    # - preprocess data
    #   - standard scaling
    # - train model
    # - evaluate model
    #   - metrics: r2, MAE for regression, accuracy for classification
    #   - for both train and test dataset
    # - save model
    # - save results, i.e. predictions, errors, metrics etc.
