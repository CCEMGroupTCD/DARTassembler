import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from DARTassembler.src.ligand_extraction.utilities_Molecule import stoichiometry2atomslist
from DARTassembler.src.constants.Periodic_Table import DART_Element

if __name__ == '__main__':

    csv_path = '/Users/timosommer/Downloads/test_DART_install/quickstart/DART_output/info_table.csv'
    outfig_path = 'plots/ligand_elements_without_filtering.png'

    # Load the data
    df = pd.read_csv(csv_path)

    df = df[df['success']] # Filter out failed complexes

    atom_lists_for_each_complex = df['stoichiometry'].apply(lambda x: list(set(stoichiometry2atomslist(x)))).tolist()
    # Flatten list of lists
    element_list = [item for sublist in atom_lists_for_each_complex for item in sublist]
    # Remove metals
    element_list = [element for element in element_list if not DART_Element(element).is_metal]
    element_hist = pd.Series(element_list).value_counts()

    # Plot the ligand elements as histogram
    plt.figure()
    sns.barplot(x=element_hist.index, y=element_hist.values)
    plt.xlabel('Element')
    plt.ylabel('Count')
    plt.savefig(outfig_path, dpi=300, transparent=True)