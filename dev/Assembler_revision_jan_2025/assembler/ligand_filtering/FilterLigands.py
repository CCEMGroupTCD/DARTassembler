from DARTassembler.ligandfilters import ligandfilters
import numpy as np

def test_filter_ligands(n_max=np.inf):
    filters = ligandfilters(filter_input_path="YAML/hexagonal.yml", nmax=n_max, delete_output_dir=True)
    return filters


if __name__ == "__main__":
    nmax = np.inf
    filters = test_filter_ligands(n_max=nmax)