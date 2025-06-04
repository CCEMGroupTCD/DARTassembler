from DARTassembler.ligandfilters import ligandfilters

def test_filter_ligands():
    filters = ligandfilters(filter_input_path="small_mono.yml")
    return filters


if __name__ == "__main__":
    filters = test_filter_ligands()