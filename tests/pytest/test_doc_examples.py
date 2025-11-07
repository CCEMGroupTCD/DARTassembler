"""
Integration test for the assembly of complexes.
"""
import sys; print('Python %s on %s' % (sys.version, sys.platform))
sys.path.extend(['/Users/timosommer/PhD/projects/DARTassembler'])
from pathlib import Path
from DARTassembler import Assembler, LigandFilters
from DARTassembler.src.misc.tests import integration_test

@integration_test(name='doc_examples')
def test_doc_examples(outdir: Path) -> Assembler:
    n = 1000    # n max ligands, takes precedence over yaml

    C6H6_filters_path = outdir.parent / 'data_input' / 'C6H6_ligandfilters.yml'
    HMDS_filters_path = outdir.parent / 'data_input' / 'amide_ligandfilters.yml'
    assembly_input = outdir.parent / 'data_input' / 'test_doc_examples.yml'

    LigandFilters.run_from_yaml(C6H6_filters_path, n=n)
    LigandFilters.run_from_yaml(HMDS_filters_path, n=n)
    assembly = Assembler.run_from_yaml(assembly_input, n_max_ligands=n)

    return assembly


if __name__ == "__main__":
    assembly = test_doc_examples()