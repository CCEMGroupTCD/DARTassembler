"""
Integration test for the filtering and assembly of Au(III) complexes (case study for the DART paper).
"""
import sys
print('Python %s on %s' % (sys.version, sys.platform))
sys.path.extend(['/Users/timosommer/PhD/projects/DARTassembler'])
from pathlib import Path
from DARTassembler import Assembler, LigandFilters
from DARTassembler.src.misc.tests import integration_test

@integration_test(name='redoxFe')
def test_redoxFe(outdir: Path) -> tuple[LigandFilters, Assembler]:
    lf = LigandFilters.run_from_yaml(outdir.parent / 'data_input' / 'redoxFe_ligandfilters.yml')

    # Run the assembler and save output to new assembler
    assembler = Assembler.run_from_yaml(outdir.parent / 'data_input' / 'redoxFe_assembly_input.yaml')

    return lf, assembler

if __name__ == "__main__":

    lf, assembler = test_redoxFe()


