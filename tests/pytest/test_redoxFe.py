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
def test_redoxFe(outdir: Path) -> Assembler:
    lf = LigandFilters.run_from_yaml(outdir.parent / Path('data_input') / 'redoxFe_ligandfilters.yml')

    return lf

if __name__ == "__main__":

    lf = test_redoxFe()


