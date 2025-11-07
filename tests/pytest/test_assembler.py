"""
Integration test for the assembly of complexes.
"""
import sys
print('Python %s on %s' % (sys.version, sys.platform))
sys.path.extend(['/Users/timosommer/PhD/projects/DARTassembler'])
from DARTassembler.src.misc.tests import integration_test
from pathlib import Path
from DARTassembler.src.assembler.assembler import Assembler


@integration_test(name="assembler")
def test_assembler(outdir: Path) -> object:
    assembly_input = outdir.parent / 'data_input' / 'test_assembly_input.yml'
    return Assembler.run_from_yaml(assembly_input)

if __name__ == "__main__":
    assembly = test_assembler()