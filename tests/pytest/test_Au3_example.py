"""
Integration test for the filtering and assembly of Au(III) complexes (case study for the DART paper).
"""
import sys
print('Python %s on %s' % (sys.version, sys.platform))
sys.path.extend(['/Users/timosommer/PhD/projects/DARTassembler'])
from pathlib import Path
from DARTassembler import Assembler, LigandFilters
from DARTassembler.src.misc.tests import integration_test

@integration_test(name='Au3_example')
def test_Au3_example(outdir: Path) -> Assembler:
    n = 5000
    Cl_filters = [
        {'filter': 'composition', 'elements': 'Cl', 'instruction': 'must_contain_and_only_contain', 'only_donors': False},
        {'filter': 'property', 'name': 'archetype', 'values': ['1-mono']}
    ]
    CN_filters = [
        {'filter': 'property', 'name': 'n_denticities', 'values': [2]},
        {'filter': 'property', 'name': 'n_haptic_groups', 'values': [0]},
        {'filter': 'property', 'name': 'charge', 'values': [-1]},
        {'filter': 'property', 'name': 'archetype', 'values': ['2-cis']},    # This filter was not in the original study, we just manually removed the 2-trans ligand afterwards
        {'filter': 'composition', 'elements': 'CN', 'instruction': 'must_contain_and_only_contain', 'only_donors': True},
    ]

    # Run the filters and save filtered ligand databases to outdir
    filter = LigandFilters(db='metalig', n=n)
    filter.run(Cl_filters, outpath=Path('ligands')/'Cl.jsonlines', dbinfo=False)
    filter.run(CN_filters, outpath=Path('ligands')/'CN.jsonlines', dbinfo=False)

    # Run the assembler and save output to new assembler
    assembly_input = outdir.parent / 'data_input' / 'Au3_assembly_input.yaml'
    assembler = Assembler.run_from_yaml(assembly_input)

    return assembler


if __name__ == "__main__":
    assembler = test_Au3_example()

