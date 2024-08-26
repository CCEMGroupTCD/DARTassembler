import os

import DARTassembler
from pathlib import Path
from DARTassembler.src.constants.Paths import project_path
from distutils.dir_util import copy_tree


def copy_input_files_from_example_and_adapt_paths(example_input_dir, test_input_dir):
    copy_tree(str(example_input_dir), str(test_input_dir))

    replace_paths = {
        'ligand_db_Br.jsonlines': 'Pd_Ni_example_assembly/ligand_databases/Br/ligand_db_Br.jsonlines',
        'ligand_db_phenyl.jsonlines': 'Pd_Ni_example_assembly/ligand_databases/phenyl/ligand_db_phenyl.jsonlines',
        'ligand_db_P_N_donors.jsonlines': 'Pd_Ni_example_assembly/ligand_databases/P_N_ligands/ligand_db_P_N_donors.jsonlines',
        'assembler_output': 'Pd_Ni_example_assembly/assembled_complexes',
        'input/Pd_phenyl_geometry_modification.xyz': 'Pd_Ni_example_assembly/data_input/Pd_phenyl_geometry_modification.xyz',
        'input/Ni_phenyl_geometry_modification.xyz': 'Pd_Ni_example_assembly/data_input/Ni_phenyl_geometry_modification.xyz'
    }

    should_adapt_n_paths = {
        'Ni_phenyl_geometry_modification.xyz': 0,
        'Pd_Ni_assembly_input.yml': 6,
        'Pd_phenyl_geometry_modification.xyz': 0,
        'README': 0,
        'ligandfilters_Br.yml': 1,
        'ligandfilters_P_N_ligands.yml': 1,
        'ligandfilters_phenyl.yml': 1,
    }

    for file in sorted(os.listdir(test_input_dir)):
        n_adapted_paths = 0
        with open(Path(test_input_dir, file), 'r') as f:
            text = f.read()
        for old, new in replace_paths.items():
            if old in text:
                text = text.replace(old, new)
                n_adapted_paths += 1
        with open(Path(test_input_dir, file), 'w') as f:
            f.write(text)

        # Check if the correct number of paths was adapted
        filename = Path(test_input_dir, file).name
        should_n_paths = should_adapt_n_paths.get(filename, None)
        if should_n_paths is not None:
            assert n_adapted_paths == should_n_paths, f'Adapted {n_adapted_paths} paths in {filename}, but should have adapted {should_n_paths} paths!'
        else:
            raise ValueError(f'No information on how many paths should be adapted for {filename}!')

    return

def test_Pd_Ni_example(nmax=False, skip_filters=False):
    example_input_dir = project_path().extend('examples/Pd_Ni_Cross_Coupling/generate_complexes/input')
    indir = Path('Pd_Ni_example_assembly/data_input')
    benchmark_dir = project_path().extend('testing/integration_tests/Pd_Ni_example_assembly')

    copy_input_files_from_example_and_adapt_paths(example_input_dir, indir)

    if not skip_filters:
        Br_filter = DARTassembler.ligandfilters(Path(indir, 'ligandfilters_Br.yml'), nmax=100)
        phenyl_filter = DARTassembler.ligandfilters(Path(indir, 'ligandfilters_phenyl.yml'), nmax=100)
        P_N_filter = DARTassembler.ligandfilters(Path(indir, 'ligandfilters_P_N_ligands.yml'), nmax=nmax)
    else:
        Br_filter = phenyl_filter = P_N_filter = None
    assembly = DARTassembler.assembler(Path(indir, 'Pd_Ni_assembly_input.yml'))

    #%% ==============    Doublecheck refactoring    ==================
    from dev.test.Integration_Test import IntegrationTest
    benchmarks = ['ligand_databases', 'assembled_complexes'] if not skip_filters else ['assembled_complexes']
    for benchmark in benchmarks:
        print(f'Comparing {benchmark}...')
        old_dir = Path(benchmark_dir, f'benchmark_{benchmark}')
        new_dir = Path(assembly.output_path.parent, benchmark)
        if old_dir.exists():
            test = IntegrationTest(new_dir=new_dir, old_dir=old_dir)
            test.compare_all()
            print(f'Test {benchmark} passed!')
        else:
            print(f'ATTENTION: could not find benchmark folder "{old_dir}"!')

    return Br_filter, phenyl_filter, P_N_filter, assembly



if __name__ == '__main__':

    nmax = False
    skip_filters = False

    Br_filter, phenyl_filter, P_N_filter, assembly = test_Pd_Ni_example(nmax=nmax, skip_filters=skip_filters)



