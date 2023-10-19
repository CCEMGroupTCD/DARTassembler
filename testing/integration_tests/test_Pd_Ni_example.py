import DARTassembler
from pathlib import Path
from DARTassembler.src.constants.Paths import project_path

def test_Pd_Ni_example(nmax=False):
    indir = project_path().extend('testing/integration_tests/Pd_Ni_example_assembly/data_input')
    Br_filter = DARTassembler.filter_ligands(Path(indir, 'ligandfilters_Br.yml'), nmax=100)
    phenyl_filter = DARTassembler.filter_ligands(Path(indir, 'ligandfilters_phenyl.yml'), nmax=100)
    P_N_filter = DARTassembler.filter_ligands(Path(indir, 'ligandfilters_P_N_ligands.yml'), nmax=nmax)
    assembly = DARTassembler.assemble_complexes(Path(indir, 'Pd_Ni_assembly_input.yml'))

    #%% ==============    Doublecheck refactoring    ==================
    from dev.test.Integration_Test import IntegrationTest
    for benchmark in ['ligand_databases', 'assembled_complexes']:
        print(f'Comparing {benchmark}...')
        old_dir = Path(assembly.output_path.parent, f'benchmark_{benchmark}')
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
    Br_filter, phenyl_filter, P_N_filter, assembly = test_Pd_Ni_example(nmax=nmax)



