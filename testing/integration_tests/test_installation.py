"""
Integration test for the assembly of complexes.
"""
from DARTassembler import installtest
from pathlib import Path
import shutil

def test_assemble_complexes():
    outdir = Path('installation/data_output')

    # Delete output directory so that the test detects if files are not written.
    shutil.rmtree(outdir, ignore_errors=True)

    installtest(outdir=outdir)

    #%% ==============    Doublecheck refactoring    ==================
    from dev.test.Integration_Test import IntegrationTest
    old_dir = Path(outdir.parent, 'benchmark_data_output')
    if old_dir.exists():
        test = IntegrationTest(new_dir=outdir, old_dir=old_dir)
        test.compare_all()
        print('Test for installation passed!')
    else:
        print(f'ATTENTION: could not find benchmark folder "{old_dir}"!')

    return


if __name__ == "__main__":
    test_assemble_complexes()