"""
Integration test for the cli, running all the modules in the DARTassembler package.
"""
import os
from pathlib import Path
from DARTassembler.src.constants.paths import project_path, default_assembler_yml_path, default_ligandfilters_yml_path
from DARTassembler.src.misc import cli
import sys
from shutil import rmtree

def run_cli(argv):
    old_argv = sys.argv # Save the original sys.argv to restore it later
    try:
        sys.argv = ["DARTassembler"] + argv
        return cli.main()
    finally:
        sys.argv = old_argv

def test_cli():
    outdir = project_path().extend('tests', 'integration_tests', 'cli', 'data_output')
    # Remove the output directory if it exists to start fresh
    if outdir.exists():
        rmtree(outdir)
    arguments = {
        'dbinfo': {
            'n': 100,
        },
        'concat': {
            'dbs': 'metalig',
            'n': 100,
            'outpath': 'concat_ligand_db.jsonlines',
        },
        'configs': {
        },
        'ligandfilters': {
            'input': default_ligandfilters_yml_path,
        },
        'assembler': {
            'input': default_assembler_yml_path,
        },
    }
    old_cwd = Path.cwd()    # Save the current working directory to return to it later
    try:
        for module, args in arguments.items():
            module_outdir = outdir / module
            module_outdir.mkdir(parents=True, exist_ok=True)
            os.chdir(module_outdir)  # Change to the output directory for the test
            print(f"Testing {module} module with arguments: {args}")
            # Convert dict to command line arguments
            args_list = [module]
            for k, v in args.items():
                if isinstance(v, list):
                    v = [str(item) for item in v]
                else:
                    v = str(v)
                args_list.extend(['--'+k, v])
            run_cli(argv=args_list)
    finally:
        # Change back to the original working directory
        os.chdir(old_cwd)

    #%% ==============    Doublecheck refactoring    ==================
    from DARTassembler.src.misc.tests import IntegrationTest
    old_dir = outdir.parent / 'benchmark_data_output'
    if old_dir.exists():
        test = IntegrationTest(new_dir=outdir, old_dir=old_dir)
        test.compare_all()
        all_modules = list(arguments.keys())
        print(f'Test for cli for modules {", ".join(all_modules)} passed!')
    else:
        print(f'ATTENTION: could not find benchmark folder "{old_dir}"!')

    return


if __name__ == "__main__":
    assembly = test_cli()