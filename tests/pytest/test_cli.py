"""
Integration test for the cli, running all the modules in the DARTassembler package.
"""
import os
from pathlib import Path
from DARTassembler.src.constants.paths import default_assembler_yml_path, default_ligandfilters_yml_path
from DARTassembler.src.misc import cli
import sys
from DARTassembler.src.misc.tests import integration_test

def run_cli(argv):
    old_argv = sys.argv # Save the original sys.argv to restore it later
    try:
        sys.argv = ["DARTassembler"] + argv
        return cli.main()
    finally:
        sys.argv = old_argv

@integration_test(name='cli')
def test_cli(outdir: Path) -> None:
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
            'n': 1000
        },
        'assembler': {
            'input': default_assembler_yml_path,
            'n_max_ligands': 1000
        },
    }
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

    return


if __name__ == "__main__":
    test_cli()