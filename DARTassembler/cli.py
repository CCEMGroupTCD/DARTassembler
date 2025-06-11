import argparse
from DARTassembler.src.modules.modules import DBInfo, Concat, Configs

init_cli_output = r"""================================================================================
                            ____  ___    ____  ______
                           / __ \/   |  / __ \/_  __/
                          / / / / /| | / /_/ / / /   
                         / /_/ / ___ |/ _, _/ / /    
                        /_____/_/  |_/_/ |_| /_/     
        
          DART - Directed Assembly of Random Transition metal complexes
              Developed by the CCEM group at Trinity College Dublin
================================================================================"""


modules = ['ligandfilters', 'assembler', 'dbinfo', 'concat', 'installtest', 'configs']

def check_n_args(args, n):
    if len(args) != n:
        raise ValueError(f'Expected {n} path, got {len(args)} arguments.')

def main():
#     desc = f"""DART command-line interface for assembling novel transition metal complexes from a database of ligands. Available modules: {", ".join(modules)}.
# Usage: dart <module> --path <path>
# """
    parser = argparse.ArgumentParser(prog="DARTassembler")
    subparsers = parser.add_subparsers(dest="module", required=True, help="DART modules")

    # DBInfo command
    parser_dbinfo = subparsers.add_parser("dbinfo", help="Get info files from a ligand database")
    parser_dbinfo.add_argument("--path", type=str, default="metalig", help="Path to ligand database")
    parser_dbinfo.add_argument("--outpath", type=str, default='.csv', help="Output path for csv file")
    parser_dbinfo.add_argument("--n-max-ligands", type=int, default=None, help="Maximum number of ligands to read in from the database.")
    parser_dbinfo.add_argument("--with-metal", action="store_true", default=True, help="Include metal atoms in the output")
    parser_dbinfo.set_defaults(func=lambda args: DBInfo().run_from_cli(path=args.path, outpath=args.outpath,
                                                                       n_max_ligands=args.n_max_ligands,
                                                                       with_metal=args.with_metal))

    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()
