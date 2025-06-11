import argparse
from DARTassembler.src.modules.modules import DBInfo, Concat, Configs
from DARTassembler.src.metalig.ligandfilters import LigandFilters

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

def none_or_str(p: str):
    if p.lower() == "none":
        return None
    return p

def main():
#     desc = f"""DART command-line interface for assembling novel transition metal complexes from a database of ligands. Available modules: {", ".join(modules)}.
# Usage: dart <module> --path <path>
# """
    print(init_cli_output)
    parser = argparse.ArgumentParser(prog="DARTassembler")
    subparsers = parser.add_subparsers(dest="module", required=True, help="DART modules")

    # DBInfo command
    parser_dbinfo = subparsers.add_parser("dbinfo", help="Get info files from a ligand database")
    parser_dbinfo.add_argument("--path", type=str, default="metalig", help="Path to ligand database")
    parser_dbinfo.add_argument("--outdir", type=none_or_str, default='.', help="Output directory for csv and xyz file")
    parser_dbinfo.add_argument("--n-max-ligands", type=int, default=None, help="Maximum number of ligands to read in from the database.")
    parser_dbinfo.add_argument("--with-metal", action="store_true", default=True, help="Include metal atoms in the output")
    parser_dbinfo.set_defaults(func=lambda args: DBInfo().run_from_cli(path=args.path, outdir=args.outdir,
                                                                       n_max_ligands=args.n_max_ligands,
                                                                       with_metal=args.with_metal))

    # Concat command
    parser_concat = subparsers.add_parser("concat", help="Concatenate multiple ligand databases")
    parser_concat.add_argument("--paths", nargs="+", type=str, help="Paths to ligand databases")
    parser_concat.add_argument("--outpath", type=none_or_str, default="concat_ligand_db.jsonlines", help="Output path for concatenated database")
    parser_concat.add_argument("--n-max-ligands", type=int, default=None, help="Maximum number of ligands per database")
    parser_concat.set_defaults(func=lambda args: Concat().run_from_cli(paths=args.paths, outpath=args.outpath, n_max_ligands=args.n_max_ligands))

    # Configs command
    parser_configs = subparsers.add_parser("configs", help="Retrieve default yaml configuration files")
    parser_configs.add_argument("--outdir", type=none_or_str, default=None, help="Directory to save the configuration files")
    parser_configs.set_defaults(func=lambda args: Configs().run_from_cli(outdir=args.outdir))

    # Install cli
    args = parser.parse_args()
    output = args.func(args)

    return output

if __name__ == '__main__':
    output = main()