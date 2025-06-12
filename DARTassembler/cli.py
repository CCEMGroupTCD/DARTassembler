import argparse
from DARTassembler.src.modules.modules import DBInfo, Concat, Configs
from DARTassembler.src.metalig.ligandfilters import LigandFilters
from DARTassembler.src.assembler.assembler import Assembler

init_cli_output = r"""================================================================================
                            ____  ___    ____  ______
                           / __ \/   |  / __ \/_  __/
                          / / / / /| | / /_/ / / /   
                         / /_/ / ___ |/ _, _/ / /    
                        /_____/_/  |_/_/ |_| /_/     
        
          DART - Directed Assembly of Random Transition metal complexes
              Developed by the CCEM group at Trinity College Dublin
================================================================================"""


def none_or_str(p: str):
    if p.lower() == "none":
        return None
    return p

def main():
    print(init_cli_output)
    parser = argparse.ArgumentParser(prog="DARTassembler", description='DART command-line interface for assembling novel transition metal complexes from a database of ligands.')
    subparsers = parser.add_subparsers(dest="module", required=True, help="DART modules")

    # Assembler
    parser_assembler = subparsers.add_parser("assembler", help="Assemble complexes from ligands in a high-throughput manner")
    parser_assembler.add_argument("--input", type=str, default=None, help="Path to a yml file with assembler options")
    parser_assembler.set_defaults(func=lambda args: Assembler.run_from_cli(input=args.input))

    # LigandFilters
    parser_ligandfilters = subparsers.add_parser("ligandfilters", help="Filter ligands from a database")
    parser_ligandfilters.add_argument("--input", type=str, default=None, help="Path to a yml file with ligand filter options")
    parser_ligandfilters.set_defaults(func=lambda args: LigandFilters.run_from_cli(input=args.input))

    # DBInfo
    parser_dbinfo = subparsers.add_parser("dbinfo", help="Get info files from a ligand database")
    parser_dbinfo.add_argument("--db", type=str, default="metalig", help="Path to ligand database")
    parser_dbinfo.add_argument("--outdir", type=none_or_str, default='.', help="Output directory for csv and xyz file")
    parser_dbinfo.add_argument("--n", type=int, default=None, help="Maximum number of ligands to read in from the database.")
    parser_dbinfo.add_argument("--metal", action="store_true", default=True, help="Include metal atoms in the output")
    parser_dbinfo.set_defaults(func=lambda args: DBInfo.run_from_cli(db=args.db, outdir=args.outdir, n=args.n, metal=args.metal))

    # Concat
    parser_concat = subparsers.add_parser("concat", help="Concatenate multiple ligand databases")
    parser_concat.add_argument("--dbs", nargs="+", type=str, help="Paths to ligand databases")
    parser_concat.add_argument("--outpath", type=none_or_str, default="concat_ligand_db.jsonlines", help="Output path for concatenated database")
    parser_concat.add_argument("--n", type=int, default=None, help="Maximum number of ligands per database")
    parser_concat.set_defaults(func=lambda args: Concat.run_from_cli(dbs=args.dbs, outpath=args.outpath, n=args.n))

    # Configs
    parser_configs = subparsers.add_parser("configs", help="Retrieve default yaml configuration files")
    parser_configs.add_argument("--outdir", type=none_or_str, default='.', help="Directory to save the configuration files")
    parser_configs.set_defaults(func=lambda args: Configs.run_from_cli(outdir=args.outdir))

    # Install cli
    args = parser.parse_args()
    output = args.func(args)

    return output

if __name__ == '__main__':
    output = main()