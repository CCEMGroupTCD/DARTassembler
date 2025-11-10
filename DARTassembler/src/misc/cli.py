"""
DART implements a list of commands (or modules) for the terminal which you can use to interact with the package. For the :ref:`Assembler Module <assembler>` and the :ref:`Ligand Filters Module <ligandfilters>`, you will typically want to provide a configuration file in YAML format that specifies the options for the module. You can generate template configuration files using the ``configs`` module. You can use the ``dbinfo`` to get more information on the ligands in a ligand database file or the ``concat`` to concatenate multiple ligand database files into one.

.. argparse::
   :module: DARTassembler.src.misc.cli
   :func: get_parser
   :prog: DARTassembler
"""
import argparse
from DARTassembler import Assembler, LigandFilters, DBInfo, Concat, Configs
from importlib.metadata import version

dart_version = version("DARTassembler")
init_cli_output = rf"""================================================================================
                            ____  ___    ____  ______
                           / __ \/   |  / __ \/_  __/
                          / / / / /| | / /_/ / / /   
                         / /_/ / ___ |/ _, _/ / /    
                        /_____/_/  |_/_/ |_| /_/     
        
          DART - Directed Assembly of Random Transition metal complexes
                              Version: {dart_version}
    Developed by the CCEM group at Trinity College Dublin & CIC energiGUNE
================================================================================"""

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ["true"]:
        return True
    elif v.lower() in ["false"]:
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")

def none_or_str(p: str):
    if p.lower() == "none":
        return None
    return p

def get_parser():
    parser = argparse.ArgumentParser(prog="DARTassembler")
    subparsers = parser.add_subparsers(dest="module", required=True, help="DART modules to choose from. See each module's help for more information.")

    # Assembler
    parser_assembler = subparsers.add_parser("assembler", help="""
Execute the assembler module using the provided configuration file.

Example: ``DARTassembler assembler --input assembler.yml``
""")
    parser_assembler.add_argument("--input", type=str, required=True, help="Path to a .yml file with assembler options.")
    parser_assembler.add_argument("--n_max_ligands", type=int, default=None, help="Maximum number of ligands to read in from the ligand databases. Useful for testing purposes. If None, all ligands are read in (default). When specified here, this argument overrides the value in the .yml file.")
    parser_assembler.set_defaults(func=lambda args: Assembler.run_from_cli(input=args.input, n_max_ligands=args.n_max_ligands))

    # LigandFilters
    parser_ligandfilters = subparsers.add_parser("ligandfilters", help="""
Execute the ligandfilters module using the provided configuration file.

Example: ``DARTassembler ligandfilters --input ligandfilters.yml``
""")
    parser_ligandfilters.add_argument("--input", type=str, required=True, help="Path to a .yml file with ligand filter options.")
    parser_ligandfilters.add_argument("--n", type=int, default=None, help="Maximum number of ligands to read in from the database. Useful for testing purposes. If None, all ligands are read in (default). When specified here, this argument overrides the value in the .yml file.")
    parser_ligandfilters.set_defaults(func=lambda args: LigandFilters.run_from_cli(input=args.input, n=args.n))

    # DBInfo
    parser_dbinfo = subparsers.add_parser("dbinfo", help="""
Get info files (.csv & .xyz) from a .jsonlines ligand database. This can be used to browse the ligands in a database. The .csv file contains information about each ligand in the database such as stoichiometry, donors or formal charge. The concatenated .xyz file contains the 3D structures of all ligands in the database. To browse through the structures, you can use ``ase gui concat_LIGAND_DB.xyz``.

Example: ``DARTassembler dbinfo --db ligand_db.jsonlines --outdir . --n 1000 --metal False``
    """)
    parser_dbinfo.add_argument("--db", type=str, default="metalig", help="Path to a .jsonlines ligand database such as ``ligand_db.jsonlines``. Use the string ``'metalig'`` instead for the entire MetaLig database (default).")
    parser_dbinfo.add_argument("--outdir", type=none_or_str, default='.', help="Output directory for the .csv and .xyz files. By default, the current directory is used.")
    parser_dbinfo.add_argument("--n", type=int, default=None, help="Maximum number of ligands to read in from the database. Useful for testing purposes. If None, all ligands are read in (default).")
    parser_dbinfo.add_argument("--metal", type=str2bool, nargs='?', const=True, default=True, help="If ``True`` (default), a Hg pseudo metal center is included in the concatenated .xyz file. If ``False``, only the ligand is included.")
    parser_dbinfo.set_defaults(func=lambda args: DBInfo.run_from_cli(db=args.db, outdir=args.outdir, n=args.n, metal=args.metal))

    # Concat
    parser_concat = subparsers.add_parser("concat", help="""
Concatenate multiple ligand databases into one.

Example: ``DARTassembler concat --dbs ligand_db1.jsonlines ligand_db2.jsonlines --outpath concat_ligand_db.jsonlines``
""")
    parser_concat.add_argument("--dbs", nargs="+", required=True, type=str, help="List of paths to .jsonlines ligand databases to concatenate.")
    parser_concat.add_argument("--outpath", type=none_or_str, required=True, help="Output filepath for the concatenated database .jsonlines file, e.g., ``concat_ligand_db.jsonlines``.")
    parser_concat.add_argument("--n", type=int, default=None, help="Maximum number of ligands to read in from each input ligand database. Useful for testing purposes. If None, all ligands are read in (default).")
    parser_concat.set_defaults(func=lambda args: Concat.run_from_cli(dbs=args.dbs, outpath=args.outpath, n=args.n))

    # Configs
    parser_configs = subparsers.add_parser("configs", help="""
Retrieve two .yml configuration files as useful templates for the assembler and ligandfilters module.

Example: ``DARTassembler configs --outdir .``
""")
    parser_configs.add_argument("--outdir", type=none_or_str, default='.', help="Directory to save the configuration files to. By default, the current directory is used.")
    parser_configs.set_defaults(func=lambda args: Configs.run_from_cli(outdir=args.outdir))

    return parser

def main():
    print(init_cli_output)
    parser = get_parser()

    # Install cli
    args = parser.parse_args()
    args.func(args)

    return

if __name__ == '__main__':
    main()