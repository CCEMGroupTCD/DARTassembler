import argparse
from DARTassembler import filter_ligands, assemble_complexes, save_dbinfo, concatenate_ligand_databases, run_installation_test

modules = ['ligandfilters', 'assembler', 'dbinfo', 'concat', 'test']

def check_n_args(args, n):
    if len(args) != n:
        raise ValueError(f'Expected {n} path, got {len(args)} arguments.')

def main():
    desc = f"""DART command-line interface for assembling novel transition metal complexes from a database of ligands. Available modules: {", ".join(modules)}.
Usage: dart <module> --path <path>
Examples: 
    -> dart assembler --path assembly_input.yml
    -> dart ligandfilters --path ligandfilters_input.yml
    -> dart dbinfo --path ligand_db.jsonlines
    -> dart concat --path ligand_db1.jsonlines ligand_db2.jsonlines
    -> dart test
"""
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('module', choices=modules, help='DART module that you want to use')
    parser.add_argument('--path', required=True, help='Path to the input file(s)', nargs='*')

    args = parser.parse_args()

    print(f'Executing module {args.module} with input file {args.path}')
    if args.module == 'ligandfilters':
        check_n_args(args.path, 1)
        filter_ligands(args.path[0])
    elif args.module == 'assembler':
        check_n_args(args.path, 1)
        assemble_complexes(args.path[0])
    elif args.module == 'dbinfo':
        check_n_args(args.path, 1)
        save_dbinfo(args.path[0])
    elif args.module == 'concat':
        concatenate_ligand_databases(args.path)
    elif args.module == 'test':
        check_n_args(args.path, 1)
        run_installation_test(args.path[0])
    else:
        raise ValueError(f'Unknown module {args.module}.')

if __name__ == '__main__':
    main()
