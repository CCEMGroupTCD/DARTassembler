import argparse
from DARTassembler import filter_ligands, assemble_complexes, make_ligand_db_csv


def main():
    modules = ['ligandfilters', 'assembler', 'dbinfo']
    desc = f"""DART command-line interface for assembling novel transition metal complexes from a database of ligands. Available modules: {", ".join(modules)}.
Usage: dart <module> --path <path>
Examples: 
    -> dart assembler --path assembly_input.yml
    -> dart ligandfilters --path ligandfilters_input.yml
    -> dart dbinfo --path ligand_db.jsonlines
"""
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('module', choices=modules, help='DART module that you want to use')
    parser.add_argument('--path', required=True, help='Path to the input file')

    args = parser.parse_args()

    print(f'Executing module {args.module} with input file {args.path}')
    if args.module == 'ligandfilters':
        filter_ligands(args.path)
    elif args.module == 'assembler':
        assemble_complexes(args.path)
    elif args.module == 'dbinfo':
        make_ligand_db_csv(args.path)
    else:
        raise ValueError(f'Unknown module {args.module}.')

if __name__ == '__main__':
    main()
