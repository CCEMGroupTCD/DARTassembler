import argparse
from DARTassembler import ligandfilters, assembler, dbinfo, concat, installtest, configs

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
    desc = f"""DART command-line interface for assembling novel transition metal complexes from a database of ligands. Available modules: {", ".join(modules)}.
Usage: dart <module> --path <path>
"""
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('module', choices=modules, help='DART module that you want to use')
    parser.add_argument('--path', required=True, help='Path to the input file(s)', nargs='*')

    args = parser.parse_args()

    print(init_cli_output)
    module_title = f'  Execute DART module: {args.module}  '
    print(f'{module_title:#^80}')
    if args.path:
        print(f'Input path: {args.path[0] if len(args.path) == 1 else args.path}')

    if args.module == 'ligandfilters':
        check_n_args(args.path, 1)
        ligandfilters(args.path[0])
    elif args.module == 'assembler':
        check_n_args(args.path, 1)
        assembler(args.path[0])
    elif args.module == 'dbinfo':
        check_n_args(args.path, 1)
        dbinfo(args.path[0])
    elif args.module == 'concat':
        concat(args.path)
    elif args.module == 'installtest':
        if len(args.path) == 0:
            path = None
        else:
            check_n_args(args.path, 1)
            path = args.path[0]
        installtest(path)
    elif args.module == 'configs':
        check_n_args(args.path, 1)
        configs(args.path[0])
    else:
        raise ValueError(f'Unknown module {args.module}.')

if __name__ == '__main__':
    main()
