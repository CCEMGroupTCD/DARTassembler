"""
This module concatenates multiple ligand databases into one.
"""
from pathlib import Path
from shutil import copyfile
from typing import Union
import textwrap
import pandas as pd
from DARTassembler.src.constants.paths import default_assembler_yml_path, default_ligandfilters_yml_path
from DARTassembler.src.metalig.db import LigandDB
from DARTassembler.src.misc.io import get_correct_ligand_db_path_from_input, read_yaml


class BaseModule(object):
    """
    Base class for all modules in the DARTassembler package. Implements the basic structure for running modules from the command line interface (CLI).
    """
    def __init__(self):
        self.module_name = self.__class__.__name__.lower() # Name of the module, e.g. 'concat', 'ligandfilters', etc.
        self.desc = self.__doc__.strip().strip("\n")

    def _before_run_from_cli(self) -> None:
        """Base method for running the module."""
        title = f'     {self.module_name.upper()} MODULE    '
        print(f'{title:=^80}')
        print('\n'.join(textwrap.wrap(self.desc, 80)))

    def _after_run_from_cli(self) -> None:
        print(f"Done! Exiting {self.module_name} module.")

    def _print_cli_input(self, **kwargs) -> None:
        """
        Print the input parameters for the module.
        """
        print(f"Input options:")
        for key, value in kwargs.items():
            print(f"  - {key}: {value}")
        print(f'Starting {self.module_name} module with the above parameters...')

    def _run_from_cli(self, **kwargs):
        self._before_run_from_cli()
        self._print_cli_input(**kwargs)
        output = self.run(**kwargs)
        self._after_run_from_cli()

        return output

    def run(self, *args, **kwargs):
        """
        This method should be implemented in the subclass.
        """
        raise NotImplementedError("This method should be implemented in each subclass.")


class Concat(BaseModule):
    """
    This module concatenates multiple ligand databases into one.
    """
    def __init__(self) -> None:
        super().__init__()

    def run_from_cli(self, paths: list[Union[str,Path]], outpath: Union[str,Path,None]='concat_ligand_db.jsonlines', n_max_ligands: Union[int, None] = None) -> LigandDB:
        """
        Run the concatenation of ligand databases from the command line interface.
        :param paths: Paths to the ligand databases.
        :param outpath: Path to the output ligand database. If None, no output file is saved.
        """
        return super()._run_from_cli(paths=paths, outpath=outpath, n_max_ligands=n_max_ligands)

    def run(self, paths: list[Union[str,Path]], outpath: Union[str,Path,None]=None, n_max_ligands: Union[int, None] = None) -> LigandDB:
        """
        Concatenate multiple ligand databases into one.
        :param paths: Paths to the ligand databases.
        :param outpath: Path to the output ligand database. If None, no output file is saved.
        :param n_max_ligands: Maximum number of ligands to be read in from each ligand database. If None, all ligands are read in. This is useful for testing purposes.
        """
        paths = [get_correct_ligand_db_path_from_input(path) for path in paths] # Ensure paths are correct

        # Load all ligand databases
        ligand_dbs = [LigandDB.from_json(path, n_max=n_max_ligands) for path in paths]

        # Print number of ligands in each database
        for i, db in enumerate(ligand_dbs):
            print(f"Ligand database {i + 1} contains {len(db.db)} unique ligands.")

        # Concatenate ligand databases
        full_db = {}
        for db in ligand_dbs:
            full_db.update(db.db)
        full_db = LigandDB(full_db)

        print(f"The final concatenated ligand database contains {len(full_db.db)} unique ligands.")

        # Save concatenated ligand database
        if outpath is not None:
            full_db.save_to_file(outpath)
            print(f"Saved concatenated ligand database to `{outpath}`.")

        return full_db


class DBInfo(BaseModule):
    """
    This module reads in a ligand database (.jsonlines) and saves a .csv file with an overview of the ligands and a concatenated .xyz file with the structures of the ligands.
    """
    def __init__(self) -> None:
        super().__init__()

    def run_from_cli(self, path: Union[str, Path,None]='metalig', outpath: Union[str, Path, None] = '.csv', n_max_ligands: Union[int, None] = None, with_metal: bool=True):
        """
        Run the DBInfo module from the command line interface.
        """
        return super()._run_from_cli(
            path=path,
            outpath=outpath,
            n_max_ligands=n_max_ligands,
            with_metal=with_metal
        )

    def run(self, path: Union[str, Path,None]='metalig', outpath: Union[str, Path, None] = None, n_max_ligands: Union[int, None] = None, with_metal: bool=True) -> tuple[LigandDB, pd.DataFrame, str]:
        """
        Reads in the given ligand database and saves a .csv file and a concatenated .xyz file with an overview of the ligands.
        :param path: Path to the ligand database. The default path is 'metalig', which points to the full ligand database.
        :param outpath: Path to the output .csv file. If None, the output file will be saved in the same directory as the input file with the same name as the input file but with the .csv extension.
        :param n_max_ligands: Maximum number of ligands to be read in from the initial full ligand database. If None, all ligands are read in. This is useful for testing purposes.
        :param with_metal: If True, the metal atom is included in the concatenated .xyz file. If False, only the ligand is included.
        :return: Tuple of (LigandDB, DataFrame, concatenated xyz string) of the ligands.
        """
        path = get_correct_ligand_db_path_from_input(path)
        db = LigandDB.from_json(path, n_max=n_max_ligands)

        print('Saving ligand info and structures...')

        if outpath == '.csv':
            outpath = path.with_suffix('.csv')
        if outpath is not None:
            outpath = Path(outpath)
            outpath.parent.mkdir(parents=True, exist_ok=True)

        # Save to csv
        if outpath is None:
            df = db.get_df()
        else:
            df = db.save_to_csv(outpath)
            print(f'  - Saved .csv to `{outpath.name}`.')

        # Save to concatenated xyz file
        if outpath is None:
            xyz_string = db.get_concat_xyz_string(with_metal=with_metal)
        else:
            xyz_filename = str(Path(f'concat_{outpath.with_suffix("").name}.xyz'))
            xyz_outpath = outpath.parent.joinpath(xyz_filename)
            xyz_string = db.save_to_concat_xyz(xyz_outpath, with_metal=with_metal)
            print(f'  - Saved .xyz to `{xyz_filename}`.')

        return db, df, xyz_string

class Configs(BaseModule):
    """
    This module returns the default .yml configuration files for the assembler and the ligandfilters. Use these files as templates to adapt to your specific needs.
    """
    def __init__(self) -> None:
        super().__init__()

    def run_from_cli(self, outdir: Union[str, Path]=None) -> tuple[dict, dict]:
        """
        Run the Configs module from the command line interface.
        """
        return super()._run_from_cli(outdir=outdir)

    def run(self, outdir: Union[str, Path,None]=None) -> tuple[dict, dict]:
        """
        Get the default yaml configuration files for the assembler and the ligandfilters and optionally save them to the specified output path.
        :param outdir: Output directory where the configuration files will be saved. If None, the files are not saved and only the dictionaries are returned.
        :return: A tuple containing the assembler options and the ligandfilters options as dictionaries.
        """
        if outdir is not None:
            outdir = Path(outdir)
            outdir.mkdir(parents=True, exist_ok=True)

        # Copy assembler.yml
        filename = default_assembler_yml_path.name
        print(f'\t- get {filename}')
        assembler_options = read_yaml(default_assembler_yml_path)
        if outdir is not None:
            dest = Path(outdir, filename)
            copyfile(default_assembler_yml_path, dest)

        # Copy ligandfilters.yml
        filename = default_ligandfilters_yml_path.name
        print(f'\t- get {filename}')
        ligandfilters_options = read_yaml(default_ligandfilters_yml_path)
        if outdir is not None:
            dest = Path(outdir, filename)
            copyfile(default_ligandfilters_yml_path, dest)

        if outdir is not None:
            print(f"Saved config files to `{outdir}`.")

        return assembler_options, ligandfilters_options

if __name__ == "__main__":
    n_max = 100  # Set a maximum number of ligands for testing purposes

    # Try out the modules.
    assembler_dict, ligandfilters_dict = Configs().run_from_cli()
    ligand_db, df_ligands, xyz_ligands = DBInfo().run_from_cli(n_max_ligands=n_max, outpath=None)
    out_db = Concat().run_from_cli(paths=['metalig'], outpath=None, n_max_ligands=n_max)
