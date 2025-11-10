"""
This file contains the base module class and several small modules for the DARTassembler package. The :ref:`assembler <assembler>` and :ref:`ligandfilters <ligandfilters>` modules are not included here.
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
    def __init_subclass__(cls):
        super().__init_subclass__()
        # Set class attributes on each subclass immediately when it’s defined so you can access them without instantiating the class (as a cls object).
        cls._module_name = cls.__name__.lower()
        cls._desc        = (cls.__doc__ or "").strip()

    def __init__(self):
        pass

    @classmethod
    def _before_run_from_cli(cls) -> None:
        """
        Base method for running the module.
        """
        title = f'     {cls._module_name.upper()} MODULE    '
        cls._print(f'{title:=^80}')
        cls._print(f'{cls._module_name}: {cls._desc}')

    @classmethod
    def _after_run_from_cli(cls) -> None:
        cls._print(f"Done! Exiting {cls._module_name} module.")

    @classmethod
    def _print_cli_input(cls, **kwargs) -> None:
        """
        Print the input parameters for the module.
        """
        cls._print(f"Input parameters:")
        for key, value in kwargs.items():
            cls._print(f"  - {key}: {value}")
        cls._print(f'Starting {cls._module_name} module with the above parameters...')

    @classmethod
    def run_from_cli(cls, **kwargs):
        cls._before_run_from_cli()
        cls._print_cli_input(**kwargs)
        module = cls()  # Create an instance of the module class
        module.run(**kwargs)
        module._after_run_from_cli()

        return module

    @staticmethod
    def _print(text: str) -> None:
        """
        Print text to the console.
        """
        print(textwrap.fill(text=text, width=80))


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

    def run(self, dbs: list[Union[str,Path]], outpath: Union[str,Path,None]=None, n: Union[int, None] = None) -> LigandDB:
        """
        Concatenate multiple ligand databases into one.

        :param dbs: Paths to the ligand databases.
        :param outpath: Path to the output ligand database. If None, no output file is saved.
        :param n: Maximum number of ligands to be read in from each ligand database. If None, all ligands are read in. This is useful for testing purposes.
        """
        dbs = [get_correct_ligand_db_path_from_input(db) for db in dbs] # Ensure paths are correct

        # Load all ligand databases
        ligand_dbs = [LigandDB.from_json(db, n_max=n) for db in dbs]

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
            outpath = Path(outpath).resolve()
            full_db.save_to_file(outpath)
            print(f"Saved concatenated ligand database to `{outpath.name}`.")

        return full_db


class DBInfo(BaseModule):
    """
    This module reads in a ligand database (.jsonlines) and saves a .csv file with an overview of the ligands and a concatenated .xyz file with the structures of the ligands.
    """
    def __init__(self) -> None:
        super().__init__()

    def run(self, db: Union[str, Path,None]='metalig', outdir: Union[str, Path, None] = None, n: Union[int, None] = None, metal: bool=True) -> tuple[LigandDB, pd.DataFrame, str]:
        """
        Reads in the given ligand database and saves a .csv file and a concatenated .xyz file with an overview of the ligands.

        :param db: Path to the ligand database. The default path is 'metalig', which points to the full ligand database.
        :param outdir: Path to the output .csv file. If None, no output file is saved. If '.csv', the output file is saved in the same directory as the input file with the same name but with the .csv extension.
        :param n: Maximum number of ligands to be read in from the initial full ligand database. If None, all ligands are read in. This is useful for testing purposes.
        :param metal: If True, the metal atom is included in the concatenated .xyz file. If False, only the ligand is included.
        :return: Tuple of (LigandDB, DataFrame, concatenated xyz string) of the ligands.
        """
        db = get_correct_ligand_db_path_from_input(db)
        ligand_db = LigandDB.from_json(db, n_max=n)

        if outdir is None:
            df = ligand_db.get_df()
            xyz_string = ligand_db.get_concat_xyz_string(with_metal=metal)
        else:
            # Handle default output path
            if outdir == '.':
                outdir = Path.cwd()
            outdir = Path(outdir).resolve()
            outdir.parent.mkdir(parents=True, exist_ok=True)

            # Save to csv
            stem = db.name.removesuffix('.bz2').removesuffix('.jsonlines').removesuffix('.json')
            print(f'Saving ligand info and structures to `{outdir.name}`...')
            csv_filename = stem + '.csv'
            df = ligand_db.save_to_csv(Path(outdir, csv_filename))
            print(f'  - Saved .csv to `{csv_filename}`.')
            # Save to concatenated xyz file
            xyz_filename = 'concat_' + stem + '.xyz'
            xyz_string = ligand_db.save_to_concat_xyz(Path(outdir, xyz_filename), with_metal=metal)
            print(f'  - Saved .xyz to `{xyz_filename}`.')

        return ligand_db, df, xyz_string


class Configs(BaseModule):
    """
    This module returns the default .yml configuration files for the assembler and the ligandfilters. Use these files as templates to adapt to your specific needs.
    """
    def __init__(self) -> None:
        super().__init__()

    def run(self, outdir: Union[str, Path,None]=None) -> tuple[dict, dict]:
        """
        Get the default yaml configuration files for the assembler and the ligandfilters and optionally save them to the specified output path.

        :param outdir: Output directory where the configuration files will be saved. If None, the files are not saved and only the dictionaries are returned.
        :return: A tuple containing the assembler options and the ligandfilters options as dictionaries.
        """
        # Read yaml files for output
        assembler_options = read_yaml(default_assembler_yml_path)
        ligandfilters_options = read_yaml(default_ligandfilters_yml_path)

        if outdir is not None:
            outdir = Path(outdir).resolve()
            outdir.mkdir(parents=True, exist_ok=True)

            # Copy assembler.yml
            filename = default_assembler_yml_path.name
            print(f'\t- get {filename}')
            dest = Path(outdir, filename)
            copyfile(default_assembler_yml_path, dest)

            # Copy ligandfilters.yml
            filename = default_ligandfilters_yml_path.name
            print(f'\t- get {filename}')
            dest = Path(outdir, filename)
            copyfile(default_ligandfilters_yml_path, dest)

            print(f"Saved config files to `{outdir.name}`.")

        return assembler_options, ligandfilters_options


if __name__ == "__main__":
    n_max = 100  # Set a maximum number of ligands for testing purposes

    # Try out the modules without saving any output files.
    assembler_dict, ligandfilters_dict = Configs.run_from_cli(outdir=None)
    ligand_db, df_ligands, xyz_ligands = DBInfo.run_from_cli(n=n_max, outdir=None)
    out_db = Concat.run_from_cli(dbs=['metalig'], outpath=None, n=n_max)
