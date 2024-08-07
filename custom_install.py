# from pathlib import Path
# import zipfile
# import shutil
# import os
# from DARTassembler.src.constants.Paths import default_ligand_db_path
# from setuptools import Command
#
# class CustomInstallCommand(Command):
#     description = "Unzip MetaLig .jsonlines.zip file after installation"
#     user_options = []
#
#     def initialize_options(self):
#         pass
#
#     def finalize_options(self):
#         pass
#
#     def run(self):
#         check_if_MetaLig_exists_else_uncompress_from_zip(delete_zip=True)
#
# def check_if_MetaLig_exists_else_uncompress_from_zip(delete_zip=False):
#     """
#     Checks if the MetaLig database exists as uncompressed file, and if not uncompresses it.
#     """
#     install.run(self)
#     if not Path(default_ligand_db_path).exists():
#
#         zip_file = str(default_ligand_db_path) + ".zip"
#         if not Path(zip_file).exists():
#             raise FileNotFoundError(f"Could not find MetaLig database zip file at {Path(zip_file).resolve()}. Please download it from the DART github repository and place it there.")
#
#         db_dir = Path(default_ligand_db_path).parent
#         uncompress_file(zip_file, db_dir)
#         assert Path(default_ligand_db_path).exists(), f"Could not find MetaLig database at {Path(default_ligand_db_path).resolve()}. Please download it from the DART github repository and place it there."
#         print(f"Uncompressed MetaLig database to {Path(default_ligand_db_path).resolve()}.\n")
#
#         if delete_zip:
#             Path(zip_file).unlink()
#
#     return
#
# def uncompress_file(zip_file_path, output_dir):
#     """
#     Uncompress a zip file into a file.
#     """
#     # Create a temporary directory to extract the files
#     temp_dir = Path(output_dir, 'temp_extract')
#
#     # Extract all files to the temporary directory
#     with zipfile.ZipFile(zip_file_path, 'r') as zipf:
#         zipf.extractall(temp_dir)
#
#     # Move the contents of the temporary directory to the output directory
#     for filename in os.listdir(temp_dir):
#         shutil.move(Path(temp_dir, filename), Path(output_dir, filename))
#
#     # Remove the temporary directory
#     Path(temp_dir).rmdir()
#
#     return


from setuptools import setup
from setuptools.command.develop import develop
from setuptools.command.install import install


def friendly(command_subclass):
    """A decorator for classes subclassing one of the setuptools commands.

    It modifies the run() method so that it prints a friendly greeting.
    """
    orig_run = command_subclass.run

    def modified_run(self):
        print("Hello, developer, how are you? :)")
        orig_run(self)

    command_subclass.run = modified_run
    return command_subclass

...

@friendly
class CustomDevelopCommand(develop):
    pass

@friendly
class CustomInstallCommand(install):
    pass