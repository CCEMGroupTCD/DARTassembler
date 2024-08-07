from setuptools import setup, find_packages
from setuptools.command.install import install
from setuptools.command.develop import develop
from setuptools.command.egg_info import egg_info
from pathlib import Path
import os
import zipfile
import sysconfig
import glob

def extract_metalig_zipfile() -> None:
    """
    Unzips the MetaLig database zip file after installation. Removes the zip file after unzipping.
    :return: None
    """
    print("Executing post-installation tasks...")
    # Get the path to the metalig zip files
    install_path = sysconfig.get_paths()['purelib']
    data_dir = Path(install_path, 'DARTassembler', 'data', 'metalig')
    metalig_files_to_unzip = glob.glob('MetaLigDB_v*.jsonlines.zip', root_dir=data_dir)

    for filename in metalig_files_to_unzip:
        file_path = Path(data_dir, filename)
        if file_path.exists():
            print(f"Unzipping {file_path}...")
            with zipfile.ZipFile(file_path, 'r') as zip_ref:
                zip_ref.extractall(os.path.dirname(file_path))
            file_path.unlink() # Remove the zip file after unzipping
            print(f"{file_path} unzipped successfully.")
        else:
            print(f"{file_path} not found.")


class CustomInstallCommand(install):
    def run(self):
        install.run(self)
        extract_metalig_zipfile()


class CustomDevelopCommand(develop):
    def run(self):
        develop.run(self)
        extract_metalig_zipfile()


class CustomEggInfoCommand(egg_info):
    def run(self):
        egg_info.run(self)
        extract_metalig_zipfile()


setup(
    name="DARTassembler",
    version="0.1",
    description="Simple to use package for building 3D structures of novel transition metal complexes from a large database of ligands extracted from the Cambridge Structural Database.",
    url='https://github.com/CCEMGroupTCD/DART',
    python_requires=">=3.9",
    author="Timo Sommer, Cian Clarke, Felix Kleuker",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[],  # Requirements are handled by conda. See file conda_DART.yml.
    package_data={
        'DARTassembler': [
            'data/*.csv',
            'data/metalig/*.jsonlines.zip',
            'data/metalig/test*.jsonlines',
            'data/metalig/*.csv',
            'data/tests/test_installation/*'
        ]
    },
    cmdclass={
        'install': CustomInstallCommand,
        'develop': CustomDevelopCommand,
        'egg_info': CustomEggInfoCommand,
    },
    entry_points={
        'console_scripts': [
            'DARTassembler=DARTassembler.dart_cli:main',
        ],
    },
)
