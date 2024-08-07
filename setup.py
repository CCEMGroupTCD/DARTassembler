from setuptools import setup, find_packages
from setuptools.command.install import install
from setuptools.command.develop import develop
from setuptools.command.egg_info import egg_info
import subprocess
import os
import zipfile

def custom_command():
    print('FANCY SHIIIIIT YESS LEETTTTSSS GOGOGGGGGGOO!!!!!!')
    print("Executing post-installation tasks...")
    data_dir = os.path.join(os.path.dirname(__file__), 'DARTassembler', 'data')
    files_to_unzip = ['metalig/MetaLigDB_v1.0.0.jsonlines.zip']  # Update with actual file names

    for filename in files_to_unzip:
        file_path = os.path.join(data_dir, filename)
        if os.path.isfile(file_path):
            print(f"Unzipping {file_path}...")
            with zipfile.ZipFile(file_path, 'r') as zip_ref:
                zip_ref.extractall(os.path.dirname(file_path))
            print(f"{file_path} unzipped successfully.")
        else:
            print(f"{file_path} not found.")


class CustomInstallCommand(install):
    def run(self):
        install.run(self)
        custom_command()


class CustomDevelopCommand(develop):
    def run(self):
        develop.run(self)
        custom_command()


class CustomEggInfoCommand(egg_info):
    def run(self):
        egg_info.run(self)
        custom_command()


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
