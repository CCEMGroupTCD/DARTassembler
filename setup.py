from setuptools import setup, find_packages

version = "1.0.4.1"

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

with open('README.md') as f:
    long_description = f.read()

setup(
    name="DARTassembler",
    version=version,
    description="Simple to use package for building 3D structures of novel transition metal complexes from a large database of ligands extracted from the Cambridge Structural Database.",
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/CCEMGroupTCD/DARTassembler',
    python_requires=">=3.9",
    author="Timo Sommer, Cian Clarke, Felix Kleuker, Max García-Melchor",
    packages=find_packages(),
    install_requires=requirements,
    package_data={
        'DARTassembler': [
            'data/*.csv',
            'data/metalig/MetaLigDB*.bz2',
            'data/default/*'
        ]
    },
    entry_points={
        'console_scripts': [
            'DARTassembler=DARTassembler.src.misc.cli:main',
        ],
    },
    classifiers=[
        # See https://pypi.org/classifiers/ for the full list
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Programming Language :: Python :: 3.14",
        "Operating System :: OS Independent",
    ],
)
