from setuptools import setup, find_packages

setup(
    name="DARTassembler",
    version="1.0.0",
    description="Simple to use package for building 3D structures of novel transition metal complexes from a large database of ligands extracted from the Cambridge Structural Database.",
    url='https://github.com/CCEMGroupTCD/DART',
    python_requires=">=3.9",
    author="Timo Sommer, Cian Clarke, Felix Kleuker",
    packages=find_packages(),
    install_requires=['ase', 'jsonlines', 'networkx', 'numpy', 'openbabel-wheel', 'pandas', 'pyyaml', 'rdkit', 'scipy', 'scikit-learn', 'stk', 'sympy', 'tqdm', 'pysmiles'],
    package_data={
        'DARTassembler': [
            'data/*.csv',
            'data/metalig/*.zip',
            'data/metalig/*.bz2',
            'data/tests/test_installation/*',
            'data/default/*'
        ]
    },
    entry_points={
        'console_scripts': [
            'DARTassembler=DARTassembler.dart_cli:main',
        ],
    },
)
