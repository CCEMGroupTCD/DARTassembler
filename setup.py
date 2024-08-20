from setuptools import setup, find_packages

setup(
    name="DARTassembler",
    version="0.1",
    description="Simple to use package for building 3D structures of novel transition metal complexes from a large database of ligands extracted from the Cambridge Structural Database.",
    url='https://github.com/CCEMGroupTCD/DART',
    python_requires=">=3.9",
    author="Timo Sommer, Cian Clarke, Felix Kleuker",
    packages=find_packages(),
    include_package_data=True,
    install_requires=['ase', 'jsonlines', 'networkx', 'numpy', 'openbabel-wheel', 'pandas', 'pbr', 'pyyaml', 'rdkit', 'scipy', 'scikit-learn', 'stk', 'sympy', 'tqdm', 'pysmiles'],  # Requirements are handled by conda. See file conda_DART.yml.
    package_data={
        'DARTassembler': [
            'data/*.csv',
            'data/metalig/*.jsonlines.zip',
            'data/metalig/test*.jsonlines',
            'data/metalig/*.csv',
            'data/tests/test_installation/*'
        ]
    },
    entry_points={
        'console_scripts': [
            'DARTassembler=DARTassembler.dart_cli:main',
        ],
    },
)
