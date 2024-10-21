from setuptools import setup, find_packages

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name="DARTassembler",
    version="1.0.0",
    description="Simple to use package for building 3D structures of novel transition metal complexes from a large database of ligands extracted from the Cambridge Structural Database.",
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/CCEMGroupTCD/DART',
    python_requires=">=3.9",
    author="Timo Sommer, Cian Clarke, Felix Kleuker",
    packages=find_packages(),
    install_requires=requirements,
    package_data={
        'DARTassembler': [
            'data/*.csv',
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
