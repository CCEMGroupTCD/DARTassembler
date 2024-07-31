from setuptools import setup, find_packages

setup(
    name="DARTassembler",
    version="0.1",
    description="A software for chemists to build novel metal complexes from random ligands on the computer.",
    url='https://github.com/CCEMGroupTCD/DART',
    python_requires=">=3.9",
    author="Timo Sommer, Cian Clarke, Felix Kleuker",
    packages=find_packages(),
    install_requires=[],    # requirements are handled by conda
    entry_points={
        'console_scripts': [
            'dart=dart_cli:main',
        ],
    },
)
