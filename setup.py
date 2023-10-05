from setuptools import setup, find_packages

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name="DART-Assembly",
    version="0.1",
    description="A software for chemists to build novel metal complexes from random ligands on the computer.",
    url='https://github.com/CCEMGroupTCD/DART',
    python_requires=">=3.9.5",
    author="Cian Clarke, Timo Sommer, Felix Kleuker",
    packages=find_packages(),
    install_requires=requirements,
)
