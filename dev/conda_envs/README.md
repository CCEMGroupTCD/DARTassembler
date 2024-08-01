=== Overview of conda environments ===

Important conda environments:
* conda_DART_dev_exact.yml: To use for all developers. This contains the exact packages to install so that we should all be able to run everything. Will need to be updated regularly.
* conda_DART_dev.yml: For developers, but here the package versions are not specified. Additionally to the user environment, contains packages needed for the documentation and the build of the DART package. Will need to be updated sometimes.
* conda_DART.yml: For users. Do not use for development because this env is as minimalistic as possible, just what you need to run DART.

Old conda environments:
These are supposed to stay as backups and to be able to reproduce early development environments, in case something goes wrong, like e.g. the fact that the charge assignment needs scipy<=1.10.1.
* conda_env_conda_mols_rdkit.yml: Original development environment for this project. Should have mostly correct versions of all packages, including scikit-learn for the charge assignment.
* conda_env_try_to_minimize.yml: Attempt to minimize the environment to only the packages needed for the main DART functionality, i.e. ligand filters and assembler.
