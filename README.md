
# DART (Directed Assembly of Random Transition metal complexes)
This repository contains the Python code for DART, a chemistry software for easily assembling novel metal complexes from a large database of ligands extracted from the Cambridge Structural Database. In its current implementation, DART supports the assembly of a wide range of octahedral and square-planar geometries from 41,018 different ligands with assigned charges. A major focus of DART lies on targeting specific chemical spaces by offering multiple filters to condense the number of ligands used for the assembly.

## Using DART
todo

### Installing DART
We are currently working on providing DART as python package on conda-forge to make the installation a matter of seconds. In the meantime, it takes a few minutes to set up DART. We recommend to install DART by downloading the repository from github and installing it as local python package with the correct dependencies as specified in the conda environment file. In the following we will go through these steps one-by-one. Please open the terminal on your computer and copy and execute the following commands:

1. Install conda on your computer. Conda is a safe everyday tool for Python users and can easily be installed as explained in `https://conda.io/projects/conda/en/latest/user-guide/install/index.html`. Download Miniconda by following the corresponding link to the installation guide depending on if you use Windows, Mac or Linux.
2. Download the DART repository into the current directory
   ```sh
   git clone todo
   ```
3. Set up the conda environment for DART downloading all necessary dependencies
   ```sh
   conda env create -f DART\conda\DART_env.yml      <-- for Windows users
   conda env create -f DART/conda/DART_env.yml      <-- for Mac/Linux users
   ```
   Note: The conda environment will be linked to the absolute path to the cloned version of this repo. If you move the downloaded DART folder, you need to redo this step to link the conda environment again to the correct path.
4. Activate the conda environment:
   ```sh
   conda activate DART
   ```
Now DART is installed as python package in the conda environment DART. Whenever you want to use DART, you need to be in this conda environment. In doubt, just type 'conda activate DART' into your terminal, this brings you into the correct conda environment. Then, you can use all DART commands.


## How to cite DART
Please cite our paper as given in [1].


## Reproducing the DART paper
todo


## License
... TODO
The 3DSC<sub>MP</sub> database is subject to the Creative Commons Attribution 4.0 License, implying that the content may be copied, distributed, transmitted, and adapted, without obtaining specific permission from the repository owner, provided proper attribution is given to the repository owner. All software in this repository is subject to the MIT license. See `LICENSE.md` for more information.


## Origin of data
We are grateful to the providers of the Cambridge Structural Database[2], which is the source of all ligands in the ligand database.


## References
todo
[1] our own paper
[2] Cambridge Structural Database

