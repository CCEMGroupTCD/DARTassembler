
# DART - Directed Assembly of Random Transition metal complexes
Welcome to the DART platform, a cutting-edge suite of tools for the exploration of coordination chemistry! Developed by the CCEM group at Trinity College Dublin in Ireland & CIC energiGUNE in Spain, DART is designed as an accessible and simple-to-use software to generate transition metal complexes based on ligands from decades of crystallographic data.

DART integrates a collection of several modules:

- **MetaLig Ligand Database :**
    Explore the comprehensive MetaLig database with 41,018 ligands extracted from the Cambridge Structural Database, complete with high-quality formal charge and ligand coordination archetype assignments.

- **Assembler :**
    Assemble novel transition metal complexes in seconds from 22 different ligand coordination archetypes, supporting even haptic and multi-metallic systems.

- **LigandFilters :**
    Assemble complexes with exactly defined sub-structures by applying advanced ligand filters for each binding site.

Using DART is simple. After download, just run the DART assembler and start generating complexes by executing the following command in your terminal:

    DARTassembler assembler --path assembler.yml

## Documentation and Examples
The documentation of DART under [https://dartassembler.readthedocs.io](https://dartassembler.readthedocs.io) will show you how to install and use DART. It contains a quickstart guide and an advanced example, walking you through the main features of DART:

 - Assemble arbitrary transition metal complexes from 22 different ligand coordination archetypes.
 - Use advanced ligand filters to select ligands with specific sub-structures at each binding site.
 - Generate multi-metallic complexes.
 - Assemble complexes with haptic ligands.
 - Browse and search 41,018 ligands in the MetaLig ligand database.

## How to cite DART
Please cite [our paper](https://chemrxiv.org/engage/chemrxiv/article-details/6717eb4e83f22e4214d2b98b).

## License
DART is subject to the GPL-3.0-or-later license. See `LICENSE` for more information.

## Origin of data
We are grateful to the providers of the [Cambridge Structural Database](https://www.ccdc.cam.ac.uk/structures/), which is the source of all ligands in the ligand database.

