
# DART - Directed Assembly of Random Transition metal complexes
Welcome to the DART platform, a cutting-edge suite of tools for the generation and exploration of mono-metallic transition metal complexes! Developed by the CCEM group at Trinity College Dublin, DART is engineered to facilitate the design and analysis of molecular complexes for chemistry research.

DART integrates a collection of modules, each serving a unique function in the assembly process:

- **MetaLig** :
    Explore the comprehensive MetaLig database with 41,018 ligands extracted from the Cambridge Structural Database, complete with high-quality formal charge assignments.

- **Assembler** :
    Assemble novel transition metal complexes in a matter of seconds, guided by a simple configuration file for precise control over the resulting structures.

- **Ligand Filters** :
    Tailor assembled ligands to your research needs with a wide range of chemical and data-driven filters.

Using DART is simple. After download, just run the DART assembler and start generating complexes by executing the following command in your terminal:

    DARTassembler assembler --path assembly_input.yml

## Documentation
The documentation of DART (TODO: LINK) will show you how to install and use DART.

## How to cite DART
Please cite our paper as given in [1].

## Reproducing the DART paper
The documentation contains a section about the Pd/Ni cross coupling example that was used in the paper. It will show you how to reproduce the results of the paper.

## License
DART is subject to the MIT license. See `LICENSE.md` for more information.

## Origin of data
We are grateful to the providers of the [Cambridge Structural Database](https://www.ccdc.cam.ac.uk/structures/)[2], which is the source of all ligands in the ligand database.

## References
1. our own paper # TODO
2. Groom, C. R., Bruno, I. J., Lightfoot, M. P. & Ward, S. C. The Cambridge Structural Database. Acta Crystallogr B Struct Sci Cryst Eng Mater 72, 171–179 (2016) 

