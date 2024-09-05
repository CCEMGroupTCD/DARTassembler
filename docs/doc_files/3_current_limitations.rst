.. _current_limitations:

Future Work and Limitations
=============================

In the current version of DART, we have identified a number of DART features that are potentially interesting but not yet implemented. If any of them are important for your research, please `let us know <https://github.com/CCEMGroupTCD/DART/issues>`_ and we will try to address them asap.

**Assembler :**
  - Ligands which cannot yet be used in the assembler:

    - Hydride (H)
    - Tridentate/tetradentate ligands with non-planar coordinating atoms
    - Haptic ligands (ligands with neighboring coordinating atoms)

  - We have implemented the most important octahedral and square planar geometries. Let us know if you need any other geometries for your research and we will add them (e.g. tetrahedral).

**Ligand Filters :**
  - Currently, DART supports a lot of ligand filters via pre-defined filters with options. If you need more filters, please let us know.

**MetaLig Ligand Database :**
  - The MetaLig database might be updated with more functionalities for exploring it's ligands and their properties.
    - It might be possible to add support for extending the MetaLig database with user-defined ligands.

Other limitations are more fundamental and not likely to be addressed in the near future. These include:

**MetaLig Ligand Database :**
  - Unique ligands are defined by their graph when artificially connected to a Hg metal center. Therefore stereo-isomers are not distinguished in the ligand extraction process and only one stereo-isomer would be part of the MetaLig database.



