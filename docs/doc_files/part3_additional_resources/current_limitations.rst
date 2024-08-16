.. _current_limitations:

Future Work and Limitations
=============================

In the current version of DART, we have identified a few areas for future work. If any of them are important for your research, please `let us know <https://github.com/CCEMGroupTCD/DART/issues>`_ and we will try to prioritize them asap.

Assembler:
    - Ligands which cannot yet be assembled
        - Hydride (H)
        - Tridentate/tetradentate ligands with non-planar coordinating atoms
        - Haptic ligands (ligands with neighboring coordinating atoms)
    - We have implemented the most important octahedral and square planar geometries. Let us know if you need any other geometries for your research and we will add them (e.g. octahedral (2-2-2, fac-3-2-1, 1-1-1-1-1-1), square planar 1-1-1-1, tetrahedral, ...).

Ligand Filters:
    - Currently, DART supports a lot of ligand filters via pre-defined filters with options. If you need more filters, please let us know.
    - We are also thinking about adding the possibility to write a customizable python function for filtering ligands.

MetaLig Ligand Database:
    - The MetaLig database might be updated with more functionalities for exploring it's ligands and their properties.
    - It might be possible to add support for extending the MetaLig database with user-defined ligands.

Other limitations are more fundamental:

MetaLig Ligand Database:
    - Ligands might have missing Hydrogen atoms or wrong charges. This should be extremely rare, but cannot be fully excluded because of the large number of ligands in the MetaLig database.
    - Unique ligands are defined by their graph when connected to a (pseudo) metal center. Therefore stereo-isomers are not distinguished in the ligand extraction process and only one stereo-isomer would be part of the MetaLig database.



