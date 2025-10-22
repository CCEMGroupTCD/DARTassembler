.. _current_limitations:

Limitations and Recommendations
================================

In the current version of DART, we have identified a number of DART features that are potentially interesting but not yet implemented. While we continue working on further developing DART, please `let us know by creating a new issue on GitHub <https://github.com/CCEMGroupTCD/DART/issues>`_ if any of them are important for your research and we will try to address them ASAP.

**Assembler :**
  - Some types of ligands cannot yet be used in the assembler:

    - Hydride (H)
    - Tridentate/tetradentate ligands with non-planar coordinating atoms
    - Ligands with haptic interactions

  - We have implemented the most important octahedral and square planar geometries. Let us know if you need any other geometries for your research and we will add them (e.g. tetrahedral).
  - Currently, ligands for the assembler are sourced from the extensive MetaLig ligand database, but it is possible to add support for user-provided ligands in the future.

**Ligand Filters :**
  - Currently, DART supports a lot of ligand filters via pre-implemented filter parameters. You can also :ref:`write customized filters via Python <metalig_python_filtering>`, but if you see use cases for additional pre-implemented filters, please let us know.

Other limitations are more fundamental and not likely to be addressed in the near future. These include:

**MetaLig Ligand Database :**
   - Unique ligands are defined by their graph when artificially connected to a Hg metal center. Therefore stereo-isomers are not distinguished in the ligand extraction process and only one stereo-isomer would be part of the MetaLig database.
   - Ligands will never contain any transition metals from the d- or f-block since these would have been considered as bi-metallic complexes in the ligand extraction process.
   - Strongly redox active ligands are likely to not be included in the MetaLig database due to the charge assignment process in which ligands are excluded if the calculated charge is inconsistent between multiple instances of the same ligand.
   - The formal charge assignment process described in the DART paper was benchmarked to be highly accurate, but due to the large size of the MetaLig database it is possible that in very rare occasions ligands would be assigned an incorrect formal charge.

**Assembler :**
   - The assembler generates only mono-metallic complexes.
   - The assembler generates geometric isomers, but not other stereo-isomers or conformers.




