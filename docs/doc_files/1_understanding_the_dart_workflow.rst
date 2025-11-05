.. _dart_workflow:

Understanding DART
===============================

.. contents:: :local:

DART in a Nutshell
---------------------------------

Here is a brief overview of how DART works:

**1. The MetaLig ligand database :** The MetaLig contains 41,018 ligands extracted from the Cambridge Structural Database (CSD) and includes accurately determined formal charges and :ref:`ligand coordination archetypes <ligand_archetypes>`. It is a great resource for exploring ligands and generating novel complexes.

**2. LigandFilters :** The LigandFilters Module enables you to refine the MetaLig to target ligands that express chemical properties important to your personal research.

**3. Assembler :** The Assembler module generates novel complexes by selecting ligands from the MetaLig and assembling them to a metal center. It supports 22 different ligand coordination archetypes and both haptic and multi-metallic systems. You can filter ligands, specify a different set of ligands for each binding site, and enforce chemical equivalence across multiple binding sites. The assembler configuration file allows you to specify many options such as the geometry and charge of the complex, the metal center, the ligands to use, how to manage geometric isomers and more.

The DART Files
-------------------------------------------

DART has a few different types of files:

    - **.jsonlines :** These files are ligand databases. You need to provide one or more of these to the :ref:`DART assembler module <assembler>`.
    - **.yml :** These files are configuration files for the Ligand Filters and Assembler Module. They allow you to specify a wide range of ligand filters and assembly options to customize your output to your needs.
    - **.xyz :** These files are outputted by the dbinfo and the assembler module. They contain structures of ligands or complexes respectively. If they start with "concat", they contain multiple structures. All .xyz files can easily be visualized for example with ase with the command ``ase gui FILE.xyz``.
    - **.csv :** These files are outputted by the dbinfo and the assembler module. They contain information about ligands and complexes respectively.
    - **.json :** These files are outputted by the assembler. They contain machine-readable information about each complex.

.. _how_assembler_works:

How the DART Assembler Module works
-------------------------------------------

When running the DART Assembler Module, it will first read in the provided configuration file and load the specified ligand database files. Then, for each batch specified in the configuration file the assembler starts to generate complexes:

Ligand sampling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First, the assembler samples one ligand for each binding site in the specified geometry from all ligands in the specified ligand database with the correct denticity. Depending on the configuration, the ligand choice is either random or iteratively. If the user specified a different ligand database for each binding site, the ligand for each site is chosen from only that database. If the user specified to enforce chemical equivalence across multiple binding sites using the keyword ``same_as_previous``, the ligand for the second binding site will simply be copied from the previous binding site. Then, all chosen ligands undergo a first check if the sum of formal charges agrees with the specified metal oxidation state and total charge of the complex and if DART :ref:`can assemble all of the ligands <current_limitations>`. If the check fails, the assembler will discard this combination of ligands and continue to the next one.

Ligand placement and rotation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The set of ligands is then coordinated to each binding site of the specified metal center atom via placement and rotation of the ligands. DART does not use quantum-mechanical methods or force-fields, but the empirical ligand-metal orientation from the CSD source structures. Each ligand in the MetaLig is normalized such that its original metal center is at the origin. Therefore, to place a ligand at a binding site, the assembler first translates the ligand such that its metal center is at the binding site coordinates. Then, the assembler rotates the ligand such that its metal-donor vectors align with the user-provided :ref:`target_vectors <target_vectors>` for that binding site. This procedure ensures that the ligand is placed in a chemically reasonable orientation, bypassing the inaccuracies and high computational cost of force-field or quantum-mechanical based geometry optimizations.

Geometrical isomers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For each complex, the assembler automatically generates all geometric isomers. One can also choose to permute ligands with the same ligand archetype.

AssembledIsomer naming
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After assembling all geometric isomers, the assembler automatically generates an 8-letter pronounceable name for each isomer. Per default, this name will be the same for all isomers but appended with an increasing number, or if only the lowest energy isomer is returned, the number will be left away. Alternatively, users can specify in the configuration file to generate a completely new name for each isomer.

Post-assembly check
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All isomers undergo a post-assembly check which looks out for steric clashes between ligands by checking if any two ligands are closer than a threshold distance. If the check fails, the assembler will discard the isomer. It will also check whether two isomers are duplicates of each other (e.g. due to symmetric ligand arrangements) and discard duplicates.

Saving the complex
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Finally, the assembler will save all information of the generated complex, including the 3D structure of all of it's isomers, to a .json file. After assembling all required complexes, the assembler will save some overview files, such as a .csv file with information about all attempted complexes and concatenated .xyz files for all successful and failed complexes. The assembler will also print a summary of the assembly process to the console.


