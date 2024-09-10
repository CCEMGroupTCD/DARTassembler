.. _dart_workflow:

Understanding DART
===============================

.. contents:: :local:

DART in a Nutshell
---------------------------------

Here is a brief overview of how DART works:

**1. The MetaLig ligand database :** The MetaLig contains 41,018 ligands extracted from the Cambridge Structural Database (CSD) and includes accurately determined formal charges. It is a great resource for exploring ligands and generating novel complexes.

**2. Ligand Filters :** The Ligand Filter Module enables you to refine the MetaLig to target ligands that express chemical properties important to your personal research.

**3. Assembler :** The Assembler module generates novel complexes by selecting ligands from the MetaLig and assembling them to a metal center. You can filter ligands, specify a different set of ligands for each binding site, and enforce chemical equivalence across multiple binding sites. The assembler configuration file allows you to specify many options such as the geometry and charge of the complex, the metal center, the ligands to use, how to manage geometric isomers and more.

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

**1. Ligand choice :**
    First, the assembler selects one ligand for each binding site in the specified geometry from all ligands in the specified ligand database with the correct denticity. Depending on the configuration, the ligand choice is either random or iteratively. If the user specified a different ligand database for each binding site, the ligand for each site is chosen from only that database. If the user specified to enforce chemical equivalence across multiple binding sites, the ligand for the second binding site will simply be copied from the previous binding site. Then, all chosen ligands undergo a first check if the sum of formal charges agrees with the specified metal oxidation state and total charge of the complex and if DART :ref:`can assemble all of the ligands <current_limitations>`. If the check fails, the assembler will discard this combination of ligands and continue to the next one.

**2. Ligand placement and rotation :**
    The set of ligands is then coordinated to each binding site of the specified metal center atom via placement and rotation of the ligands implemented using the Python package `stk <https://github.com/lukasturcani/stk>`_. DART has a different algorithm of placement and rotation implemented for each denticity and depending on the ligand shape, such as if a bidentate ligand in it's original CSD complex had a planar or non-planar metallacycle.

**3. Geometrical isomers :**
    For each complex, the assembler automatically generates all geometric isomers. If the user specified to only return the lowest energy isomer, the assembler will calculate the energy of each isomer using a single point calculation with a universal force field and return only the one with the lowest-energy, otherwise all isomers will be processed.

**4. Isomer naming :**
    After assembling all geometric isomers, the assembler automatically generate an 8-letter pronounceable name for each isomer. Per default, this name will be the same for all isomers but appended with an increasing number, or if only the lowest energy isomer is returned, the number will be left away. Alternatively, users can specify in the configuration file to generate a completely new name for each isomer.

**5. Post-assembly check :**
    All isomers undergo a post-assembly check which looks our for steric clashes between ligands by checking if two ligands are more than a certain distance away from each other. If the check fails, the assembler will discard the isomer. If all isomers fail the check, the assembler will continue to the next combination of ligands.

**6. Saving the complex :**
    Finally, the assembler will save the geometry of each isomer to an .xyz file, some info about its ligands to a .csv file and extensive machine-readable information about the isomer to a .json file. The assembler will then proceed to assemble the next complex, until the specified maximum number of complexes is reached or all combinations of ligands are exhausted.

After assembling all required complexes, the assembler will save some overview files, such as a .csv file with information about all attempted complexes and concatenated .xyz files for all successful and failed complexes. The assembler will also print a summary of the assembly process to the console.


