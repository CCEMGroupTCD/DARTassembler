.. _dart_workflow:

Understanding the DART Workflow
===============================

DART is an advanced computational platform designed for the rapid generation and exploration of novel transition metal complexes from a database of 41,018 unique ligands. In this section, we will walk you through the ideas and concepts behind DART and it's implementation.

The DART Principles - Users First
----------------------------------

DART is designed to be accessible and user-focused without requiring Python programming skills. Instead, DART works by calling modules from the command line, such as ``DARTassembler dbinfo --path metalig``, which will immediately :ref:`output a .csv file and an .xyz file <metalig_ligand_properties>`  about the MetaLig database into your current directory.

To keep things simple, DART is all about files. You can move them around, copy them, and delete them as you see fit. DART has a few different types of files:

    - **.jsonlines :** These files are ligand databases. You need to provide one or more of these to the :ref:`DART assembler module <assembly_input>`.
    - **.yml :** These files are configuration files for Ligand Filters or Assembler. They allow you to specify a wide range of ligand filters and assembly options to customize your output to your needs.
    - **.xyz :** These files are outputted by the dbinfo and the assembler module. They contain structures of ligands or complexes respectively. If they start with "concat", they contain multiple structures. All .xyz files can easily be visualized for example with ase with the command ``ase gui FILE.xyz``.
    - **.csv :** These files are outputted by the dbinfo and the assembler module. They contain information about ligands and complexes respectively.
    - **.json :** These files are outputted by the assembler. They contain machine-readable information about each complex.

The DART Workflow
-----------------

Here is a brief overview of the DART workflow:

**1. The MetaLig ligand database :** The MetaLig contains 41,018 ligands extracted from the Cambridge Structural Database (CSD) and includes accurately determined formal charges. It is the starting point for all DART workflows.

**2. Ligand Filters :** The Ligand Filter module helps you refine the MetaLig to get a smaller ligand database with exactly defined ligand properties. Nobody needs the entire MetaLig database for their research project, so the Ligand Filter module helps you to select the ligands relevant for you.

**3. Assembler :** The Assembler module generates complexes from a specified ligand dataset. You can even specify a different ligand dataset for each binding site or enforce chemically equivalent ligands across multiple binding sites with the same denticity. The assembler configuration file allows you to specify loads of options such as the geometry and charge of the complex, the metal center, the ligands to use, how to handle isomers and more.
