.. _dart_workflow:

Understanding the DART Workflow
===============================

.. contents:: :local:

How DART helps researchers
----------------------------------------------

DART is designed as an accessible and user-focused platform for researchers without requiring Python programming skills. The simplicity of DART is due the integrated MetaLig ligand database:

- If you're interested in a certain type of transition metal complex, you can very easily generate many different instances of these complexes and explore their chemical space.
- It also means that we could make DART very simple to utilize, because users can just submit a file with their preferred settings and DART will do the rest.
- Users can explore and leverage the MetaLig database for their own research projects.

The DART Workflow
-----------------

Here is a brief overview of the DART workflow:

**1. The MetaLig ligand database :** The MetaLig contains 41,018 ligands extracted from the Cambridge Structural Database (CSD) and includes accurately determined formal charges. It is the starting point for all DART workflows.

**2. Ligand Filters :** The Ligand Filter module helps you refine the MetaLig to get a smaller ligand database with exactly defined ligand properties. Nobody needs the entire MetaLig database for their research project, so the Ligand Filter module helps you to select the ligands relevant for you.

**3. Assembler :** The Assembler module generates complexes from a specified ligand dataset. You can even specify a different ligand dataset for each binding site or enforce chemically equivalent ligands across multiple binding sites with the same denticity. The assembler configuration file allows you to specify loads of options such as the geometry and charge of the complex, the metal center, the ligands to use, how to handle isomers and more.

The DART Files
-------------------------------------------

To keep things simple, DART is all about files on your computer which you can move around as you want. DART has a few different types of files:

    - **.jsonlines :** These files are ligand databases. You need to provide one or more of these to the :ref:`DART assembler module <assembler>`.
    - **.yml :** These files are configuration files for Ligand Filters or Assembler. They allow you to specify a wide range of ligand filters and assembly options to customize your output to your needs.
    - **.xyz :** These files are outputted by the dbinfo and the assembler module. They contain structures of ligands or complexes respectively. If they start with "concat", they contain multiple structures. All .xyz files can easily be visualized for example with ase with the command ``ase gui FILE.xyz``.
    - **.csv :** These files are outputted by the dbinfo and the assembler module. They contain information about ligands and complexes respectively.
    - **.json :** These files are outputted by the assembler. They contain machine-readable information about each complex.


