.. _module_overview:

Module Overview
================

DART implements a list of commands for the terminal which you can use to interact with the package:

.. option:: DARTassembler --help

    Displays a list of available commands. Useful to check if DART is available.

.. option::  DARTassembler assembler --path assembler_input.yml

    Executes the :ref:`Assembler Module <assembler>` using the provided configuration file.

.. option:: DARTassembler ligandfilters --path ligandfilters_input.yml

    Executes the :ref:`Ligand Filters Module <ligandfilters>` using the provided configuration file.

.. option:: DARTassembler dbinfo --path ligand_db.jsonlines

    Run this command to get more information on the ligands in the ligand database file ``ligand_db.jsonlines``. It will immediately save two files in the current directory, which here would be called ``ligand_db.csv`` and ``concat_ligand_db.xyz``. The .csv file contains information about each ligand in the database such as stoichiometry, donors or formal charge. The concatenated .xyz file contains the 3D structures of all ligands in the database, coordinated to a Cu metal center for visualization purposes. To browse through the structures, you can use the command ``ase gui concat_ligand_db.xyz``.

    Instead of providing a path to a ligand database file, you can also use ``--path metalig`` to generate the overview for the entire :ref:`MetaLig <metalig>` database or ``--path test_metalig`` to generate the overview for a test database including only the first 1000 ligands of the MetaLig database.

.. option:: DARTassembler concat --path ligand_db1.jsonlines ligand_db2.jsonlines

    Concatenates two or more ligand database files into one. The output will be saved to a file called ``concat_ligand_db.jsonlines``.

.. option:: DARTassembler installtest --path testdir

    Runs a quick test of the ``assembler``, ``ligandfilters`` and ``dbinfo`` module and saves the results in the directory ``testdir``. If no directory is specified, the results will not be saved.

.. option:: DARTassembler configs --path .

    Saves two files to the specified directory, ``assembler.yml`` and ``ligandfilters.yml``. These files are templates for the configuration files of the assembler and ligandfilters module. They are handy to get started with a new project.
