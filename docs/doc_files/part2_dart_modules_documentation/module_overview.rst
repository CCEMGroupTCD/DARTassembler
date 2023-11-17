.. _module_overview:

Module Overview
================

DART implements a list of commands for the terminal which you can use to interact with the package:

``dart --help`` :
    Displays a list of available commands. Useful to check if DART is available.

``dart assembler --path assembler_input.yml`` :
    Initiates the :ref:`assembler module <assembly_input>` using the specified configuration file which is named `assembler_input.yml` here. This will run the assembler and save the output to the directory specified in the configuration file.

``dart ligandfilters --path ligandfilters_input.yml`` :
    Executes the :ref:`ligand filters module <ligandfilters>` based on the provided configuration file which is named `ligandfilters_input.yml` here. This will run the filter and save the output to the directory specified in the configuration file.

``dart dbinfo --path ligand_db.jsonlines``:
    Saves an overview table of the specified ligand database file which is named `ligand_db.jsonlines` here. The table will be saved as csv file with the same filename, e.g. here `ligand_db.csv`. If the input file is the word `MetaLig`, the table will be generated for the entire :ref:`MetaLig <metalig>` database.