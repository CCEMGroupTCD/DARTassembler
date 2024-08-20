.. _installation_guide:

Installation Guide
======================

The DART installation process is very simple and should take only about 2 minutes. DART requires about 400MB of space due to the MetaLig database.

First, open your terminal and make sure conda is installed. If these instructions are unclear, please refer to the FAQ at the bottom of this page.

Installation Steps
-------------------

Step 1: Download the DART conda file from GitHub
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Run the following command in the terminal to download the DART conda file from GitHub:

.. code-block:: sh

    curl -O https://raw.githubusercontent.com/CCEMGroupTCD/DART/master/conda_DART.yml

If you don't have the ``curl`` command installed, just download the file manually from the link in the command above.

Step 2: Install DART into a Conda Environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Create a new Conda environment with the DART conda file by running the following command:

    .. code-block:: sh

       conda env create --file conda_DART.yml

Activate the environment with ``conda activate DART``. Congratulations! You have successfully installed DART.

.. note::
    The DART installation is now linked to this conda environment. You must activate this environment every time you want to use DART. If you don't want to do that, you can instead install DART into the conda base environment.

Ensure that DART is Installed Correctly
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Run the following command to verify the installation:

.. code-block:: sh

   DARTassembler test --path

This command will run a brief test of the ``ligandfilters`` and ``assembler`` module of DART, without saving any files.

.. note::
    If this command didn't work, you probably forgot to activate the DART environment with ``conda activate DART``.

Usage
-----

After installation, DART is ready for use. Here are some basic terminal commands for interaction:

- ``DARTassembler --help``: Displays available commands.
- ``DARTassembler assembler --path assembler.yml``: Starts the :ref:`assembler module <assembly_input>` with the specified configuration file.
- ``DARTassembler ligandfilters --path ligandfilters_input.yml``: Runs the :ref:`ligand filters module <ligandfilters>` with the specified configuration file.
- ``DARTassembler dbinfo --path ligand_db.jsonlines``: Generates an overview table of the specified ligand database file and saves it as a CSV file. Use the command ``DARTassembler dbinfo --path metalig`` to generate a table for the entire MetaLig database.
- ``DARTassembler test --path outdir``: Runs a brief test of the ``ligandfilters`` and ``assembler`` modules of DART. If no path is specified, no files will be saved. Otherwise, the output will be saved in the specified directory (here in the directory ``outdir``).

Installation FAQ
------------------

What's the terminal?
^^^^^^^^^^^^^^^^^^^^

The terminal is a text-based interface that allows you to interact with your computer using commands. It is also known as the command line or shell. The terminal is a powerful tool that allows you to perform a wide range of tasks, from navigating your file system to installing software.

How to open your terminal?

- On Windows, search for Command Prompt or PowerShell in the Start menu.
- On macOS, open the Applications folder, then Utilities, and double-click on Terminal.
- On Linux, you can usually find the Terminal in your applications menu or you can use a keyboard shortcut like Ctrl+Alt+T.

Install Conda
^^^^^^^^^^^^^^^^^^^^

Conda is an environment and package manager that simplifies Python installations. If you don't have it already, you can install Conda by following the instructions on the `official Conda installation guide <https://conda.io/projects/conda/en/latest/user-guide/install/index.html>`_. Choose the Miniconda version appropriate for your operating system (Windows, Mac, or Linux). During the installation, you will be asked if you want to add Conda to your PATH variable. Make sure to select 'yes' to this option.

**Ensure that Conda is installed correctly:**
    Run ``conda --version`` to confirm that Conda is installed and working correctly. If this outputs the version number, conda is installed.