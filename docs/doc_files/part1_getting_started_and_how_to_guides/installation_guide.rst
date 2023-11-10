.. _installation_guide:

Installation Guide
======================

Easy Installation (Coming Soon)
-------------------------------

We're in the process of making DART available as a Python package on conda-forge for quick and easy installation. Stay tuned for updates!

Manual Installation
-------------------

For now, you can install DART manually by downloading the code from GitHub and setting it up as a local Python package. Below are simple step-by-step instructions to guide you through this process.

.. note:: Before you begin, make sure to open the terminal application on your computer. On Windows, you can use either the Command Prompt or PowerShell. On macOS or Linux, use the Terminal.

Step 1: Install Conda
^^^^^^^^^^^^^^^^^^^^^^

Conda is an environment and package manager that simplifies Python installations. If you don't have it already, you can install Conda by following the instructions on the `official Conda installation guide <https://conda.io/projects/conda/en/latest/user-guide/install/index.html>`_. Choose the Miniconda version appropriate for your operating system (Windows, Mac, or Linux).

Step 2: Download DART from GitHub
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- **If you have Git installed:**
  Open your terminal application, navigate to your home directory, and run the following command:

  .. code-block:: sh

     git clone git@github.com:CCEMGroupTCD/DART.git

- **If you don't have Git installed:**
  Visit the repository URL ``git@github.com:CCEMGroupTCD/DART`` in your web browser and download the code as a ZIP file. Unzip the file in your home directory.

  .. note:: For Windows users: If you don't have Git, you can download it from the `official site <https://git-scm.com/download/win>`_.

Step 3: Navigate to the DART Directory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In your terminal, change to the DART directory with this command:

.. code-block:: sh

   cd DART

Step 4: Set Up the Conda Environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. todo: add that if a user is not very experienced in conda, they should install DART and its dependencies in the base environment

Run the following command to set up a new Conda environment and install all the necessary dependencies for DART:

.. code-block:: sh

   conda env create -f DART_env.yml

.. note:: The Conda environment will link to the specific location where you've downloaded DART. If you move the DART folder later, you'll need to repeat this step.

Step 5: Activate the Conda Environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To activate the Conda environment, run:

.. code-block:: sh

   conda activate DART

Usage
-----

After completing the installation, DART will be available within the Conda environment named ``DART``. Here are some basic commands for the terminal you can use to interact with the package:

- ``dart --help``: Displays a list of available commands.
- ``dart assembly --path assembly_input.yml``: Initiates the :ref:`assembler module <assembly_input>` using the specified configuration file.
- ``dart ligandfilters --path ligandfilters_input.yml``: Executes the :ref:`ligand filters module <ligandfilters>` using the specified configuration file.
- ``dart dbinfo --path ligand_db.jsonlines``: Generates an overview table of the specified ligand database file and saves it to a csv file. If the file is the word `MetaLig`, the table will be generated for the entire MetaLig database.

Troubleshooting: 'DART cannot be found' Error
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you encounter issues running DART commands, make sure that you've activated the correct Conda environment. You can activate the DART environment by entering ``conda activate DART`` in your terminal.

To verify that DART is correctly installed and accessible, run the command ``dart --help`` in your terminal. If the installation was successful, this will display a list of available DART commands.

For more information on troubleshooting, please see the :ref:`troubleshooting section <troubleshooting>`.




