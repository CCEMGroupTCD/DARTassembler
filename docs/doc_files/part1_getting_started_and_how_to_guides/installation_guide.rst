.. _installation_guide:

Installation Guide
======================

DART requires about 400MB of space, mostly due to the MetaLig database. The installation process should take about 5-10 minutes.

Easy Installation (Coming Soon)
-------------------------------

We are in the process of making DART available as a Python package on conda-forge for quick and easy installation. Stay tuned for updates!

Manual Installation
-------------------

For now, you can install DART manually by downloading the code from GitHub and setting it up as a local Python package. Below are simple step-by-step instructions to guide you through this process.

Step 1: Open Terminal and Navigate to Installation Directory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First, open your terminal application:

- On Windows, search for Command Prompt or PowerShell in the Start menu.
- On macOS, open the Applications folder, then Utilities, and double-click on Terminal.
- On Linux, you can usually find the Terminal in your applications menu or you can use a keyboard shortcut like Ctrl+Alt+T.

Next, navigate to the directory where you want to install DART. Your home directory is a good choice if you don't have a specific location in mind. To navigate to the home directory use the command ``cd $HOME`` on Linux and macOS, or ``cd %USERPROFILE%`` on Windows.

Step 2: Install Conda
^^^^^^^^^^^^^^^^^^^^^^

Conda is an environment and package manager that simplifies Python installations. If you don't have it already, you can install Conda by following the instructions on the `official Conda installation guide <https://conda.io/projects/conda/en/latest/user-guide/install/index.html>`_. Choose the Miniconda version appropriate for your operating system (Windows, Mac, or Linux). During the installation, you will be asked if you want to add Conda to your PATH variable. Make sure to select 'yes' to this option.

**Ensure that Conda is installed correctly:**
    Run ``conda --version`` to confirm that Conda is installed and working correctly. If this outputs the version number, conda is installed.

Step 3: Download DART from GitHub
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Try to download DART into the current folder using git:

.. code-block:: sh

    git clone https://github.com/CCEMGroupTCD/DART.git

If this gives an error message that git is not installed:
    You can download DART directly from the GitHub website. Visit the `DART repository <https://github.com/CCEMGroupTCD/DART>`_ in your web browser and download the code as a ZIP file by clicking on the big green button that says ``<> Code`` and then ``Download ZIP``. Unzip the file into your current directory.

Step 4: Navigate to the DART Directory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now, change to the DART directory with this command:

.. code-block:: sh

   cd DART

If the command fails, it could be due to an unsuccessful download or DART being in a different directory. Please verify the download location and repeat step 3 if necessary.

Step 5: Install DART into a Conda Environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- **For inexperienced Conda users:**
    Install DART in the base environment of Conda for simplicity. Run the following command:

    .. code-block:: sh

       conda env update --name base --file DART_env.yml

- **For experienced Conda users:**
    If you wish to create a new conda environment for DART, run this command and replace `YOUR_ENV` with your desired environment name:

    .. code-block:: sh

       conda env create --name YOUR_ENV --file DART_env.yml

    Activate the environment with ``conda activate YOUR_ENV``.

.. caution:: Be aware that the Conda environment will link to the absolute path where you've downloaded DART. Avoid moving the DART folder after installation, as this would require reinstalling the environment.

Step 6: Ensure that DART is Installed Correctly
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Run the following command to verify the installation:

.. code-block:: sh

   dart --help

If this command is successful, it will display a list of available DART commands. If not, please refer to the :ref:`troubleshooting section <troubleshooting>` for further assistance.

Usage
-----

After installation, DART is ready for use. Here are some basic terminal commands for interaction:

- ``dart --help``: Displays available commands.
- ``dart assembler --path assembly_input.yml``: Starts the :ref:`assembler module <assembly_input>` with the specified configuration file.
- ``dart ligandfilters --path ligandfilters_input.yml``: Runs the :ref:`ligand filters module <ligandfilters>` with the specified configuration file.
- ``dart dbinfo --path ligand_db.jsonlines``: Generates an overview table of the specified ligand database file and saves it as a CSV file. Use the command ``dart dbinfo --path metalig`` to generate a table for the entire MetaLig database.