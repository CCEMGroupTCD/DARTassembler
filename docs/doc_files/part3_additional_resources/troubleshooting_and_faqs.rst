.. _troubleshooting:

Troubleshooting And FAQs
============================

Please check out the FAQs below for common issues and solutions. If you encounter a problem that is not listed here, please `let us know.<https://github.com/CCEMGroupTCD/DART/issues>`_

Potential issues and solutions
------------------------------

DART cannot be found :
    First, check if conda is installed by running the command ``conda --help``. If this is not successful, refer to the next section.

    If the command ``dart --help`` doesn't display a list of possible DART commands, then either DART is not installed or the correct conda environment is not activated. Run ``conda activate base`` (or replace `base` with the name of the conda environment in which you installed DART) and try again. To see a list of all conda environments, run ``conda env list``. Please refer to the :ref:`installation_guide` for more details.

    It can also happen that you moved the DART folder or that one of the parent folders changed its name. In this case, you need to re-install DART since the path to the DART folder is hard-coded in the installation process.

conda cannot be found :
    If the command ``conda --help`` doesn't print a list of available conda commands, then either conda is not installed or it wasn't added to the PATH. A common issue is that conda is downloaded, but it was not added to the PATH, which means it won't be activated by default in a new terminal session. Please retry the installation of conda or google how to install conda correctly. If you re-install conda, remember to choose 'yes' when it asks if you want to add conda to the PATH.

The input yaml file cannot be read:
    Yaml files are quite sensitive to indentation. Make sure that the indentation is correct. Ideally, copy paste the example yaml file and modify it.

The input yaml file cannot be found:
    Make sure that the path to the yaml file is correct.

The input yaml file is not valid:
    Verify that the options in the YAML file are set up correctly to match the criteria desired. Ensure that the syntax of the yaml file, including lists and nested elements, is correct.

FAQs
----

Updating DART:
    If DART is not performing as documented, ensure that you have the latest version. Updates may include bug fixes and performance improvements.
