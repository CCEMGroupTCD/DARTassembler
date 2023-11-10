.. _troubleshooting:

Troubleshooting And FAQs
============================

Potential issues and solutions
------------------------------

DART cannot be found :
    Either DART is not installed or the correct conda environment is not activated. Run ``conda activate YOUR_ENV`` (where YOUR_ENV is the name of the conda environment in which you installed DART) and try again. To see a list of all conda environments, run ``conda env list``. Please refer to the :ref:`installation_guide` for more details.

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
