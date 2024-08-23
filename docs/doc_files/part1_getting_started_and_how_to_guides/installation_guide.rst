.. _installation_guide:

Installation
======================

To install DART please run the following command in your terminal:

.. code-block:: sh

    pip install "git+https://github.com/CCEMGroupTCD/DART@master#egg=DARTassembler"

DART requires about 400MB of space due to the MetaLig database. If you're not sure what is meant with the terminal or pip, please refer to the :ref:`troubleshooting`.

To verify that DART is installed correctly, run the following command:

.. code-block:: sh

   DARTassembler test --path

This command will run a brief test of the ``assembler``, ``ligandfilters`` and ``dbinfo`` module of DART, without saving any files. If you run DART for the very first time, you will also notice that DART takes around 20 seconds to extract the MetaLig database from the .zip file that was downloaded during installation. That is normal and only happens once.

Ready to go? Let's start off with the :ref:`quickstart` or check out the :ref:`DART modules <module_overview>`!
