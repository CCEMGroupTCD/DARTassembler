.. _installation_guide:

Installation
======================

DART requires about 470MB of space due to the integrated MetaLig database. At the moment, DART supports Linux and MacOs on Python 3.9 to Python 3.12, and support for Windows will be released in the following days. To easily install DART with ``pip``, please run the following command in your terminal:

.. code-block:: sh

    pip install DARTassembler

If you're not sure what is meant with the terminal or ``pip``, or if ``pip`` is not yet installed on your computer, please google or refer to the :ref:`troubleshooting`.

.. tip::

    It is recommended to install DART as regular user, not as sudo user, i.e. please use ``pip install`` instead of ``sudo pip install``.

To verify that DART is installed correctly, run the following command:

.. code-block:: sh

   DARTassembler installtest --path

This command will run a brief test of the ``assembler``, ``ligandfilters`` and ``dbinfo`` module of DART, without saving any files. If this command displays any errors, please `let us know immediately so we can fix it. <https://github.com/CCEMGroupTCD/DART/issues>`_ Also, if you run DART for the very first time, it may be that there is a short delay of around 30s while ``pip`` sets up everything. That is normal and only happens once.

Ready to go? Let's start off with the :ref:`quickstart` or check out the :ref:`DART modules <module_overview>`!