.. _installation_guide:

Installation
======================

DART requires about 38MB of space due to the integrated MetaLig database. DART supports Linux, macOS and Windows on Python >=3.9. To easily install DART with ``pip``, please run the following command in your terminal:

.. code-block:: bash

    pip install DARTassembler

If you're not sure what is meant with the terminal or ``pip``, or if ``pip`` is not yet installed on your computer, please google or refer to the :ref:`troubleshooting`.

.. tip::

    It is recommended to install DART as regular user, not as sudo user, i.e. please use ``pip install`` instead of ``sudo pip install``.

Verify the installation
------------------------------------------

After the installation is complete, please verify that DART was installed correctly by running the following command in your terminal:

.. code-block:: bash

   DARTassembler --help

This command will display the help message of DART, confirming that the installation was successful. Please check that the installed DART version is at least ``1.1.0``: this version introduced major new features such as multi-metallic complexes, haptic ligands, and support for arbitrary complex geometries from 22 different ligand coordination archetypes. We have also made DART significantly faster and removed several dependencies, so ``pip`` should usually have no problem installing the latest DART version.

If this command displays any errors, please check the :ref:`FAQs <troubleshooting>` and if you don't find the answer, please `let us know so we can fix it. <https://github.com/CCEMGroupTCD/DART/issues>`_ Also, if you run DART for the very first time, it may be that there is a short delay of around 30s while ``pip`` sets up everything. That is normal and only happens once.

Ready to go? Let's start off with the :ref:`quickstart` or check out the :ref:`DART modules <module_overview>`!