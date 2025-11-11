=============================================================
DART: Directed Assembly of Random Transition Metal Complexes
=============================================================

Welcome to the DART platform, a cutting-edge suite of tools for the exploration of coordination chemistry! Developed by the CCEM group at Trinity College Dublin in Ireland & CIC energiGUNE in Spain, DART is designed as an accessible and simple-to-use software to generate transition metal complexes based on ligands from decades of crystallographic data. If you use DART in your research, please cite the following paper:

    Clarke, C., Sommer, T., Kleuker, F., García-Melchor, M. DART: Unlocking Coordination Chemistry Beyond the Cambridge Structural Database. Preprint at https://doi.org/10.26434/chemrxiv-2024-tljj9 (2024).

DART integrates a collection of several modules:

- :ref:`MetaLig <metalig>` :
    Explore the comprehensive MetaLig database with 41,018 ligands extracted from the Cambridge Structural Database, complete with high-quality formal charge and ligand coordination archetype assignments.

- :ref:`Assembler <assembler>` :
    Assemble novel transition metal complexes in seconds from 22 different ligand coordination archetypes, supporting even haptic and multi-metallic systems.

- :ref:`LigandFilters <ligandfilters>` :
    Assemble complexes with exactly defined sub-structures by applying advanced ligand filters for each binding site.

**Your role:** DART is designed for researchers in chemistry. While you have to think about the chemistry yourself, DART will allow you to explore the chemical space you want with the ligands you need.

If you use DART in your research, please cite the following paper:

    Clarke, C., Sommer, T., Kleuker, F., García-Melchor, M. DART: Unlocking Coordination Chemistry Beyond the Cambridge Structural Database. Preprint at https://doi.org/10.26434/chemrxiv-2024-tljj9 (2024).

.. figure:: /_static/homepage_picture.png
   :width: 90%
   :align: center

   A DART user generates novel molecular complexes for their research on square planar Ni(II) complexes.

Who?
=====================

- **Experimental Chemists**: Explore the effect of varying ligands in your complexes.
- **Computational Chemists**: Model complexes from a user-defined chemical space.
- **Cheminformatics Developers**: Generate novel complexes for high-throughput screening.

Why?
=============

- **Advanced Complex Control**: DART supports a near-infinite variety of transition metal complexes, including geometrical isomers, haptic ligands, and multi-metallic systems.

- **Intuitive Interface**: Users can interact with DART on the terminal via simple configuration files, no Python required.

- **Open-Source**: DART is an open-source tool developed in Python, making it entirely free to use.

Getting Started
===============

Are you ready to get started? :ref:`Install DART <installation_guide>`, explore our hands-on :ref:`quickstart guide <quickstart>` or read through the :ref:`ideas and concepts <dart_workflow>` behind DART.

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   ./doc_files/1_installation_guide
   ./doc_files/1_quickstart
   ./doc_files/1_NaFe_example
   ./doc_files/1_understanding_the_dart_workflow

.. toctree::
   :maxdepth: 2
   :caption: DART Modules

   ./doc_files/2_metalig_database
   ./doc_files/2_ligand_filters
   ./doc_files/2_assembler
   ./doc_files/2_module_overview

.. toctree::
   :maxdepth: 2
   :caption: Advanced DART API

   ./doc_files/2_ligands
   ./doc_files/2_isomers
   ./doc_files/2_metalig_api
   ./doc_files/2_ligand_filters_api
   ./doc_files/2_other_modules

.. toctree::
   :maxdepth: 1
   :caption: Additional Resources

   ./doc_files/3_troubleshooting_and_faqs
   ./doc_files/3_tips_and_tricks
   ./doc_files/3_current_limitations
   ./doc_files/3_version_history

