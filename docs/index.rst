=============================================================
DART: Directed Assembly of Random Transition Metal Complexes
=============================================================

Welcome to the DART platform, a cutting-edge suite of tools for the exploration of coordination chemistry! Developed by the CCEM group at Trinity College Dublin, DART is designed as an accessible and simple-to-use software to generate mono-metallic transition metal complexes based on ligands from decades of crystallographic data.

DART integrates a collection of several modules:

- :ref:`MetaLig <metalig>` :
    Explore the comprehensive MetaLig database with 41,018 ligands extracted from the Cambridge Structural Database, complete with high-quality formal charge assignments.

- :ref:`Assembler <assembler>` :
    Assemble novel transition metal complexes in a matter of seconds, guided by a single configuration file for precise control over the resulting structures.

- :ref:`Ligand Filters <ligandfilters>` :
    Assemble complexes with exactly defined sub-structures by applying advanced ligand filters for each binding site.

**Your role:** DART is designed for researchers in chemistry. While you have to think about the chemistry yourself, DART will allow you to explore the chemical space you want with the ligands you need.

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

- **Advanced Complex Control**: DART supports the assembly of a wide range of octahedral and square-planar geometries. Want more? Leave us a `feature request! <https://github.com/CCEMGroupTCD/DART/issues>`_

- **Intuitive Interface**: Users interact with DART via simple configuration files, no Python required.

- **Open-Source**: DART is an open-source tool developed in Python, making it entirely free to use.

Getting Started
===============

Are you ready to get started? :ref:`Install DART <installation_guide>`, explore our hands-on :ref:`quickstart guide <quickstart>` or read through the :ref:`ideas and concepts <dart_workflow>` behind DART.

.. toctree::
   :maxdepth: 2
   :caption: How to Use DART

   ./doc_files/part1_getting_started_and_how_to_guides/installation_guide
   ./doc_files/part1_getting_started_and_how_to_guides/examples/simple_octahedral_assembly/simple_octahedral_assembly
   ./doc_files/part1_getting_started_and_how_to_guides/tutorials_and_examples
   ./doc_files/part1_getting_started_and_how_to_guides/understanding_the_dart_workflow

.. toctree::
   :maxdepth: 2
   :caption: DART Modules

   ./doc_files/part2_dart_modules_documentation/module_overview
   ./doc_files/part2_dart_modules_documentation/metalig/metalig_database
   ./doc_files/part2_dart_modules_documentation/ligand_filters/ligand_filters
   ./doc_files/part2_dart_modules_documentation/assembler/assembler
.. todo: ./doc_files/part2_dart_modules_documentation/api/api_documentation



.. toctree::
   :maxdepth: 1
   :caption: Additional Resources

   ./doc_files/part3_additional_resources/troubleshooting_and_faqs
   ./doc_files/part3_additional_resources/tips_and_tricks
   ./doc_files/part3_additional_resources/current_limitations
   ./doc_files/part3_additional_resources/how_to_cite_dart
   ./doc_files/part3_additional_resources/version_history



