=============================================================
DART: Directed Assembly of Random Transition Metal Complexes
=============================================================

Welcome to the DART platform, a cutting-edge suite of tools for the generation and exploration of mono-metallic transition metal complexes! Developed by the CCEM group at Trinity College Dublin, DART is engineered to facilitate the design and analysis of molecular complexes for chemistry research.

DART integrates a collection of modules, each serving a unique function in the assembly process:

- **MetaLig** :
    Explore the comprehensive MetaLig database with 41,018 ligands extracted from the Cambridge Structural Database, complete with high-quality formal charge assignments.

- **Assembler** :
    Assemble novel transition metal complexes in a matter of seconds, guided by a simple configuration file for precise control over the resulting structures.

- **Ligand Filters** :
    Tailor assembled ligands to your research needs with a wide range of chemical and data-driven filters.

.. tip::

    Using DART is simple. After download, just run the DART assembler and start generating complexes by executing the following command in your terminal:

    .. code-block:: bash

        DARTassembler assembler --path assembly_input.yml


**Your role:** While DART is a powerful new tool for chemistry researchers to streamline the exploration of novel molecular complexes, it is not a substitute for the expertise of a practiced chemist. Therefore, DART facilitates the inspection of generated complexes by providing various output formats, empowering the expert analysis and evaluation that is the hallmark of a skilled chemist.

.. figure:: /_static/homepage_picture.png
   :width: 70%
   :align: center

   DART workflow for generating molecular complexes.

Who Should Use DART?
=====================

- **Experimental Chemists**: Get suggestions for new complexes with exactly defined chemical attributes.
- **Computational Researchers**: Generate novel transition metal complexes for high-throughput screening.
- **Software Developers**: Integrate DART into your chemistry workflow with our open-source code.

Why Use DART?
=============

Proposing new transition metal complexes can be a meticulous and lengthy endeavor. DART addresses this challenge by providing a streamlined, automated solution:

- **Powerful Assembler**: Generate novel transition metal complexes in a matter of seconds with DART's automated processes and comprehensive user options.

- **Advanced Complex Control**: DART supports the assembly of a wide range of octahedral and square-planar geometries. Want more? Leave us a `feature request! <https://github.com/CCEMGroupTCD/DART/issues>`_

- **Comprehensive Ligand Database**: Tap into a diverse array of 41,018 ligands from real-world structures for innovative materials discovery.

- **Sophisticated Filters**: Design specific complex types with advanced ligand filtering capabilities, targeting exactly defined chemical attributes.

- **Intuitive Interface**: Even without programming expertise, users can engage with DART through simple configuration files.

- **High Synthetic Viability**: DART supports the generation of complexes with high synthetic viability, thanks to its data-driven ligand filters based on actual crystal structures from the Cambridge Structural Database.

- **Open-Source**: DART is an open-source tool developed in Python, making it entirely free to use and customizable for your research.

Getting Started
===============

Eager to begin? Install DART or explore our tutorials for a hands-on demonstration of its capabilities.

.. toctree::
   :maxdepth: 1
   :caption: How to Use DART

   ./doc_files/part1_getting_started_and_how_to_guides/installation_guide
   ./doc_files/part1_getting_started_and_how_to_guides/examples/simple_octahedral_assembly/simple_octahedral_assembly
   ./doc_files/part1_getting_started_and_how_to_guides/understanding_the_dart_workflow
   ./doc_files/part1_getting_started_and_how_to_guides/tutorials_and_examples

.. toctree::
   :maxdepth: 1
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



