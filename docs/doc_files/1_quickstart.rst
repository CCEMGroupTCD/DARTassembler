.. _quickstart:

Quickstart Guide
=================================

.. contents:: :local:

Welcome to the quickstart guide for DART!

As an introductory example, we will walk through the process of assembling 100 octahedral Pd(II) complexes with neutral formal charge. Each complex will feature one `mer`-tridentate, one bidentate, and one monodentate ligand, which is referred to as ``mer-3-2-1`` geometry in DART.

DART is based around the :ref:`MetaLig database <metalig>`, featuring 41,018 ligands as source for the assembly of novel complexes. In this tutorial, we will first assemble complexes using a randomly selected subset of 1000 ligands, without targeting any particular chemical space. Then, we will learn how to filter down the input ligands in order to generate complexes targeted to your own field of research and to generate those that are more likely to form stable complexes.

Confirm DART Installation
----------------------------

Before starting, ensure DART is correctly installed and configured:

1. Open your terminal.
2. Type ``DARTassembler installtest --path`` and press Enter.

This command runs a quick test of DART. If it throws an error, please consult the :ref:`Troubleshooting section<troubleshooting>`.

Inspect the Ligand Database
-------------------------------

As the source of ligands we will specify ``test_metalig``, which is a subset of 1000 ligands from the large MetaLig database. The full MetaLig database can be specified with the keyword ``metalig`` and needs about 1 minute to be read in on a standard laptop.

To inspect the ``test_metalig`` ligand database, use the ``dbinfo`` module:

.. code-block:: bash

    DARTassembler dbinfo --path test_metalig

This will immediately save two files, a concatenated .xyz file to inspect the structures and a .csv file to inspect the properties of the ligands. You can visualize and browse through the structures of the ligands by typing ``ase gui concat_test1000_MetaLigDB_v1.0.0.xyz`` in your terminal. By opening the .csv file with a program like Excel, you can inspect the properties of each ligands such as stoichiometry, denticity, donor atoms, and formal charge.

Assemble Novel Complexes
--------------------------------

To use the :ref:`Assembler Module <assembler>`, we need to provide an input file which outlines all settings for the assembly. Please create a new file called ``assembler.yml`` and copy-paste the following settings:

.. code-block::

       # File: assembler.yml

       output_directory: DART_output         # Path to a directory for the output files.
       batches:                              # List of batches to generate.
         - name: Octahedral_Pd(II)           # Name of the batch.
           metal_center: Pd                  # Chemical symbol of the desired metal center.
           metal_oxidation_state: 2          # Oxidation state of the desired metal center.
           total_charge: 0                   # Total charge of the complex.
           geometry: mer-3-2-1               # Geometry of the complexes. Options: 2-1-1, 2-2, mer-3-2-1, mer-4-1-1, 5-1
           ligand_db_file: test_metalig      # Path to the ligand db file. Options: metalig, test_metalig, filepath or list of paths/keywords (see documentation).
           max_num_complexes: 100            # Maximum number of complexes/isomers to generate.
           isomers: lowest_energy            # Which isomers to generate. Options: lowest_energy, all
           random_seed: 0                    # Optional. Random seed for reproducibility of results. Choose any integer.

Now execute the following command in your terminal:

.. code-block:: bash

    DARTassembler assembler --path assembler.yml

You will see that the assembler module prints the progress to the terminal and after around 2 minutes saves the output files in the ``DART_output`` folder. One thing to notice from the output is that while DART assembled 100 complexes successfully, there were 226 complexes which had to be discarded due to clashing ligands. This is a common and understandable issue when assembling complexes with completely random ligands. Soon, we will show you how to filter ligands and you will see that this often leads to a much higher success rate in assembling complexes.

Let's go into the ``DART_output`` folder and examine the generated complexes. First, we can browse through all successfully assembled geometries by typing ``ase gui concat_passed_complexes.xyz``. This command is useful for a quick visual inspection of the complexes.

We can also inspect the properties of assembled complexes by opening the file ``info_table.csv`` with a program such as Excel. This file displays information on all complexes which DART tried to assemble, including the one which failed for various reasons. If we look at the column ``note``, we see that many complexes failed due to clashing ligands, which is a common issue when assembling complexes with completely random ligands.

The following complexes are randomly picked from the output of the assembler module. The complexes are named ``AMAPIHOJ``, ``YOZEPECO`` and ``HOKEHAYA``, these names are randomly generated by DART for each assembled complex.

.. tip::

    All complexes displayed in the DART documentation are interactive 3D structures and can be rotated, translated and zoomed, thanks to `moldoc <https://github.com/lukasturcani/moldoc>`_.

.. moldoc::

    from DARTassembler.src.docs.docs_utils import get_moldoc_molecule

    run = 'quickstart_without_filters'
    complex_name = 'AMAPIHOJ'

    moldoc_display_molecule = get_moldoc_molecule(run=run, complex_name=complex_name)

.. moldoc::

    from DARTassembler.src.docs.docs_utils import get_moldoc_molecule

    run = 'quickstart_without_filters'
    complex_name = 'YOZEPECO'

    moldoc_display_molecule = get_moldoc_molecule(run=run, complex_name=complex_name)

.. moldoc::

    from DARTassembler.src.docs.docs_utils import get_moldoc_molecule

    run = 'quickstart_without_filters'
    complex_name = 'HOKEHAYA'

    moldoc_display_molecule = get_moldoc_molecule(run=run, complex_name=complex_name)

As you can see, using the entire MetaLig database without any filters results in a very diverse chemical space. In the following section, we will learn how to filter the ligands to generate complexes with more chemically uniform structures.


Target Chemical Space
------------------------

You can achieve a more targeted exploration of TMC chemical space by employing the :ref:`Ligand Filters Module <ligandfilters>`. This module allows you to filter the MetaLig by providing an input file with configurations for each pre-implemented filter. For example, let's suppose an inorganic chemist wants to explore Pd(II) octahedral complexes that occupy a `mer`-3-2-1 geometry. Additionally, they are only interested in ligands which meet the following criteria:

  - Monodentate: Neutral charge, composed only of C, H and N
  - Bidentate: N-N donor, composed only of C, H, N, O
  - Tridentate: Composed only of C, H, N, O
  - All ligands should have

    - no donor atoms with haptic interactions
    - no CH\ :sub:`2` units
    - specified bond orders
    - less than 30 atoms
    - been observed to coordinate to Ni, Pd or Pt in the Cambridge Structural Database

The last filter does not specify physical properties, but it is very useful to increase the likelihood that our Pd complexes will be chemically viable, since the ligands have precedent coordinating to a metal center from the same group. Helpfully, the MetaLig database contains not only physical ligand properties but also statistical information from the Cambridge Structural Database.

The following file translates these requirements into a set of filters that DART can understand. Please create a new file called ``ligandfilters.yml`` and copy-paste the following filters:

.. code-block::

    # File: ligandfilters.yml

    input_db_file: test_metalig
    output_db_file: filtered_ligand_db.jsonlines

    filters:

      # Keep only monodentates, bidentates and tridentates
      - filter: denticities
        denticities: [1, 2, 3]

      # Keep only monodentates which are neutral. Other denticities will be ignored by this filter.
      - filter: ligand_charges
        ligand_charges: [0]
        apply_to_denticities: [1]

      # The bi- and tridentate should be composed only of C, H, N, O
      - filter: ligand_composition
        elements: [C, H, N, O]
        instruction: must_only_contain_in_any_amount
        apply_to_denticities: [2, 3]

      # The monodentate should be composed only of C, H, N
      - filter: ligand_composition
        elements: [C, H, N]
        instruction: must_only_contain_in_any_amount
        apply_to_denticities: [1]

      # The bidentate must be an N-N donor
      - filter: coordinating_atoms_composition
        elements: [N, N]
        instruction: must_contain_and_only_contain
        apply_to_denticities: [2]

      # Remove ligands with likely haptic interactions
      - filter: remove_ligands_with_adjacent_coordinating_atoms
        remove_ligands_with_adjacent_coordinating_atoms: true
        apply_to_denticities:

      # Remove ligands with CH2 units
      - filter: smarts
        smarts: '[C&H2]'
        should_contain: false
        include_metal: false
        apply_to_denticities:

      # Remove ligands with missing bond orders. Should be used together with the smarts filter.
      - filter: remove_ligands_with_missing_bond_orders
        remove_ligands_with_missing_bond_orders: true
        apply_to_denticities:

      # All ligands should be relatively small with less than 30 atoms
      - filter: number_of_atoms
        min:
        max: 30
        apply_to_denticities:

      # Only allow ligands which have been observed to coordinate to Ni, Pd or Pt
      - filter: metal_ligand_binding_history
        metal_ligand_binding_history: [Ni, Pd, Pt]
        apply_to_denticities:

Now, run the ligand filters module:

.. code-block:: bash

    DARTassembler ligandfilters --path ligandfilters.yml

You will see that the ``test_metalig`` database is filtered down to 81 ligands that meet the above criteria, including 43 monodentates, 28 bidentates and 10 tridentates. This is already quite an interesting chemical space of ligands, but remember we are working only with a small test set of ligands. If we would have used the entire MetaLig, the numbers would be much higher: 699 ligands with 264 monodentates, 311 bidentates and 124 tridentates.

The Ligand Filters Module outputs a new ligand database file (``filtered_ligand_db.jsonlines``) and a folder with additional information about the filtering process (``info_filtered_ligand_db``). By now, you probably expect to find a concatenated .xyz file to inspect the ligand structures and a .csv file to inspect the ligand properties. And of course you're right!

First, you can check that all passed ligands have no CH\ :sub:`2` units with ``ase gui concat_Passed.xyz`` in the directory ``info_filtered_ligand_db/concat_xyz``. Furthermore, you will find one concatenated .xyz file for each filter, containing all ligands which were filtered out in this step. This is very useful to make sure that the filters are working exactly as you intended. For example, you can check all ligands with CH\ :sub:`2` units that were filtered out in the file ``concat_Filter07_smarts.xyz``.

**Assembling Complexes with Targeted Chemical Space:**

Now, we will redo the assembly process with the refined ligand database. First, update the ``ligand_db_file`` in the ``assembler.yml`` file so that it specifies the path to your newly filtered database. Also, change the ``output_directory`` to prevent overwriting previous results.

.. code-block::

    # update assembler.yml
    output_directory: DART_output_targeted
    ...
    batches:
          ...
          ligand_db_file: filtered_ligand_db.jsonlines
          ...

Now, run the assembler module again:

.. code-block:: bash

    DARTassembler assembler --path assembler.yml

The assembler will now draw all it's ligands from the 81 ligands that match the criteria we specified earlier. The resulting complexes will have a more uniform chemistry, while still covering a wide chemical space within the defined boundaries. This method is excellent for generating a diverse set of complexes with well defined chemical properties for your research.

Below you can see three randomly picked complexes from the output of the assembler module with the filtered ligands. Their names are ``CIGIBENA``, ``ITUXUGEB`` and ``WOREDEKI``. In contrast to before, they show a more well defined chemical space. In the same way, you can rapidly generate complexes for your own field of research by editing the ligand filters configuration file.

.. moldoc::

    from DARTassembler.src.docs.docs_utils import get_moldoc_molecule

    run = 'quickstart_with_filters'
    complex_name = 'CIGIBENA'

    moldoc_display_molecule = get_moldoc_molecule(run=run, complex_name=complex_name)

.. moldoc::

    from DARTassembler.src.docs.docs_utils import get_moldoc_molecule

    run = 'quickstart_with_filters'
    complex_name = 'ITUXUGEB'

    moldoc_display_molecule = get_moldoc_molecule(run=run, complex_name=complex_name)


.. moldoc::

    from DARTassembler.src.docs.docs_utils import get_moldoc_molecule

    run = 'quickstart_with_filters'
    complex_name = 'WOREDEKI'

    moldoc_display_molecule = get_moldoc_molecule(run=run, complex_name=complex_name)


Understand the Output of the Assembler Module
------------------------------------------------

The ``DART_output_targeted`` directory holds all the output files from the assembly module. For an in-depth explanation of each file, see the :ref:`assembly_output` section. The assembled complexes can be found in ``batches/Octahedral_Pd(II)/complexes``. Each complex is stored in a separate folder, named after the complex.

Let's examine the complex named ``CIGIBENA`` to understand the range of information provided:

**CIGIBENA_structure.xyz:**
    This file describes the geometry of the complex, showcasing an octahedral configuration with a Pd center and three distinct ligands. The file is formatted in XYZ, a simple and widely used format containing the 3D coordinates of all atoms.

**CIGIBENA_ligandinfo.csv:**
    This file provides a quick reference of all ligands in the complex, specifying properties like stoichiometry, denticity, donor atoms, and formal charge.

**CIGIBENA_data.json:**
    This comprehensive file offers detailed data on the complex, like structure, molecular graph and ligands, in a machine-readable format suitable for further processing.

Explore Your Complexes
----------------------------

The folder ``DART_output_targeted`` now contains a rich spectrum of complexes, all adhering to the parameters you specified earlier. This approach enables DART users to do a a deep dive into well-defined chemical spaces, bringing forward potentially interesting complexes for various applications. We encourage you to explore the DART output and discover the wealth of information it provides.

Want to learn more? Dive into a :ref:`case study using advanced DART features <pd_ni_cross_coupling>` or read more about the :ref:`ideas behind DART <dart_workflow>`.





