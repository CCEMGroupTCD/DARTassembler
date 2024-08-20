.. _quickstart:

Quickstart Guide
=================================

Welcome to the quickstart guide for DART!

As an introductory example, we will walk through the process of assembling 100 octahedral Pd(II) complexes with neutral charge. Each complex will feature one `mer`-tridentate, one bidentate, and one monodentate ligand, which is referred to as ``mer-3-2-1`` geometry in DART.


First, we will use the entire MetaLig database with 41,018 ligands as source of ligands for the assembled complexes, without any filters. Then, we will learn how to filter down the input ligands in order to generate complexes targeted to your own field of research and to generate those that are more likely to form stable complexes.

Confirming DART Installation
----------------------------

Before starting, ensure DART is correctly installed and configured:

1. Open your terminal.
2. Type ``DARTassembler --help`` and press Enter.

This command should display information how to use DART. If it throws an error, please consult the :ref:`Troubleshooting section<troubleshooting>`.

Inspecting the Ligand Database
-------------------------------

As the source of ligands, we will use the ``test_metalig`` for now, which is a subset of 1000 ligands from the quite large MetaLig database. This is just to speed up the assembly process for this quickstart guide and can be useful for your development phase as well. The full MetaLig database can be specified with the keyword ``metalig`` and needs about 2 minutes to be read in on a standard laptop.

To inspect the ``test_metalig`` ligand database, use the ``dbinfo`` module:

.. code-block:: bash

    DARTassembler dbinfo --path test_metalig

This will immediately save two files, a concatenated .xyz file to inspect the structures and a .csv file to inspect the properties of the ligands. You can visualize and browse through the structures of the ligands by typing ``ase gui concat_test1000_MetaLigDB_v1.0.0.xyz`` in your terminal. By opening the .csv file with a program like Excel, you can inspect the properties of each ligands such as stoichiometry, denticity, donor atoms, and formal charge.

Assembling Novel Complexes
--------------------------------

Using the :ref:`Assembler Module <assembler>` is as easy as executing the following command in your terminal:

.. code-block:: bash

    DARTassembler assembler --path assembler.yml

The file ``assembler.yml`` outlines the user-defined options for the assembly. Make a new file and copy in the following contents:

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

The command starts the assembler module, which creates complexes based on the instructions in the file ``assembler.yml``. The assembler will print the progress to the terminal and after around 2 minutes save the output files in the ``DART_output`` folder.

Let's go into the ``DART_output`` folder and examine the generated complexes. First, we can browse through all successfully assembled geometries by typing ``ase gui concat_passed_complexes.xyz``. This command is useful for a quick visual inspection of the complexes.

We can also inspect the properties of assembled complexes in the file ``info_table.csv``. This file displays information on all complexes which DART tried to assemble, including the one which failed for various reasons. If we look at the column ``note``, we see that many complexes failed due to clashing ligands, which is a common issue when assembling complexes with completely random ligands.

**Intermediary conclusion:** Using the entire MetaLig database without any filters clearly results in many complexes that are neither useful nor very likely to form stable complexes. In the following section, we will learn how to filter the ligands to generate complexes with a more realistic chemistry and with specific kinds of ligands.

Targeting Chemical Space
------------------------

To achieve complexes with more realistic and stable chemistry targeted to your own field of research, it is essential to filter the ligands used for the assembler. The :ref:`Ligand Filters Module <ligandfilters>` makes this easy:

.. code-block:: bash

    DARTassembler ligandfilters --path ligandfilters.yml

The file ``ligandfilters.yml`` contains all the filter options that we want to set. Let's stay with assembling octahedral Pd(II) complexes with a `mer`-3-2-1 geometry, but let's restrict the ligands used for each binding site:

- Monodentate: Neutral, composed only of C, H and N
- Bidentate: N-N donor, composed only of C, H, N, O, P and S
- Tridentate: Composed only of C, H, N, O, P and S
- All ligands should have less than 30 atoms
- All ligands should have been observed to coordinate to Ni, Pd or Pt in the Cambridge Structural Database

The last filter does not specify physical properties, but it is very useful to increase the likelihood that our Pd complexes will be stable, since the ligands have precedent coordinating to a group 10 transition metal. Helpfully, the MetaLig database contains this kind of statistical information from the CSD as well as 'normal' ligand properties.

The following file translates these requirements into a set of filters that DART can understand:

.. code-block::

    # File: ligandfilters.yml

    input_db_file: test_metalig
    output_db_file: filtered_ligand_db.jsonlines

    filters:

      # Keep only monodentates, bidentates and tridentates
      - filter: denticities
        denticities: [1, 2, 3]

      # Keep only monodentates which are neutral, ignore other denticities
      - filter: ligand_charges
        ligand_charges: [0]
        apply_to_denticities: [1]

      # For all ligands, keep only those with the following elements or subsets of these elements
      - filter: ligand_composition
        elements: [C, H, N, O, P, S]
        instruction: must_only_contain_in_any_amount
        apply_to_denticities:

      # Only the monodentate should be only composed of C, H, N though
      - filter: ligand_composition
        elements: [C, H, N]
        instruction: must_only_contain_in_any_amount
        apply_to_denticities: [1]

      # The bidentate must be an N-N donor
      - filter: coordinating_atoms_composition
        elements: [N, N]
        instruction: must_contain_and_only_contain
        apply_to_denticities: [2]

      # All ligands should be relatively small with less than 30 atoms
      - filter: number_of_atoms
        min:
        max: 30
        apply_to_denticities:

      # Only allow ligands which have been observed to coordinate to Ni, Pd or Pt
      - filter: metal_ligand_binding_history
        metal_ligand_binding_history: [Ni, Pd, Pt]
        apply_to_denticities:

After running the above command, the ``test_metalig`` database will be filtered down to 116 ligands that meet the above criteria, including 53 monodentates, 42 bidentates and 21 tridentates. If we would have used the entire MetaLig, the numbers would be much higher: 1,561 ligands with 427 monodentates, 615 bidentates and 519 tridentates.

The Ligand Filters module outputs a new ligand database file (``filtered_ligand_db.jsonlines``) and a folder with additional information about the filtering process (``info_filtered_ligand_db``). By now, you probably expect to find a concatenated .xyz file to inspect structures, and a .csv file to inspect properties. And of course you're right!

Additionally to the structures of the passed ligands (``concat_Passed.xyz``), you will find one concatenated .xyz file for each filter, containing all ligands which were filtered out in this step. This is very handy to make sure that the filters are working exactly as you intended.

**Assembling Complexes with Targeted Chemistry:**

Now, we will redo the assembly process with the refined ligand database. First, update the ``ligand_db_file`` in the ``assembler.yml`` file so that it specifies the path to your freshly filtered database. Also, change the output directory to prevent overwriting previous results.

.. code-block::

    # update assembler.yml
    output_directory: DART_output_targeted
    ...
    batches:
          ...
          ligand_db_file: filtered_ligand_db.jsonlines
          ...

The assembler will now draw all it's ligands from the 116 ligands that match the criteria we specified earlier. The resulting complexes will have a more uniform chemistry, while still covering a wide chemical space within the defined boundaries. This method is excellent for generating a diverse set of complexes with realistic and targeted chemical properties for your research.

.. figure:: DART_output_targeted/quickstart_figure2.png
   :width: 100%
   :align: center

   Figure 2: (left) a histogram of the elements in the 300 complexes generated with the above input file. The chemistry is now confined to the six organic elements we specified. (right) 2 examples of the complexes generated with the above input file for targeted complexes. Note the N-N bidentate ligand and the neutral monodentate ligand. (The shown complexes might be different if DART has been updated since this guide was written.)

Understanding the Output of the Assembler Module
------------------------------------------------

The ``DART_output_targeted`` directory holds all the output files from the assembly module. For an in-depth explanation of each file, see the :ref:`assembly_output` section. The assembled complexes can be found in ``batches/Octahedral_Pd(II)/complexes``. Each complex is stored in a separate folder, named after the complex.

Let's examine the complex named ADINOBUX to understand the range of information provided:

**ADINOBUX_structure.xyz:**
    This file describes the geometry of the complex, showcasing an octahedral configuration with a Pd center and three distinct ligands.

**ADINOBUX_ligandinfo.csv:**
    This file offers a snapshot of the MetaLig database, detailing the ligands in this complex. It provides a quick reference for properties like stoichiometry, denticity, donor atoms, and formal charge.

**ADINOBUX_data.json:**
    This comprehensive file offers detailed data on the complex, like structure, molecular graph and ligands, in a machine-readable format suitable for further processing.

Explore Your Complexes
----------------------

The folder ``DART_output_targeted`` now contains a rich spectrum of complexes, all adhering to the parameters you specified earlier. This approach enables DART users to do a a deep dive into well-defined chemical spaces, bringing forward potentially interesting complexes for various applications. We encourage you to explore the DART output and discover the wealth of information it provides.

Keen to learn more? Dive into a :ref:`case study using advanced DART features <Pd_Ni_Cross_Coupling>` or read more about the :ref:`DART philosophy <dart_workflow>`.





