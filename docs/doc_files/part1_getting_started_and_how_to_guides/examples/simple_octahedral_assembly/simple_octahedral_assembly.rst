.. _quickstart:

Quickstart Guide
=================================

Welcome to the Quickstart Guide for DART! As an introductory example, we will walk through the process of assembling 100 octahedral Pd(II) complexes with a neutral charge. Each complex will consists of one tridentate, one bidentate, and one monodentate ligand, which is referred to as a 3-2-1 topology in DART.

For now, we will use the entire MetaLig database with 41,018 ligands to assemble the complexes. Then, we will learn how to filter down the input ligands to target complexes suitable for your own field of research and to generate those that are more likely to form stable complexes.


Confirming DART Installation
----------------------------

Before starting, ensure DART is correctly installed and configured:

1. Open your terminal.
2. Type ``dart --help`` and press Enter.

This command should display information how to use DART. If it throws an error, please consult the :ref:`Troubleshooting section<troubleshooting>`.

Preparing for Assembly
----------------------

1. **Creating an Output Directory**: This directory will store all generated files from the assembly process.

   .. code-block:: bash

       mkdir DART_simple_example
       cd DART_simple_example

2. **Configuring the Assembly Input**: The file ``assembly_input.yml`` outlines the options for the assembler. Create this file and input the following contents. For more information on each parameter, please refer to the section :ref:`assembly_input`.

.. code-block::

    # File: assembly_input.yml

    output_path: DART_output              # Path to the output folder.
    overwrite_output: false               # Whether to overwrite the output if it already exists. Recommended: false.
    ffmovie: true                         # Whether to output a movie of the optimization process. Set to false to save disk space.
    concatenate_xyz: true                 # Whether to concatenate the xyz files of the optimization process.
    verbosity: 2                          # How much output to print. Options: (0, 1, 2, 3), recommended is 2.
    complex_name_length: 8                # Length of the complex name. Recommended: 8.

    batches:                              # List of batches to generate. The first option in each batch needs a hyphen ('-') in front of it to mark the start of the batch.
      - name: Octahedral_Pd(II)           # Name of the batch. Note the hyphen in front of it to mark the start of the batch.
        topology: 3-2-1                   # Topology of the complexes. Options: 2-1-1, 2-2, 3-2-1, 4-1-1, 5-1
        ligand_db_paths:                  # Path to the ligand database. Either single path or list of [path, 'same_ligand_as_previous'].
        ligand_choice: random             # How to choose the ligands. Options: random, all
        max_num_complexes: 100            # Maximum number of complexes to generate.
        metal_center: Pd                  # Chemical symbol of the desired metal center.
        metal_oxidation_state: 2          # Oxidation state of the desired metal center.
        total_charge: 0                   # Total charge of the complex.
        forcefield: true                  # Whether to optimize the structures after generation with a force field.
        isomers: lowest_energy            # Which isomers to generate. Options: lowest_energy, all
        bidentate_rotator: auto           # How to rotate the bidentate ligands. Options: horseshoe, slab, auto
        geometry_modifier_filepath:       # Path to the geometry modifier file. If not given, no geometry modification is performed.
        random_seed: 0                    # Random seed for the generation of the complexes.
        complex_name_appendix:            # String to append to the randomly generated complex name.

Running the Assembler
---------------------

To use the assembler, execute the following command in your terminal:

.. code-block:: bash

    dart assembler --path assembly_input.yml

This command starts the assembler module, which creates complexes based on the ``assembly_input.yml`` file's instructions. The assembler will print the progress to the terminal and after a few minutes save the output files in the ``DART_output`` folder. Upon examining the generated complexes, you'll notice some wild chemistry. Using the entire MetaLig database with 41,018 different ligands clearly results in many complexes that are not very likely to form stable complexes. In fact, Figure 1 shows that there are 14 different elements in only these 300 ligands. It also shows two examples of complexes assembled in this run. In the following section, we will learn how to filter the ligands to generate complexes with a more realistic chemistry.

.. figure:: DART_output/quickstart_figure1.png
   :width: 100%
   :align: center

   Figure 1: (left) a histogram of the elements in the 300 complexes generated with the above input file. Clearly, the chemistry contained in just these few ligands is quite wild. (right) 2 examples of the complexes generated with the above input file.

Targeting Chemical Space
------------------------

To achieve complexes with more realistic and stable chemistry targeted to your own field of research, it is essential to filter the ligands used for the assembler. The ligand filter module makes this easy:

.. code-block:: bash

    dart ligandfilters --path ligandfilter_input.yml

The file ligandfilter_input.yml contains all the filter options that we want to set. For example, let's generate complexes in which the monodentate is neutral and only composed of C, H and N. Both the bidentate and the tridentate should be composed only of C, H, N, O, P and S. The bidentate should exclusively be an N-N donor. Additionally, we want to keep the ligands relatively small and set an upper limit of 30 atoms per ligand. Finally, we restrict our ligands to those that have been observed coordinating to either Ni, Pd or Pt in the Cambridge Structural Database. This is helpful to increase the likelihood that our Pd complexes will be stable, since the ligands have precedent coordinating to a group 10 transition metal. Helpfully, the MetaLig database contains a lot of this extrinsic information.

The input file for these filters looks like this:

.. code-block::

    # File: ligandfilter_input.yml

    input_ligand_db_path:
    output_ligand_db_path: filtered_ligand_db.json

    filters:
      # Keep only monodentates which are neutral
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

      # Keep only monodentates, bidentates and tridentates, since others will be ignored anyway for our 3-2-1 complexes
      - filter: denticities
        denticities: [1, 2, 3]


After running the above command, the MetaLig database will be filtered down to 1,561 ligands that meet the above criteria, including 427 monodentates, 615 bidentates and 519 tridentates.

**Inspecting the ligand database with 'dbinfo':**

To view the filtered ligands in a table format, execute this command:

.. code-block:: bash

    dart dbinfo --path filtered_ligand_db.json

This will generate a .csv file listing all the ligands in ``filtered_ligand_db.json``, which you can review in Excel or a similar program to ensure they meet your specifications.

**Assembling Complexes with Targeted Chemistry:**

To redo the assembly using the refined ligand database, update the `ligand_db_paths` in the assembly input file to the path of your filtered database. Also, change the output directory to prevent overwriting previous results.

.. code-block::

    # update assembly_input.yml
    output_path: DART_output_targeted
    ...
    batches:
          ...
          ligand_db_paths: filtered_ligand_db.json
          ...

The assembler will now draw from the 1,423 ligands that have been filtered to match our criteria. The resulting complexes will have a more uniform chemistry, while still covering a wide chemical space within the defined parameters. A histogram of the elements and two example complexes are shown in Figure 2. This method is excellent for generating a diverse set of complexes with realistic and targeted chemical properties for your research.

.. figure:: DART_output_targeted/quickstart_figure2.png
   :width: 100%
   :align: center

   Figure 2: (left) a histogram of the elements in the 300 complexes generated with the above input file. The chemistry is now confined to the six organic elements we specified. (right) 2 examples of the complexes generated with the above input file for targeted complexes. Note the N-N bidentate ligand and the neutral monodentate ligand. (The shown complexes might be different if DART has been updated since this guide was written.)

Understanding the Output of the Assembler Module
------------------------------------------------

The `DART_output_targeted` directory holds all the output files from the assembly module. For an in-depth explanation of each file, see the :ref:`assembly_output` section.

Let's examine the complex named IKIDAMIG to understand the range of information provided:

**IKIDAMIG_structure.xyz:**
    This file describes the geometry of the complex, showcasing an octahedral configuration with a Pd center and three distinct ligands. The structure of ICIDAMIG is shown in Figure 2.

**IKIDAMIG_ligandinfo.csv:**
    .. .. csv-table:: IKIDAMIG_ligandinfo.csv
       :file: IKIDAMIG_ligandinfo.csv
       :widths: 10, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9
       :header-rows: 1

    This file offers a snapshot of the MetaLig database, detailing the ligands in the complex. It provides a quick reference for properties like stoichiometry, denticity, donor atoms, and formal charges.

    Additionally, it includes data from complexes in the Cambridge Structural Database (CSD) that incorporate these ligands. These extrinsic properties include the IDs of each complex, the number of occurrences and all metal centers the ligand was found with. This  information can guide ligand selection and synthesis efforts.

**IKIDAMIG_ffmovie.xyz:**
    The file shows the forcefield relaxation process for the complex, indicating minor adjustments from the initial DART assembly.

**IKIDAMIG_data.json:**
    This comprehensive file offers detailed data on the complex, like the molecular graph, in a format suitable for further processing with DART modules or other applications.


Explore Your Complexes
----------------------

After the assembly, the folder ``DART_output_targeted`` will contain a rich array of complexes, all adhering to the specified chemical parameters. This targeted approach allows for a deep dive into a specific chemical space, bringing forward potential research candidates. We encourage you to explore the output and use the `dbinfo` module to gain more insight into the ligand database.





