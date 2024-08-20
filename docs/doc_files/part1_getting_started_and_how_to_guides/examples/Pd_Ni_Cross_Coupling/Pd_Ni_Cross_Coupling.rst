.. _Pd_Ni_Cross_Coupling:

Pd/Ni Cross Coupling
--------------------

This example introduces advanced features of DART through the construction of Pd(II) and Ni(II) square-planar complexes as intermediates of the oxidative addition step of the transition metal catalyzed C-C cross-coupling reaction. Users new to DART should first review the :ref:`quickstart`. The following capabilities will be highlighted:

- Filter ligands for each binding site individually
- Fix a certain ligand to always be present in each complex
- Generate complexes with all possible combinations of ligands
- Shift individual atoms in a ligand to a specified position for fine-grained control
- Generate all possible isomers of a complex
- Choose the best rotator for bidentate ligands

Utilizing just four commands, we will generate neutral complexes in a ``2-1-1`` geometry in which the monodentate ligands are fixed to be a phenyl group and a bromine substrate while the bidentate ligand iterates through 173 possible P-N donors:

.. code-block:: bash

    DARTassembler ligandfilters --path input/ligandfilters_phenyl.yml
    DARTassembler ligandfilters --path input/ligandfilters_Br.yml
    DARTassembler ligandfilters --path input/ligandfilters_P_N_ligands.yml
    DARTassembler assembler --path input/Pd_Ni_assembler.yml

Subsequently, we present results from DFT calculations performed with Gaussian16 on all the assembled complexes, highlighting a diverse range of properties for these seemingly similar complexes.

Planning the Workflow
^^^^^^^^^^^^^^^^^^^^^^

To assemble the intermediates with DART, we will follow these steps:

1. **Ligand Selection**: We will produce three different ligand databases—just bromine, just phenyl, and all suitable P-N donor ligands—to precisely control the ligands at each binding site in the ``2-1-1`` geometry.

2. **Assembler Configuration**:

   - Assemble all possible ligand combinations.
   - Fix bromine and phenyl ligands within each complex, varying only the P-N donor ligands.
   - Implement a user-defined rotation of the phenyl ligand to prevent collisions with the P-N donor ligands.
   - Generate all possible isomers of each complex.

Preparation
^^^^^^^^^^^^^^^^^^^^^^

First, let's make a new directory called ``DART_cross_coupling`` and go there. Then, let's download all input files for this example from our GitHub repository so that we can explore the input files without worrying about mistakes:

.. code-block:: bash

    curl -O https://raw.githubusercontent.com/CCEMGroupTCD/DART/master/examples/Pd_Ni_Cross_Coupling/generate_complexes/input/ligandfilters_Br.yml
    curl -O https://raw.githubusercontent.com/CCEMGroupTCD/DART/master/examples/Pd_Ni_Cross_Coupling/generate_complexes/input/ligandfilters_phenyl.yml
    curl -O https://raw.githubusercontent.com/CCEMGroupTCD/DART/master/examples/Pd_Ni_Cross_Coupling/generate_complexes/input/ligandfilters_P_N_ligands.yml
    curl -O https://raw.githubusercontent.com/CCEMGroupTCD/DART/master/examples/Pd_Ni_Cross_Coupling/generate_complexes/input/Pd_Ni_assembler.yml
    curl -O https://raw.githubusercontent.com/CCEMGroupTCD/DART/master/examples/Pd_Ni_Cross_Coupling/generate_complexes/input/Ni_phenyl_geometry_modification.xyz
    curl -O https://raw.githubusercontent.com/CCEMGroupTCD/DART/master/examples/Pd_Ni_Cross_Coupling/generate_complexes/input/Pd_phenyl_geometry_modification.xyz


This directory includes six essential files:

- `ligandfilters_Br.yml` : filters to extract the bromine ligand
- `ligandfilters_phenyl.yml` : filters to extract the phenyl ligand
- `ligandfilters_P_N_ligands.yml` : filters to extract the P-N donor ligands
- `Pd_Ni_assembler.yml` : configuration file for the DART assembler module
- `Ni_phenyl_geometry_modification.xyz` : file specifying how to shift the phenyl ligand for the Ni complexes
- `Pd_phenyl_geometry_modification.xyz` : file specifying how to shift the phenyl ligand for the Pd complexes

We will go through the files one by one and explain what they do.

Running the Ligand Filters
^^^^^^^^^^^^^^^^^^^^^^^^^^

To initiate the assembly process within DART, you'll create three distinct ligand databases, each corresponding to a specific site of your complex:

1. **Bromine Database:** contains the bromine ligand.
2. **Phenyl Database:** contains the phenyl ligand.
3. **P-N Donor Database:** contains all the P-N donor ligands of interest.

Filtering is straightforward when you're dealing with unique stoichiometries like 'Br' and 'C6H5'. If you need to narrow down your selection further, you can pinpoint the desired ligands using their Graph ID.

**Setting Up the Phenyl Ligand Database**

Let's start with creating the phenyl database:

.. code-block:: bash

    DARTassembler ligandfilters --path input/ligandfilters_phenyl.yml

This command will output a new file named ``ligand_db_phenyl.jsonlines``. At the end of the filtering process, you'll see a summary like this:

.. code-block::

    Number of ligands before filtering: 41018
    Number of ligands filtered out: 41017
    Number of ligands after filtering: 1
    Number of ligands per denticity: 1: 1
      --> The selected ligand: C6H5

This confirms that your phenyl ligand database is now ready, containing just the one ligand you need.

**Isolating the Bromine Ligand**

To prepare the bromine database, the steps are similar:

.. code-block:: bash

    DARTassembler ligandfilters --path input/ligandfilters_Br.yml

Executing this will create the ``ligand_db_Br.jsonlines`` file, reserved for the bromine ligand.

**Creating the P-N Donor Ligand Database**

Finally, let's compile the database for your P-N donors:

.. code-block:: bash

    DARTassembler ligandfilters --path input/ligandfilters_P_N_ligands.yml

Upon completion, you'll have the ``ligand_db_P_N_donors.jsonlines`` file. Here's what the output will look like:

.. code-block::

    Number of ligands before filtering: 41018
    Number of ligands filtered out: 40845
    Number of ligands after filtering: 173
    Number of ligands per denticity: 2: 173

This indicates you have successfully filtered down to 173 bidentate ligands. To check the ligands, use:

.. code-block:: bash

    DARTassembler dbinfo --path ligand_db_P_N_donors.jsonlines

With these three ligand databases in hand, you're all set to move on to the assembly module.

Running the Assembler
^^^^^^^^^^^^^^^^^^^^^

The assembler module is configured by the file ``input/Pd_Ni_assembler.yml``. The documentation for all these options can be found at :ref:`assembler`. Let us go through the file and look at the important options:

.. code-block::

    verbosity: 2
    ffmovie: true
    concatenate_xyz: true
    overwrite_output: true
    output_directory: assembler_output
    complex_name_length: 8

    batches:
      - ligand_db_file: [ligand_db_P_N_donors.jsonlines, ligand_db_Br.jsonlines, ligand_db_phenyl.jsonlines]
        isomers: all
        max_num_complexes: 999999999999999999
        ligand_choice: all
        geometry_modifier_filepath: input/Pd_phenyl_geometry_modification.xyz
        bidentate_rotator: slab
        metal_center: Pd
        metal_oxidation_state: 2
        name: P_N_Donors_Pd_Metal_Centre
        forcefield: false
        random_seed: 0
        geometry: 2-1-1
        total_charge: 0
        complex_name_appendix: _PN_Pd


      - ligand_db_file: [ligand_db_P_N_donors.jsonlines, ligand_db_Br.jsonlines, ligand_db_phenyl.jsonlines]
        isomers: all
        max_num_complexes: 999999999999999999
        ligand_choice: all
        geometry_modifier_filepath: input/Ni_phenyl_geometry_modification.xyz
        bidentate_rotator: slab
        metal_center: Ni
        metal_oxidation_state: 2
        name: P_N_Donors_Ni_Metal_Centre
        forcefield: false
        random_seed: 0
        geometry: 2-1-1
        total_charge: 0
        complex_name_appendix: _PN_Ni


At the beginning we define global preferences, such as enabling forcefield trajectories and XYZ file concatenation. Then, two batches are set up, with identical options apart from the metal center. Let us go through the important options:

1. All combinatorial possible ligand combinations will be assembled (`ligand_choice: all`), and to prevent premature halting, `max_num_complexes` is set to a very high number.

2. The `geometry_modifier_filepath` and `bidentate_rotator` options are for advanced control over the assembly process and are explained below in :ref:`optimizing_geometry`.

3. To explore all isomeric forms, we opt for `isomers: all`.

4. Most importantly, `geometry: 2-1-1` instructs DART to craft complexes with one bidentate and two monodentate ligands, i.e. three different ligand sites. The three ligand databases in `ligand_db_file` are provided in the same order, so that the first ligand database is used for the first ligand site, the second ligand database for the second ligand site, and so on. This allows us to fix the bromine and phenyl ligands to always be present in the complex, while varying the P-N donor ligands.

Now that we have configured the assembler, we can run it:

.. code-block:: bash

    DARTassembler assembler --path input/Pd_Ni_assembler.yml

This will generate a new folder ``assembler_output`` which contains the generated complexes. To get an understanding of the output of the assembler module please refer the section :ref:`assembly_output`. The output of the assembler module concludes with the following lines:

.. code-block::

    ============  Total summary of DART assembly  ============
      - 692 complexes tried, 620 complexes successfully assembled.
      - 72 complexes failed because of post-filters:
        - clashing ligands: 72
    DART Assembler output files saved to your_path/assembler_output
    Total runtime for assembling 620 complexes: 0:00:51.236134
    Done! All complexes assembled. Exiting DART Assembler.

A total of 620 complexes were assembled successfully, while 72 complexes failed the post-filters because of clashing ligands. Figure 1 showcases a subset of the assembled complexes.


.. figure:: /_static/part1/examples/Pd_Ni_Cross_Coupling/assembled_complexes.png
   :width: 100%
   :align: center

   Figure 1: A selection of assembled Pd(II) complexes with fixed bromine and phenyl ligands and varying P-N donor ligands.


Exploring Molecular Properties with DART
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

DART excels in the automatic assembly of novel complexes, leveraging a diverse ligand pool to pave the way for innovative complex design. This approach shines in our example, where even a limited ligand selection shows a broad spectrum of complex properties. Here we show a quick overview of the results of DFT calculations for the P-N bite angle and the HOMO-LUMO gap for all 620 complexes, performed using Gaussian16 as detailed in our DART publication.

The data presented in Figure 2 underscores the extensive range of properties achievable by modifying even a single ligand type within the complexes.

.. figure:: /_static/part1/examples/Pd_Ni_Cross_Coupling/dft_figure_reduced.png
   :width: 100%
   :align: center

   Figure 2: DFT calculated properties of the assembled complexes.

.. _optimizing_geometry:

Optimizing the Output Geometry in DART
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Proper geometry optimization is crucial for the successful assembly of transition metal complexes. DART provides three key options to ensure optimal geometry: `forcefield`, `geometry_modifier_filepath`, and `bidentate_rotator`. These options are documented in detail in section :ref:`assembler`.

The `forcefield` option leverages a UFF forcefield to relax the output structure before going through the post-filter.

The `geometry_modifier_filepath` is a powerful tool for manual geometry correction. It's particularly useful when standard optimization methods fail to prevent certain ligand collisions, which might be the case with complex ligand structures or when specific orientations are required. With this option, users can input a file detailing the desired adjustments, and DART will reposition the atoms accordingly.

The `bidentate_rotator` setting controls the internal rotation of bidentate ligands. It can be set to `auto`, which allows DART to choose the rotation, or to specific modes (`slab` or `horseshoe`), giving users control over the bidentate ligand orientation.

To evaluate the efficacy of these optimization tools, we conducted an experiment focused on the assembly success rate—a key indicator of optimal geometry. The experiment involved multiple assembly runs, each varying the optimization method:

.. csv-table::
    :header: "Bidentate Rotator", "Without Optimization", "With Forcefield", "With Geometry Modifier"
    :widths: 25, 25, 25, 25

    "Auto", "68.4%", "55.2%", "72.7%"
    "Slab", "56.8%", "66.2%", "89.6%"
    "Horseshoe", "49.4%", "47.7%", "62.4%"

Our results show that manual intervention via the `geometry_modifier_filepath` significantly increases the success rate, particularly when the `slab` option is employed for the `bidentate_rotator`. However, these results are completely dependent on which kind of ligands to assemble. While DART's default settings provide satisfactory results for many cases, these tools offer valuable avenues for optimization, enhancing the likelihood of successful complex assembly.





