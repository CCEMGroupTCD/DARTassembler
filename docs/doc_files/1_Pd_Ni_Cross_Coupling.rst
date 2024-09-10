.. _pd_ni_cross_coupling:

Advanced Example
---------------------------------------

.. contents:: :local:

In this example we will introduce some advanced features of DART with a case study of the Pd/Ni catalyzed C-C cross coupling reaction. You will be guided step-by-step through the generation of Pd(II) and Ni(II) square-planar complexes which represent the oxidative addition intermediate, a key step in transition metal catalyzed cross coupling. Users new to DART should first review the :ref:`quickstart`. The following capabilities will be highlighted:

  - Filter ligands for each binding site individually
  - Fix a certain ligand to always be present in each complex
  - Generate complexes with all possible combinations of ligands
  - Generate all possible geometric isomers of a complex
  - Choose the best rotator for bidentate ligands
  - Systematically modify and customize the geometry of selected ligands

Utilizing just four commands, we will generate neutral complexes in a square-planar ``2-1-1`` geometry in which the monodentate ligands are fixed to be a phenyl (Ph) group and a Br substrate while the bidentate ligand iterates through 173 possible P,N-donors:

.. code-block:: bash

    DARTassembler ligandfilters --path ligandfilters_phenyl.yml
    DARTassembler ligandfilters --path ligandfilters_Br.yml
    DARTassembler ligandfilters --path ligandfilters_P_N_ligands.yml
    DARTassembler assembler --path Pd_Ni_assembler.yml

These commands expect four different configuration files which we will create and explain in the following sections.

Plan the Workflow
^^^^^^^^^^^^^^^^^^^^^^

To assemble the intermediates with DART, we will follow these steps:

**1. Ligand Selection**: We will first create a separate ligand set for each binding site so that we can control where different ligands coordinate in the square-planar ``2-1-1`` geometry. The first set contains a single Br ligand, the second one a single Ph ligand, and the third one all suitable P,N-donor ligands.

**2. Assembler Configuration**:

   - Assemble all possible ligand combinations.
   - Fix Ph and Br ligands within each complex, varying only the P,N-donor ligands.
   - Instruct DART to do a user-defined rotation of the Ph ligand to prevent collisions with the P,N-donor ligands.
   - Generate all possible geometric isomers of each complex.

Target Ligands via Ligand Filters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To initiate the assembly process within DART, we will create three distinct ligand databases, each corresponding to a specific binding site of our square-planar ``2-1-1`` complexes:

1. ``ligand_db_Br.jsonlines`` : containing 1 ligand (Br).
2. ``ligand_db_phenyl.jsonlines`` : containing 1 ligand (Ph).
3. ``ligand_db_P_N_donors.jsonlines`` : containing 173 neutral P,N-donor ligands.


**Creating the Br Dataset**

Let's start with creating the Br database. First, make a new file called ``ligandfilters_Br.yml`` with the following content:

.. code-block:: yaml

    # File: ligandfilters_Br.yml

    input_db_file: metalig
    output_db_file: ligand_db_Br.jsonlines
    output_ligands_info: False

    filters:
      - filter: denticities
        denticities: [1]

      - filter: ligand_composition
        elements: Br
        instruction: must_contain_and_only_contain
        apply_to_denticities:

Now let's run the ligand filters to obtain the file ``ligand_db_Br.jsonlines``. This process takes around 1-2 minutes since the entire MetaLig database with 41,018 ligands has to be read in:

.. code-block:: bash

    DARTassembler ligandfilters --path ligandfilters_Br.yml

**Creating the Ph Dataset**

The Ph ligand database is created in the same way. The input file ``ligandfilters_phenyl.yml`` is as follows:

.. code-block:: yaml

    # File: ligandfilters_phenyl.yml

    input_db_file: metalig
    output_db_file: ligand_db_phenyl.jsonlines
    output_ligands_info: False

    filters:
      - filter: denticities
        denticities: [1]

      - filter: ligand_composition
        elements: C6H5
        instruction: must_contain_and_only_contain
        apply_to_denticities:

To get the Ph ligand database, run:

.. code-block:: bash

    DARTassembler ligandfilters --path ligandfilters_phenyl.yml

Because we filter the Ph ligand simply by its composition, it would be possible that there are other monodentate ligands with the same composition. Yet, we can see from the printed output of the ligand filters that this is not the case:

.. code-block:: bash

    ===========   TOTAL   ===========
    Before filtering:  41018 ligands
    Filtered out:      41017 ligands
    Passed:            1 ligands
    Denticities:       1: 1
    Passed ligands:    C6H5

If there would be other ligands that we don't want in the database, you could simply add more filters. You can also pinpoint individual ligands by specifying their :ref:`Graph ID <filter_graph_IDs>` or :ref:`write customized filters using simple Python code <metalig_python_filtering>`.

**Creating the P,N-donor Ligand Dataset**

Finally, let's compile the dataset for our neutral P,N-donors. Please create the file ``ligandfilters_P_N_ligands.yml`` with the following content:

.. code-block:: yaml

    # File: ligandfilters_P_N_ligands.yml

    input_db_file: metalig
    output_db_file: ligand_db_P_N_donors.jsonlines
    output_ligands_info: False

    filters:
        # Keep only bidentate ligands
      - filter: denticities
        denticities: [2]

        # Keep only neutral ligands
      - filter: ligand_charges
        ligand_charges: [0]
        apply_to_denticities:

        # Keep only P,N-donors
      - filter: coordinating_atoms_composition
        elements: [P, N]
        instruction: must_contain_and_only_contain
        apply_to_denticities:

        # Keep only ligands that contain C, H, N, and P
      - filter: ligand_composition
        elements: [C, H, N, P]
        instruction: must_only_contain_in_any_amount
        apply_to_denticities:

        # Keep only ligands that have a history of binding to Pd or Ni
      - filter: metal_ligand_binding_history
        metal_ligand_binding_history: [Pd, Ni]
        apply_to_denticities:

        # Remove haptic ligands because DART cannot assemble those yet and skips them.
      - filter: remove_ligands_with_adjacent_coordinating_atoms
        remove_ligands_with_adjacent_coordinating_atoms: True


To generate the P,N-donor ligand database, run:

.. code-block:: bash

    DARTassembler ligandfilters --path ligandfilters_P_N_ligands.yml

To check the ligands in this database, you can run the ``dbinfo`` command to get a concatenated .xyz file of all ligand structures and a .csv file with information on all ligands:

.. code-block:: bash

    DARTassembler dbinfo --path ligand_db_P_N_donors.jsonlines

To browse through these ligands, run ``ase gui concat_ligand_db_P_N_donors.xyz``. Due to the filtering process, all 173 ligands are neutral P,N-donors.


Assemble the Intermediates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now that we have generated the ligand datasets, we can assemble the Pd/Ni square-planar complexes. To use the :ref:`Assembler Module <assembler>` we create a new input file called ``Pd_Ni_assembler.yml``. This file specifies the various input instructions necessary to generate neutral square-planar Pd/Ni(II) complexes with one Br, one Ph, and one P,N-donor ligand:

.. code-block:: yaml

    # File: Pd_Ni_assembler.yml

    output_directory: assembler_output
    batches:
      - name: Pd
        metal_center: Pd
        metal_oxidation_state: 2
        total_charge: 0
        geometry: 2-1-1
        ligand_db_file: [ligand_db_P_N_donors.jsonlines, ligand_db_Br.jsonlines, ligand_db_phenyl.jsonlines]
        max_num_complexes: all
        isomers: all
        random_seed: 0
        complex_name_appendix: _Pd

      - name: Ni
        metal_center: Ni
        metal_oxidation_state: 2
        total_charge: 0
        geometry: 2-1-1
        ligand_db_file: [ligand_db_P_N_donors.jsonlines, ligand_db_Br.jsonlines, ligand_db_phenyl.jsonlines]
        max_num_complexes: all
        isomers: all
        random_seed: 0
        complex_name_appendix: _Ni

Let us go through the relevant options:

1. All combinatorially possible ligand combinations will be assembled because ``max_num_complexes`` = ``all``.

2. Both geometric isomers for each complex will be generated by ``isomers`` = ``all``.

3. Most importantly, ``ligand_db_file`` specifies a list of three different ligand databases, one for each binding site in the ``2-1-1`` geometry. The list of databases in ``ligand_db_file`` has to be in the same order as the denticities in ``geometry``. This instructs DART to create complexes in which the first binding site (bidentate) is populated with ligands from the first ligand database, the second binding site (monodentate) is populated with ligands from the second database etc. Since the Br and Ph databases each consist of just one ligand, this allows us to fix the Br and Ph ligands to always be present in the complex, while varying the P,N-donor ligands.

Now that we have configured the assembler, we can run it:

.. code-block:: bash

    DARTassembler assembler --path Pd_Ni_assembler.yml

This will generate a new folder ``assembler_output`` which contains the generated complexes. The output of the assembler module concludes with the following lines:

.. code-block::

    ============  Total summary of DART assembly  ============
      - 692 complexes tried, 396 complexes successfully assembled.
      - 296 complexes failed because of post-filters:
        - clashing ligands: 296
    DART Assembler output files saved to YOURPATH/DART_cross_coupling/assembler_output
    Total runtime for assembling 396 complexes: 0:00:38.598620
    Done! All complexes assembled. Exiting DART Assembler.

A total of 396 complexes were assembled successfully, while 296 complexes failed the post-filters because of clashing ligands. As always, you can browse through the successfully assembled complexes using ``ase gui concat_passed_complexes.xyz``. This will give you a good overview of the chemical space of the complexes you just assembled.

Optimize the Output Geometry
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

While DART can assemble a wide range of geometries with high success rates, there are always ways to improve the quality of the assembled complexes for selected ligands. Now we will learn how to optimize the generated geometries by tuning the assembly parameters in DART. For the most part, this is as simple as changing a single line in the input file.

As a measure for the quality of the structures, we will use the percentage of successfully assembled complexes. Intuitively, the more complexes that are successfully assembled without clashing ligands, the better the geometry. However, it is always a good idea to look through the assembled complexes to ensure that the geometries are chemically reasonable, which is why we put so much emphasis on the visualization of structures with ``ase gui``.

There are three settings that you can play with, documented in detail in the :ref:`Assembler Module <assembler>`:

1. ``forcefield``: Relax the complexes with a UFF forcefield. Because the UFF sometimes struggles to describe metals, the metal and the donor atoms are kept fixed.
2. ``bidentate_rotator``: Choose the best rotator for bidentate ligands. The default mode is ``auto``, which lets DART choose the best rotator automatically, but you can also directly specify either ``slab`` or ``horseshoe``.
3. ``geometry_modifier_filepath``: Manually shift atoms in the assembled complexes from one position to another in a semi-automated way. Very powerful but requires a little manual work from the user to specify the exact shift.

**Bidentate Rotator and Forcefield Optimization**

The first two settings, the forcefield and the bidentate rotator mode, are very simple to use since they are just a single line in the input file:

.. code-block:: yaml

    # Update file: Pd_Ni_assembler.yml

    batches:
      - name: Pd
        ...
        forcefield: true
        bidentate_rotator: slab

      - name: Ni
        ...
        forcefield: true
        bidentate_rotator: slab

It is very quick to try out which of these options gives the best results since each assembly run in this example takes less than a minute. To evaluate the effect of these settings, we conducted a simple experiment as a proof-of-concept in which we tried to optimize the number of successfully assembled complexes out of a maximum of 692 possible complexes:

.. csv-table::
    :header: "Bidentate Rotator", "Without Optimization", "With Forcefield"
    :widths: 33, 33, 33

    "auto", 396, 385
    "slab", **473**, 458
    "horseshoe", 342, 330

The results show that the ``slab`` rotator is the best choice for our P,N-donor ligands, with a maximum of 473 geometric isomers passing the post-assembly filters. The forcefield optimization had little effect and rather decreased the number of successfully assembled complexes. In general, we do not recommend to use the UFF forcefield since in most cases it does not seem to improve the output geometries, but it is an easy option for you to try out.

**Custom Rotation of the Ph Ligand**

The third option ``geometry_modifier_filepath`` is very powerful because it allows the user to automatically shift atoms in an assembled complex from one position to another. In our example, we want to rotate the Ph ligand a little in order to reduce clashing with the P,N-donor. To do this, we have to provide a concatenated .xyz file with exactly two Ph molecules at different locations, specifying origin and destination of the shift. In order to implement this, please create a new file called ``Pd_phenyl_geometry_modification.xyz`` with the following content:

.. code-block::

    11
    Origin of shift
    C       -1.37885822      -1.37885822       0.00000000
    C       -2.61882178      -1.15229953       0.58746920
    H       -2.81065540      -0.30966243       0.97820005
    C       -3.58540325      -2.16127881       0.60496110
    H       -4.42111716      -2.00722798       1.02876196
    C       -3.33994508      -3.37153613       0.01642353
    H       -3.99087572      -4.06123881       0.04508034
    C       -2.16468025      -3.57422256      -0.60075492
    H       -2.01324444      -4.39862572      -1.04742888
    C       -1.15617112      -2.60121467      -0.59980122
    H       -0.31967840      -2.78328598      -1.01291115
    11
    Destination of shift
    C       -1.35204078      -1.31992066      -0.48197010
    C       -1.63807002      -1.59921746      -1.81393823
    H       -1.15796809      -1.15933936      -2.50351594
    C       -2.63068217      -2.52632650      -2.14228996
    H       -2.83682226      -2.69439953      -3.05388981
    C       -3.31015444      -3.19615732      -1.16198741
    H       -3.99601236      -3.81199015      -1.38686634
    C       -2.99747104      -2.97010890       0.12423560
    H       -3.43865917      -3.46858118       0.80170842
    C       -2.03898543      -2.01594063       0.49125628
    H       -1.86191301      -1.84921651       1.41015686

If you check the origin molecule in this file with the position of the Ph ligand in the assembled Pd complexes, you will see that they are identical. For the destination molecule, we are providing the atomic positions such that the Ph ligand is rotated. The ase gui tool is very helpful for these kinds of manipulations of .xyz files. You can also see the rotation of the Ph ligand by running ``ase gui Pd_phenyl_geometry_modification.xyz``.

As a tip, the best way to create these files is to first assemble the complexes without any forcefield optimization to get the .xyz file of the assembled complexes. Then, extract the coordinates of the Ph ligand from any of the assembled complexes and save it as .xyz file. To get the destination coordinates, the ``ase gui`` tool is very handy to manipulate .xyz files. Just read in the origin .xyz file with ase, manipulate it and then save the new structure as another .xyz. Finally, append the origin and the destination molecules and save them as a single file. Importantly, please make sure that there are no empty lines in between the two molecules or at the end of the file. To check if the file is correct, just read it in using ``ase gui`` and you should get two frames with the same molecule, but at different positions.

On another note, it is not necessary to shift all atoms of a specific ligand. When shifting atoms, DART simply iterates through all atoms in the complex and shifts every atom with matching chemical element and coordinates from its origin coordinates to its destination coordinates.

In order to run DART with the rotated Ph ligand, we specify ``geometry_modifier_filepath`` = ``Pd_phenyl_geometry_modification.xyz`` in the Pd batch. One thing to keep in mind is that you cannot provide the same file for the Ni batch, because Ni has a different atomic radius than Pd, leading to slightly different cartesian coordinates of the Ph ligand when coordinated to Ni. Therefore, you  have to create a new file ``Ni_phenyl_geometry_modification.xyz`` with the following content:

.. code-block::

    11
    Origin of shift
    C       -1.30814755      -1.30814755       0.00000000
    C       -2.54811110      -1.08158885       0.58746920
    H       -2.73994472      -0.23895175       0.97820005
    C       -3.51469257      -2.09056813       0.60496110
    H       -4.35040648      -1.93651730       1.02876196
    C       -3.26923441      -3.30082546       0.01642353
    H       -3.92016504      -3.99052813       0.04508034
    C       -2.09396957      -3.50351188      -0.60075492
    H       -1.94253376      -4.32791504      -1.04742888
    C       -1.08546044      -2.53050399      -0.59980122
    H       -0.24896772      -2.71257530      -1.01291115
    11
    Destination of shift
    C  -1.2937473829822488  -1.2398794889499034  -0.459800567
    C  -1.601562203225019  -1.4938654694550033  -1.7919791134
    H  -1.1309417375171225  -1.0428217585005106  -2.480879663
    C  -2.60218398542082  -2.411889969420291  -2.121542410022
    H  -2.8231586547736818  -2.5626138747145775  -3.032700502
    C  -3.2683208079706407  -3.09759503004543  -1.14306436557
    H  -3.9596020908062854  -3.707191723557817  -1.3683280760
    C  -2.9347449854640226  -2.896087320001541  0.14198826503
    H  -3.3668350933519986  -3.405606940434718  0.81710224434
    C  -1.9675509241539992  -1.9516175659140573  0.5112421903
    H  -1.7755298364628027  -1.8022912010877177  1.4301324877

As before, these numbers are obtained by rotating the Ph ligand. We can now add the three settings for optimizing DART output structures to the assembler configuration file:

.. code-block:: yaml

    # Update file: Pd_Ni_assembler.yml

    batches:
      - name: Pd
        ...
        geometry_modifier_filepath: Pd_phenyl_geometry_modification.xyz

      - name: Ni
        ...
        geometry_modifier_filepath: Ni_phenyl_geometry_modification.xyz

As before, we have performed an experiment to evaluate the effect of a simple rotation of the Ph ligand. The number of successfully assembled complexes out of a maximum of 692 is shown in the table below:

.. csv-table::
    :header: "Bidentate Rotator", "Without Optimization", "With Forcefield", "With Geometry Modifier"
    :widths: 25, 25, 25, 25

    "auto", 396, 385, 506
    "slab", 473, 458, **620**
    "horseshoe", 342, 330, 432

Our results show that manual intervention via the ``geometry_modifier_filepath`` can significantly increase the success rate. Together with using ``slab`` as the ``bidentate_rotator``, it allows to generate a maximum of 620 geometric isomers. However, these results are dependent on which kind of ligands you work with in DART.

**Conclusion**

This example demonstrates how to use DART in advanced mode for assembling highly customized complexes. While DART's default settings provide very good results in most cases, DART enables users to try a range of options to further fine-tune their assembled complexes.




