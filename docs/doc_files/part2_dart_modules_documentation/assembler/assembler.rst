.. _assembly_input:

Assembler Module
====================

.. contents:: :local:


Assembler Input
""""""""""""""""


The DART Assembler Module is designed for generating complexes from a predefined ligand database. It operates by either randomly or systematically selecting ligands and assembling them into a transition metal complex geometry.

The assembler module is run from the command line by providing a single configuration file:

.. code-block::

    DARTassembler assembler --path assembler_input.yml

**Assembler Input File:**

    The input for the assembler module is structured as a YAML file, consisting of:

    - **Global Options** include settings like output path, verbosity, and complex name length. These settings apply to all batches in the file.

    - **Batch Options** specify the details of each complex generation batch, such as metal center, oxidation state, geometry, and ligand database paths. Each batch must be named, and various parameters can be set to control the assembly process.

        - The first batch option needs a hyphen ('-') in front of it to mark the start of the batch.
        - Each batch option is required, but some can be left empty
        - Multiple batches can be specified and will be assembled in sequence.


**Copy-Paste Template:**

.. code-block::

    output_path: assembly/data_output     # Output folder path.
    overwrite_output: false               # Whether to overwrite the output if it already exists. Recommended: false.
    ffmovie: true                         # Whether to output a movie of the optimization process. Set to false to save disk space.
    concatenate_xyz: true                 # Whether to concatenate the xyz files of the optimization process.
    verbosity: 2                          # Output verbosity level (0-3), recommended: 2.
    complex_name_length: 8                # Length for generated complex names, recommended: 8.

    batches:                              # List of batches to generate.
      - name: First_batch                 # Name of the batch.
        metal_center: Fe                  # Chemical symbol of the desired metal center.
        metal_oxidation_state: 2          # Oxidation state of the desired metal center.
        total_charge: 0                   # Total charge of the complex.
        geometry: 2-1-1                   # Geometry of the complexes. :options: 2-1-1, 2-2, 3-2-1, 4-1-1, 5-1
        ligand_db_paths: [bidentate_ligands.jsonlines, monodentate_ligands.jsonlines, same_ligand_as_previous] # Paths to the ligand databases. Either single path or list of [path, 'same_ligand_as_previous'].
        ligand_choice: random             # How to choose the ligands. :options: random, all
        max_num_complexes: 100            # Maximum number of complexes to generate.
        forcefield: true                  # Whether to optimize the structures after generation with a force field.
        isomers: lowest_energy            # Which isomers to generate. :options: lowest_energy, all
        bidentate_rotator: auto           # How to rotate bidentate ligands in square-planar complexes. :options: horseshoe, slab, auto
        geometry_modifier_filepath:       # Path to the geometry modifier file. If not given, no geometry modification is performed.
        random_seed: 0                    # Random seed for the generation of the complexes.
        complex_name_appendix:            # String to append to the randomly generated complex name.

      - name: Second_batch                # Here comes the second batch.
        metal_center: Pd                  # Chemical symbol of the desired metal center.
        metal_oxidation_state: 3          # Oxidation state of the desired metal center.
        total_charge: 1                   # Total charge of the complex.
        geometry: 3-2-1                   # Geometry for the second batch. :options: 2-1-1, 2-2, 3-2-1, 4-1-1, 5-1
        ligand_db_paths: [tridentate_ligands.jsonlines, bidentate_ligands.jsonlines, monodentate_ligands.jsonlines] # Path to the ligand database. Either single path or list of [path, 'same_ligand_as_previous'].
        ligand_choice: all                # How to choose the ligands. :options: random, all
        max_num_complexes: 100            # Maximum number of complexes to generate.
        forcefield: false                 # Whether to optimize the structures after generation with a force field.
        isomers: all                      # Which isomers to generate. :options: lowest_energy, all
        bidentate_rotator: horseshoe      # How to rotate bidentate ligands in square-planar complexes. :options: horseshoe, slab, auto
        geometry_modifier_filepath: geometry_modifier.xyz      # Path to the geometry modifier file. If empty, no geometry modification is performed.
        random_seed: 1                    # Random seed for the generation of the complexes.
        complex_name_appendix: _2nd_batch # String to append to the randomly generated complex name.



Global Options
-----------------------------------

The assembler requires specifying global options at the file's beginning. These settings apply universally across all batches.

.. confval:: output_path

    :options: `dirpath`

    Path to directory in which the output will be saved.

.. confval:: overwrite_output

    :options: ``true``, ``false`` (recommended: ``false``)

    Whether to overwrite the output if it already exists.

.. confval:: ffmovie

    :options: ``true``, ``false``

    Whether to output a movie of the optimization process. Set to false to save disk space.

.. confval:: concatenate_xyz

    :options: ``true``, ``false``

    Whether to concatenate the xyz files of the optimization process.

.. confval:: verbosity

    :options: ``0``, ``1``, ``2``, ``3`` (recommended: ``2``)

    How much output to print (except the progress bars, which are always printed). ``0`` means only errors, ``1`` means also warnings, ``2`` means also normal info, ``3`` means also debug info.

.. confval:: complex_name_length

    :options: `integer > 0`

    Length of the randomly generated complex name (e.g. 'ZUMUVAMI'). Recommended: ``8``.

Batch Options
-----------------------------------

Each batch in the assembler input file defines configurations for generating complexes. The first option in each batch must be preceded by a hyphen ('-') to denote the start. All batch options are required, but some can remain unspecified if defaults are to be used. Multiple batches can be defined and will be assembled sequentially.

.. confval:: name

    :options: `string`

    Unique name for the batch for easy identification.

.. confval:: geometry

    :options: ``3-2-1``, ``4-1-1``, ``5-1``, ``2-1-1``, ``2-2``

    The geometry specifies the denticities of the ligands around the complex. For example, ``3-2-1`` would generate a complex with one tridentate, one bidentate and one monodentate ligand. Currently, the following topologies are supported:

        - **Octahedral complexes:** ``3-2-1``, ``4-1-1``, ``5-1``
        - **Square-planar complexes:** ``2-1-1``, ``2-2``

.. confval:: ligand_db_paths

    :options: `empty` OR `filepath` OR list(`filepath / keyword` )

    Specifies the source databases for ligands used in complex assembly. This option can be configured in two ways:

    - **List of Filepaths and/or Keywords:**
      A list where each entry is either a path to a ligand database file or the keyword ``same_ligand_as_previous``. The list should match the number of ligand sites as defined in the :confval:`geometry` option. For instance, in a ``3-2-1`` geometry, the first database in the list supplies tridentate ligands, the second supplies bidentate, and the third supplies monodentate ligands. The ``same_ligand_as_previous`` keyword can be used in place of a path to indicate that the ligand for the current site should be identical to the one used in the previous site for each assembled complex. This feature is useful for creating complexes with symmetrical or repeating ligand structures.

    - **Single Filepath or Empty:** 
      When a single path is provided, ligands for all sites will be drawn from this database. Identical to specifying a list with the same ligand db path for each ligand site. If empty, the entire :ref:`MetaLig database <metalig>` will be used.

    Note: Ligands in the database with a denticity not matching the specified :confval:`geometry` will be ignored during the assembly process. This ensures that only compatible ligands are selected for complex formation.

.. confval:: ligand_choice

    :options: ``random``, ``all``

    If ``random``, ligands will be chosen at random from the ligand database to assemble complexes. If ``all``, every possible combination of ligands will systematically be assembled. Note that this can lead to a very large number of complexes. The option :confval:`max_num_complexes` is used either way to limit the number of complexes generated.

.. confval:: max_num_complexes

    :options: `integer > 0`

    Maximum number of complexes to generate. Notes:

    - If :confval:`isomers` is set to ``all``, each isomer is counted as different complex. Note that the actual number of complexes generated can be a little higher in this case because for the last complex, all isomers are saved, even if this exceeds :confval:`max_num_complexes`.

    - If :confval:`ligand_choice` is set to ``all``, :confval:`max_num_complexes` will still be used to limit the number of complexes generated. To stop only after every possible complex is generated, set :confval:`max_num_complexes` to a very large number.

.. confval:: metal_center

    :options: `chemical symbol`

    Chemical symbol of the desired metal center, e.g. ``Pd`` or ``Fe``.

.. confval:: metal_oxidation_state

    :options: `integer > 0`

    Oxidation state of the desired metal center, e.g. ``3``.

.. confval:: total_charge

    :options: `integer`

    Total charge of the complex. Can be positive, negative or zero.

.. confval:: forcefield

    :options: ``true``, ``false``

    Whether to relax the generated structures with a force field before the post-assembly filters. Currently, the only available force field is the Universal Force Field (UFF).

.. confval:: isomers

    :options: ``lowest_energy``, ``all``

    The assembler will always generate all possible isomers. The option :confval:`isomers` determines which isomers are saved. If ``lowest_energy``, only the lowest energy isomer is saved as determined by a UFF forcefield. If ``all``, all isomers are saved.

.. confval:: bidentate_rotator

    :options: ``auto``, ``horseshoe``, ``slab``

    How to assemble bidentate ligands in square-planar complexes. Effects only the topologies ``2-2`` or ``2-1-1``. ``horseshoe`` and ``slab`` are the shapes of the underlying potential energy surfaces. ``horseshoe`` works best for ligands with a planar metallacycle, while non-planar ligands often give better results with ``slab``. ``auto`` will choose the shape automatically based on the ligand geometry.

    Tip: This option can severely affect the quality of generated complexes and how many make it through the post-assembly filter. For serious applications we recommend to set :confval:`max_num_complexes` to ``100``, try all three options and check how many complexes fail the post-assembly filter for each option (this info is returned at the end of the assembly if :confval:`verbosity` >= ``2``). Whichever option has the least complexes failing the post-assembly filter usually gives the highest quality complexes.

.. confval:: geometry_modifier_filepath

    :options: `empty` OR `filepath`

    Path to the geometry modifier file. If left empty, no geometry modification is performed.

    The geometry modifier file allows very advanced and fine-grained control over the geometry of the generated complexes. Usually it is not needed, since a forcefield optimization will often be a better option. However, there might be cases where it is desired to move atoms in an assembled ligand from one position to another position for all complexes with this ligand. This can be achieved with the geometry modifier file as shown in the Pd/Ni cross coupling example.

    For moving an atom to another position you need to supply the chemical symbol and the coordinates of the original atom and the coordinates the atom at it's new coordinates. The geometry modifier file is an .xyz file with two sets of atoms: The first set is all atoms that should be moved, the second set is the new positions of these atoms. Both sets of atoms are provided as "molecule" in the .xyz format and concatenated. The order and the chemical elements of both sets of atoms have to match up. In the assembly, for each generated complex, the atoms with coordinates in the first set are moved to the coordinates in the second set.

.. confval:: random_seed

    :options: `integer`

    Sets a seed for the random number generator to make the assembly of complexes exactly reproducible for each individual batch.

.. confval:: complex_name_appendix

    :options: `empty` or `string`

    Appends a custom string to the randomly generated name of each assembled complex. For example, if the appendix is set to ``_charge1``, a generated complex will be named 'ZUMUVAMI_charge1' if otherwise it would have been named 'ZUMUVAMI'.



Table: overview of all available assembler options
---------------------------------------------------

.. csv-table::
   :file: assembly_options_overview.csv
   :widths: 20, 10, 20, 50
   :header-rows: 1



.. include:: DART_assembly_output_explanation.rst