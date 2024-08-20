.. _assembly_input:

Assembler Module
====================

.. contents:: :local:


Assembler Input
""""""""""""""""

The DART Assembler Module is designed for generating novel transition metal complexes by specifying a metal center and coordinating ligands from a ligand database.

The assembler module is run from the command line by providing a single configuration file:

.. code-block:: bash

    DARTassembler assembler --path assembler_input.yml

**Copy-Paste Template:**

.. code-block:: yaml

    ################## Settings for the DART assembler module. ##################
    # Everything after '#' is ignored by the program and only there for the user.

    output_directory: DART                # Path to a directory for the output files.
    batches:                              # List of batches to generate.
      - name: First_batch                 # Name of the batch.
        metal_center: Fe                  # Chemical symbol of the desired metal center.
        metal_oxidation_state: 2          # Oxidation state of the desired metal center.
        total_charge: 0                   # Total charge of the complex.
        geometry: 2-1-1                   # Geometry of the complexes. Options: `2-1-1`, `2-2`, `mer-3-2-1`, `mer-4-1-1`, `5-1`
        ligand_db_file: metalig           # Path to the ligand db file. Options: `metalig`, `test_metalig`, filepath or list of paths/keywords (see documentation).
        max_num_complexes: 100            # Maximum number of complexes/isomers to generate. Integer or `all`.
        isomers: all                      # Which isomers to save for each complex. Options: `lowest_energy`, `all`
        #random_seed: 0                   # Optional. Random seed for reproducibility of results. Choose any integer.
        #forcefield: false                # Optional. Whether to optimize the structures after generation with a UFF force field. Recommended: `false`.
        #bidentate_rotator: auto          # Optional. How to rotate bidentate ligands in square-planar complexes. Options: `auto`, `horseshoe`, `slab`. Recommended: `auto`.
        #geometry_modifier_filepath:      # Optional. Path to a geometry modifier file to shift atoms in complexes.
        #complex_name_appendix:           # Optional. String to append to each randomly generated complex name for labeling purposes.

    #ffmovie: false                       # Optional. Whether to output a movie (concatenated xyz file) of the forcefield optimization process.
    #concatenate_xyz: true                # Optional. Whether to save concatenated xyz files with all passed/failed complexes respectively.
    #verbosity: 2                         # Optional. Output verbosity level (0-3), recommended: `2`.
    #same_isomer_names: true              # Optional. Whether to give the same name to isomers of the same complex and then to number them.
    #complex_name_length: 8               # Optional. Length for generated complex names, recommended: 8.

Users can download this template into their current directory by running:

.. code-block:: bash

    DARTassembler configs --path .

Batch Options
-----------------------------------

Batch options are mostly mandatory and specify details concerning the metal center and the ligands. Multiple batches can be specified and will be assembled in sequence. For each batch, the first option must be preceded by a hyphen ('-') to denote the start of the list.

.. confval:: name

    :options: `string`
    :required: ``true``

    Unique name for the batch for easy identification. Each batch must have a different name.

.. confval:: metal_center

    :options: `chemical symbol`
    :required: ``true``

    Chemical symbol of the desired metal center, e.g. ``Pd`` or ``Fe``.

.. confval:: metal_oxidation_state

    :options: `integer > 0`
    :required: ``true``

    Oxidation state of the desired metal center, e.g. ``2``.

.. confval:: total_charge

    :options: `integer`
    :required: ``true``

    Total charge of the complex. Can be positive, negative or zero.

.. confval:: geometry

    :options: ``mer-3-2-1``, ``mer-4-1-1``, ``5-1``, ``2-1-1``, ``2-2``
    :required: ``true``

    The geometry specifies the denticities of the ligands around the complex. For example, ``mer-3-2-1`` would generate a complex with one `mer`-tridentate, one bidentate and one monodentate ligand. Currently, the following topologies are supported:

        - **Octahedral complexes:** ``mer-3-2-1``, ``mer-4-1-1``, ``5-1``
        - **Square planar complexes:** ``2-1-1``, ``2-2``

.. confval:: ligand_db_file

    :options: ``empty``/``metalig``, ``test_metalig``, `filepath` OR list(`filepath / keyword` )
    :required: ``false``
    :default: ``metalig``

    Specifies the source databases for ligands used in complex assembly. This option can be configured in two ways:

    - **List of Filepaths and/or Keywords:**
      A list where each entry is either a path to a ligand database file or the keyword ``same_ligand_as_previous``. The list should match the number of ligand sites as defined in the :confval:`geometry` option. For instance, in a ``mer-3-2-1`` geometry, the first database in the list supplies tridentate ligands, the second supplies bidentate, and the third supplies monodentate ligands. The ``same_ligand_as_previous`` keyword can be used in place of a path to indicate that the ligand for the current site should be identical to the one used in the previous site for each assembled complex. This feature is useful for creating complexes with symmetrical or repeating ligand structures.

    - **Single Filepath or Empty:** 
      When a single path is provided, ligands for all sites will be drawn from this database. Identical to specifying a list with the same ligand db path for each ligand site. If empty or ``metalig``, the entire :ref:`MetaLig database <metalig>` will be used. If ``test_metalig``, a small subset of the MetaLig database will be used to speed up testing.

    Note: Ligands in the database with a denticity not matching the specified :confval:`geometry` will be ignored during the assembly process. This ensures that only compatible ligands are selected for complex formation.

.. confval:: max_num_complexes

    :options: `integer > 0` OR ``all``
    :required: ``true``

    Maximum number of complexes to generate. If :confval:`max_num_complexes` is set to ``all``, it will generate all combinatorically possible complexes.

    Note: If :confval:`isomers` is set to ``all``, each isomer is counted as different complex. Note that the actual number of complexes generated can be a little higher in this case because for the last complex, all isomers are saved, even if this exceeds :confval:`max_num_complexes`.

.. confval:: isomers

    :options: ``lowest_energy``, ``all``
    :required: ``true``

    The assembler will always generate all possible isomers. The option :confval:`isomers` determines which isomers are saved. If ``lowest_energy``, only the lowest energy isomer is saved as determined by a UFF forcefield. If ``all``, all isomers are saved.

.. confval:: random_seed

    :options: `integer`
    :required: ``false``
    :default: Randomly chosen between 1000 and 9999

    Sets a seed for the random number generator to make the assembly of complexes exactly reproducible for each individual batch. If not set, a random seed between 1000 and 9999 is completely randomly chosen for each batch and recorded. That means the run is still reproducible by checking which random seed was chosen, but the seed is not known in advance.

.. confval:: forcefield

    :options: ``true``, ``false``
    :required: ``false``
    :default: ``false``

    Whether to relax the generated structures with a force field before the post-assembly filters. Currently, the only available force field is the Universal Force Field (UFF).

.. confval:: bidentate_rotator

    :options: ``auto``, ``horseshoe``, ``slab``
    :required: ``false``
    :default: ``auto``

    How to assemble bidentate ligands in square-planar complexes. Effects only the topologies ``2-2`` or ``2-1-1``. ``horseshoe`` and ``slab`` are the shapes of the underlying potential energy surfaces. ``horseshoe`` works best for ligands with a planar metallacycle, while non-planar ligands often give better results with ``slab``. ``auto`` will choose the shape automatically based on the ligand geometry.

    Tip: This option can severely affect the quality of generated complexes and how many make it through the post-assembly filter. For serious applications we recommend to set :confval:`max_num_complexes` to ``100``, try all three options and check how many complexes fail the post-assembly filter for each option (this info is returned at the end of the assembly if :confval:`verbosity` >= ``2``). Whichever option has the least complexes failing the post-assembly filter usually gives the highest quality complexes.

.. confval:: geometry_modifier_filepath

    :options: `empty` OR `filepath`
    :required: ``false``
    :default: ``empty``

    Path to the geometry modifier file. If left empty, no geometry modification is performed.

    The geometry modifier file allows very advanced and fine-grained control over the geometry of the generated complexes. Usually it is not needed, since a forcefield optimization will often be a better option. However, there might be cases where it is desired to move atoms in an assembled ligand from one position to another position for all complexes with this ligand. This can be achieved with the geometry modifier file as shown in the Pd/Ni cross coupling example.

    For moving an atom to another position you need to supply the chemical symbol and the coordinates of the original atom and the coordinates the atom at it's new coordinates. The geometry modifier file is an .xyz file with two sets of atoms: The first set is all atoms that should be moved, the second set is the new positions of these atoms. Both sets of atoms are provided as "molecule" in the .xyz format and concatenated. The order and the chemical elements of both sets of atoms have to match up. In the assembly, for each generated complex, the atoms with coordinates in the first set are moved to the coordinates in the second set.

.. confval:: complex_name_appendix

    :options: `empty` or `string`
    :required: ``false``
    :default: `empty`

    Appends a custom string to the randomly generated name of each assembled complex. For example, if the appendix is set to ``_charge1``, a generated complex will be named 'ZUMUVAMI_charge1' if otherwise it would have been named 'ZUMUVAMI'.

Global Options
-----------------------------------

Global options are all optional and specify settings that apply to all batches.

.. confval:: output_directory

    :options: `dirpath`
    :required: ``false``
    :default: ``DART``

    Path to directory in which the output will be saved.

.. confval:: ffmovie

    :options: ``true``, ``false``
    :required: ``false``
    :default: ``false``

    Whether to output a movie (i.e. a concatenated .xyz file displaying multiple frames) of the forcefield optimization process. Useful for visualization e.g. with ``ase gui FILE.xyz``.

.. confval:: concatenate_xyz

    :options: ``true``, ``false``
    :required: ``false``
    :default: ``true``

    Whether to save concatenated xyz files with all passed/failed complexes respectively. Useful for quick visualization and browsing of the generated complexes e.g. with ``ase gui FILE.xyz``.

.. confval:: verbosity

    :options: ``0``, ``1``, ``2``, ``3``
    :required: ``false``
    :default: ``2``

    How much output to print (except the progress bars, which are always printed). ``0`` means only errors, ``1`` means also warnings, ``2`` means also normal info, ``3`` means also debug info.

.. confval:: same_isomer_names

        :options: ``true``, ``false``
        :required: ``false``
        :default: ``true``

        Whether to give isomers of the same complex the same name and then just number them. If set to ``false``, each isomer will get a completely unique name.

.. confval:: complex_name_length

    :options: `integer > 0`
    :required: ``false``
    :default: ``8``

    Length of the randomly generated name for each generated complex (e.g. 'ZUMUVAMI').

Table: overview of all available assembler options
---------------------------------------------------

.. csv-table::
   :file: assembly_options_overview.csv
   :widths: 20, 10, 20, 50
   :header-rows: 1



.. include:: DART_assembly_output_explanation.rst