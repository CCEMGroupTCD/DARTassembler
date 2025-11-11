"""
Assembler Input
-----------------------------------

The DART Assembler Module generates 3D structures of novel transition metal complexes from a database of ligands, which can either be the full :ref:`MetaLig database <metalig>` or a subset from a user-defined chemical space. While this page focuses on the input options for the assembler, you can read more about how the DART Assembler Module works :ref:`here <how_assembler_works>`.

The assembler module is run from the command line by providing a single configuration file:

.. code-block:: bash

    DARTassembler assembler --input assembler.yml

**Copy-Paste Template:**

.. literalinclude:: ../../DARTassembler/data/default/assembler.yml
   :language: yaml
   :linenos:

Users can download this template into their current directory as ``assembler.yml`` by running:

.. code-block:: bash

    DARTassembler configs --outdir .

.. _ligand_archetypes_and_target_vectors:

Ligand Archetypes and Target Vectors
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _ligand_archetypes:

**Ligand Archetypes :**

Each ligand in the MetaLig database is classified into one of 22 ``archetypes`` based on the orientation of its donor atoms around the metal center. For example, bidentate ligands are either ``2-cis`` or ``2-trans``, and tridentate ligands either ``3-facial``, ``3-meridional`` or ``3-trigonal``. Ligands of denticity 2-6 can have one of multiple archetypes defined, while ligands of denticity 1 and 7-10 have simply all ligands collected in a single archetype. The full list of defined archetypes is shown in Figure 1.

.. figure:: /_static/archetypes.png
   :width: 100%
   :align: center

   Figure 1: Graphical overview of all 22 ligand archetypes defined in DART.

.. _target_vectors:

**Target Vectors :**

The ``target_vectors`` option specifies the orientation of each ligand's donor atoms relative to the metal center during complex assembly. For example, for a ``1-mono`` ligand one has to specify a single vector (e.g. ``['z']``), meaning the metal-donor bond will point along the +Z axis, placing the monodentate on top of the metal center. For a ``2-cis`` ligand, two orthogonal vectors must be provided (e.g. ``['x', 'y']``), meaning the two donor atoms will point along the +X and +Y axes, placing the bidentate ligand in a cis arrangement in the x-y plane. For a ``2-trans`` ligand, two opposite vectors must be provided (e.g. ``['x', '-x']``), meaning the ligand will coordinate to the metal center from opposite sides along the X axis. For a ``3-trigonal`` ligand, three vectors separated by 120° must be provided (e.g. ``['xy(0)', 'xy(120)', 'xy(240)']``), meaning the three donor atoms will point in a trigonal planar arrangement around the metal center in the x-y plane.

Target vectors can be specified in three formats:

- Symbolic axis indicators: ``'+x'``, ``'-x'``, ``'+y'``, ``'-y'``, ``'+z'``, ``'-z'`` for the Cartesian axes. For example, ``'+x'`` means the vector ``[1, 0, 0]`` and ``'-z'`` means the vector ``[0, 0, -1]``. The specification ``'x'`` is equivalent to ``'+x'``.
- Angled vectors in the x-y plane from the y-axis: ``'xy(angle_in_degrees)'``. For example, ``'xy(120)'`` means the vector ``[sin(120°), cos(120°), 0]``.
- Explicit Cartesian triplets: e.g. ``[0, 0, 1]``, ``[-1, 0, 0]``, etc.


Table 1 provides a full list of all defined ligand archetypes and an example set of target vectors for each archetype.

.. csv-table:: Table 1: Defined ligand archetypes and example target vectors
   :file: ../../../DARTassembler/data/docs/archetypes.csv
   :header-rows: 1
   :widths: 20, 80
   :align: center

In general, you can specify any set of target vectors you like, as long as they match the ligand archetype you want to use for the binding site.

For example, if a user wants to assemble an octahedral structure with one monodentate on top, one mer-tridentate in the equatorial plane and one cis-bidentate occupying the remaining two faces, the target vectors would be specified as:

.. code-block:: yaml

    # Example: Octahedral Ru complex with 1-mono, 2-cis and 3-mer ligands
    metal_centers: Ru
    ligand_archetypes:
      - '1-mono'
      - '2-cis'
      - '3-meridional'
    target_vectors:
      - ['+z']                # monodentate pointing along +Z
      - ['+x', '-z']          # cis-bidentate pointing along +X and -Z
      - ['+y', '-x', '-y']    # mer-tridentate pointing along +Y, -X and -Y

.. figure:: /_static/Ru_octahedral.png
   :width: 50%
   :align: center

   Figure 1: Octahedral complex geometry defined by the ``target_vectors`` above.

Of course, this is not the only choice of target vectors. Many other combinations are possible and will result in the same complex, just rotated.

The provision of target vectors allows the user fine-grained control over the geometry of the assembled complexes. Combined with the full list of 22 different ligand archetypes, the DART assembler can create an extensive variety of complex geometries.

Examples
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Let's try another example: a square-planar complex with two monodentates and one trans-bidentate. The ligand archetypes and target vectors could be specified as:

.. code-block:: yaml

    # Example: Square-planar complex with two 1-mono and one 2-trans ligands
    ligand_archetypes:
        - '1-mono'
        - '1-mono'
        - '2-trans'
    target_vectors:
        - ['+x']            # 1. monodentate along +X
        - ['-x']            # 2. monodentate along -X
        - ['+y', '-y']      # trans-bidentate along +Y and -Y

One more example for a trigonal-bipyramidal complex with one trigonal ligand and two monodentates:

.. code-block:: yaml

    # Example: Trigonal-bipyramidal complex with one 3-trigonal and two 1-mono ligands
    ligand_archetypes:
        - '1-mono'
        - '1-mono'
        - '3-trigonal'
    target_vectors:
        - ['+z']                            # 1. monodentate along +Z
        - ['-z']                            # 2. monodentate along -Z
        - ['xy(0)', 'xy(120)', 'xy(240)']   # trigonal ligand in x-y plane

The same trigonal-bipyramidal complex, but this time with five monodentates:

.. code-block:: yaml

    # Example: Trigonal-bipyramidal complex with five 1-mono ligands
    ligand_archetypes:
        - '1-mono'
        - '1-mono'
        - '1-mono'
        - '1-mono'
        - '1-mono'
    target_vectors:
        - ['+z']            # 1. monodentate along +Z
        - ['-z']            # 2. monodentate along -Z
        - ['xy(0)']         # 3. monodentate in x-y plane at 0°
        - ['xy(120)']       # 4. monodentate in x-y plane at 120°
        - ['xy(240)']       # 5. monodentate in x-y plane at 240°

Finally, a tetrahedral complex with four monodentates. For the tetrahedral geometry, the target vectors are a little more complicated, but they are actually given in Table 1 under the ``4-tetrahedral`` ligand archetype. Since DART has no abbreviation for tetrahedral vectors, we will provide the explicit Cartesian triplets, as one can always do:

.. code-block:: yaml

    # Example: Tetrahedral complex with 1 cis-bidentate and two monodentates
    ligand_archetypes:
        - '1-mono'
        - '1-mono'
        - '2-cis'
    target_vectors:
        - [ [1.0,1.0,-1.0] ]                    # 1. monodentate
        - [ [-1.0,1.0,1.0] ]                    # 2. monodentate
        - [ [1.0,-1.0,1.0], [-1.0,-1.0,-1.0] ]  # cis-bidentate

Haptically coordinating ligands
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :ref:`MetaLig database <metalig>` also contains many haptically coordinating ligands. DART can assemble these haptically coordinating ligands just as well as non-haptic ones. In haptic ligands, each group of haptic donor atoms is treated as a single pseudo donor atom located at the centroid of the haptic group. For example, a Cp* ligand has a single group of 5 haptic donor atoms. Thus, it is treated as a ``1-mono`` ligand and can be assembled by providing a single target vector. For example, to assemble a Cp* ligand such that it coordinates from above the metal center, one would provide the target vector ``['+z']`` (see the :ref:`advanced example <advanced_example>`). One can query haptically coordinating ligands in the :ref:`LigandFilters <ligandfilters>` by targeting the properties ``n_eff_denticities``, ``n_denticies``, ``n_haptic_groups``, and ``n_haptic_atoms``.


.. _assembly_output:

Assembler Output
--------------------

The output of the DART assembler will be saved in a specific folder. This folder is determined by the :confval:`output_directory` you set in your assembly input file. Within this folder, each generated metal complex has a unique name such as ``'ZUMUVAMI'``. This name is randomly generated by DART but you can control its length via the :confval:`complex_name_length` option in the input file. If you want to append a custom string to each complex name for labeling purposes, you can do so via the :confval:`complex_name_appendix` option.

The assembler module creates not only the xyz files for each complex but also various other files that could be of interest. Below, you'll find an overview of all files and folders generated by the DART assembler.

Folder Structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Here's what the output folder will look like::

    output_directory/
    ├── isomers.csv                     (Summary Table)
    ├── input/                          (Copy of Input .yml)
    ├── log.txt                         (Log File)
    └── batches/                        (Batch Folders)
        ├── batch_1/
        │   ├── concat_passed_complexes.xyz
        │   ├── concat_failed_complexes.xyz
        │   ├── concat_all_complexes.xyz
        │   └── complexes/              (Individual Complex Folders)
        │       └── NAME/
        │           ├── NAME.json       (Complex Data File)
        │           ├── NAME1.xyz       (Isomer 1 Structure)
        │           ├── NAME2.xyz       (Isomer 2 Structure)
        │           └── ...
        └── batch_2/
            └── ...

Output Files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

General Output Files:
    These files provide a broad overview of the assembly process:

    - ``batches/``: This is a folder that contains all the batches of assembled complexes.
    - ``isomers.csv``: This is a summary table listing all isomers of all complexes generated across all batches. It includes details such as complex names, ligand combinations, archetypes, stoichiometries, and assembly status (successful or failed).
    - ``input/``: This folder contains a copy of the input .yml file used to run the assembler for record-keeping.
    - ``log.txt``: This is the main log file that records all messages, warnings, and errors during the assembly process.

Batch-Specific Files:
    Inside the ``batches/`` folder, you'll find separate folders for each batch. These folders may contain:

    - ``complexes/``: This is a folder that contains all the complexes for that batch.
    - ``concat_passed_isomers.xyz``: Concatenated .xyz file of all successfully assembled complexes.
    - ``concat_failed_isomers.xyz``: Concatenated .xyz file of all complexes that failed assembly.
    - ``concat_all_isomers.xyz``: Concatenated .xyz file of all complexes, both successful and failed.

    If a file is missing, that means no complexes fall into that category for the batch (e.g., no failed complexes). Concatenated .xyz files can easily be browsed using the ``ase gui`` command from the ASE package:

    .. code-block:: bash

        ase gui concat_passed_isomers.xyz

Complex-Specific Files:
    Within each batch, each complex has its own folder under ``complexes/NAME/``. These folders contain:

    - ``NAME.json``: A machine-readable data file with detailed information about the complex and its isomers.
    - ``NAME1.xyz``: 3D structure of the 1. isomer of the complex.
    - ``NAME2.xyz``: 3D structure of the 2. isomer of the complex.
    - ...

.. _assembler_parameters:

Assembler Options
-----------------------------------
The DART assembler module is usually run via a .yml configuration file. This .yml file will then be passed to the following class. Please refer to the docstrings of this class and its ``run_batch()`` method for a detailed description of all available options. You will see that the options match perfectly with the options in the .yml configuration template above.

"""
# Standard library imports
from typing import Union, Dict, Any, List, Tuple, Optional
import ase.io
from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm
import datetime
import logging
import random
import sys

# DART specific imports
from DARTassembler.src.assembler.isomer import AssembledIsomer, AssembledComplex
from DARTassembler.src.assembler.output import AssemblerOutput, BatchAssemblerOutput
from DARTassembler.src.metalig.utils_molecule import get_standardized_stoichiometry_from_atoms_list
from DARTassembler.src.constants.paths import default_assembler_yml_path
from DARTassembler.src.assembler.ligands import LigandChoice
from DARTassembler.src.modules.modules import BaseModule
from DARTassembler.src.metalig.db import LigandDB
from DARTassembler.src.misc.io import read_yaml, get_correct_ligand_db_path_from_input, save_json
from DARTassembler.src.metalig.archetype import valid_archetypes

# Data processing imports
from pathlib import Path
import pandas as pd
import numpy as np

# Warnings
import warnings
warnings.simplefilter("always")

class Assembler(BaseModule):
    """
    Assemble isomers of transition-metal complexes from ligand databases such as the MetaLig.
    """

    def __init__(self,
                 output_directory: Union[str, Path] = 'DARTassembler',
                 verbosity: int = 2,
                 complex_name_length: int = 8,
                 n_max_ligands: Optional[int] = None,
                 ) -> None:
        """
        Initialize the DART Assembler module. The options set here applied to all batches.

         .. tip:: All the parameters below are available as well via the assembler .yml file as global options (i.e. without indentation).

        :param str | Path output_directory: Directory to save the DART assembler output files.
        :param int verbosity: Logging verbosity (0=errors only, 1=warnings, 2=info, 3=debug).
        :param int complex_name_length: Length of random complex names such as 'ZUMUVAMI'. Increases automatically if otherwise a name clash would occur.
        :param int | None n_max_ligands: Maximum number of ligands to load from each :confval:`ligand_db_file`. If ``None`` or left unspecified, all ligands are loaded.
        :return: None
        :rtype: None
        """
        super().__init__()
        self.output_directory = Path(output_directory).resolve()
        self.verbosity = verbosity
        self.complex_name_length = complex_name_length
        self.n_max_ligands = n_max_ligands

        self._loaded_ligand_databases = {}  # to avoid reloading the same ligand database multiple times

        # Keep track of the input arguments
        self.init_args = {**locals()}
        self.init_args.pop('self')
        self.init_args.pop('__class__')

        # Set up the output directories
        self.gbl_outcontrol = AssemblerOutput(outdir=self.output_directory)

        # Set up logging
        verbosity2logging = {0: logging.ERROR, 1: logging.WARNING, 2: logging.INFO, 3: logging.DEBUG}
        stream_handler = logging.StreamHandler(stream=sys.stdout)
        stream_handler.setLevel(verbosity2logging[self.verbosity])
        # Print to stdout
        logging.basicConfig(level=verbosity2logging[self.verbosity], format='%(message)s', handlers=[logging.FileHandler(self.gbl_outcontrol.log_path, mode='w'), stream_handler], force=True)

    def run_batch(self,
                  name: str,
                  batch_idx: int,
                  target_vectors: list[list[list[float]]],
                  metal_centers: Union[str, tuple[str, tuple[float, float, float]]],
                  n_max_complexes: Union[int, str],
                  ligand_db_files: Union[list[str], str] = 'metalig',
                  ligand_archetypes: list[str] = None,
                  ligand_origins: list[tuple[float, float, float]] = None,
                  total_ligand_charges: int = None,
                  monoaxial_optimization: bool = True,
                  permutable_ligands: list[int] = None,
                  force_all_isomers: bool = False,
                  duplicate_tolerance: float = 0.5,
                  clashing_tolerance: float = -0.3,
                  clashing_metal: bool = False,
                  complex_name_suffix: str = '',
                  random_seed: Optional[int] = None,
                  background_file: str = None,
                  background_translation: Optional[Tuple[float, float, float]] = None,
                  ) -> None:
        """
        Run the DART assembler for a single batch.

        :param int batch_idx: Index of the batch in the list of batches. This is the only option that is not available in the input .yml file because it is set automatically when running multiple batches.

        .. tip::

            All the parameters below are available as well via the assembler .yml file as batch options (i.e. indented in the ``batches:`` list).

        :param str name: Name of the batch. All batch names must be unique.
        :param target_vectors: List of target vectors for each ligand. Each entry is a list of donor vectors of length 3.
        :type target_vectors: list[list[list[float, float, float]]]  # shape: (n_ligands, n_donors_per_ligand, 3)
        :param metal_centers: Metal center specification. Either a single element symbol string (e.g. 'Ru') or a nested structure per ligand. Examples:

            - Single Ru atom at origin (0,0,0): ``Ru``
            - Single Ru atom at custom position: ``['Ru', [1.5, 0.0, -1.0] ]``
            - Ru and Cu atoms 4 Å apart: ``[ ['Ru', [-2.0, 0.0, 0.0] ], ['Cu', [2.0, 0.0, 0.0] ] ]``

        :type metal_centers: str | list[str, list[float]]] | list[list[str, list[float]]]
        :param int | str n_max_complexes: Maximum number of complexes (not isomers) to assemble or ``'all'`` to exhaust all possible combinations of ligands (respecting ``total_ligand_charges``).

         .. warning:: Due to to the combinatorial explosion of possible complexes, setting this to ``'all'`` can lead to very very many complexes and long runtimes if used with multiple ligand databases with each many ligands. Use with caution.
        :param list[str] | str ligand_db_files: If left unspecified or set to ``'metalig'``, the entire MetaLig database is used. If a single path to a ligand .jsonlines file, that database is used for all ligand sites. If a list of strings, the list must have the same length as ``target_vectors`` and each entry is the ligand database that will be sampled to populate the respective ligand site. One exception: the keyword ``'same_as_previous'`` can be used instead of a filepath to indicate that the ligand for that site will always be populated with the same ligand as the one chosen from the previous ligand database.
        :param list[str] | None ligand_archetypes: If specified, filters the ligand database at this index to only contain ligands of the specified archetype. The list must have the same length as ``target_vectors``. If left unspecified, no filtering is applied and all ligands in the database with matching ``n_eff_denticity`` are considered.
        :param ligand_origins: If specified, applies a shift to the position of the ligands. If not specified, the metal center position(s) are used for ligands coordinating to one metal and the midpoint between two (or more) metal centers for bridging ligands. Must have the same length as ``target_vectors``.
        :type ligand_origins: list[ list[float, float, float] ] | None
        :param int | None total_ligand_charges: If specified, only ligand combinations with this sum of ligand charges are considered. If left unspecified or None, no charge filtering is applied. This allows to control the overall charge and metal oxidation state of the assembled complex. For example, to generate neutral Pd(II) complexes the ``total_ligand_charges`` should be set to -2.
        :param bool monoaxial_optimization: Whether to optimize the orientation of monoaxial ligands (``1-mono`` and ``2-trans`` archetypes) around their single binding axis to minimize clashes and maximize distance to other ligands. In general recommended but significantly increases runtime. Can be turned off without problems for much faster assembly and when the orientation of monoaxial ligands relative to the other ligands is not critical, especially for small monodentate ligands.
        :param list[int] | None permutable_ligands: Groups of ligands that should be permuted when generating isomers. Must have same length as ``target_vectors``. Only ligands with the same archetype can be permuted. If None or left unspecified, no ligands are permuted. For example, if you generate octahedral complexes, the following will permute the two monodentate ligands but not permute the two bidentate ligands (think of each integer as a "color" assigned to each ligand; ligands with the same color are permuted among each other):

            .. code-block:: yaml

                    target_vectors: [ [ '+z' ], [ '+x' ], [ '-x', '-y' ], [ '+y', '-z' ] ]
                    ligand_archetypes: ['1-mono', '1-mono', '2-cis', '2-cis']
                    permutable_ligands: [1, 1, 2, 3]

        :param bool force_all_isomers: If False (default), DART does not generate the following types of isomers because they are considered symmetrically equivalent:

            - For ring-like ligands with archetypes ``3-trigonal``, ``4-tetragonal``, ``5-pentagonal``, and ``6-hexagonal``, DART removes isomers which correspond to just rotating the ligand around it's centering axis. E.g. for a trigonal ligand, rotating the ligand by 120° or 240° around the metal center is not considered a new isomer. Ligands with these archetypes will therefore return only two isomers, facing up and down the centering axis (mirror images of each other).

            - Ring-like ligands with archetypes ``4-trigonal-pyramidal``, ``5-square-pyramidal``, and ``6-pentagonal-pyramidal`` can also be rotated around the z-axis. However, this time there is an additional atom on top of the ring, breaking the z mirror symmetry, therefore returning only a single isomer.

            - For ligands with archetypes ``4-tetrahedral``, ``6-octahedral``, ``6-trigonal-prismatic``, ``7-septa``, ``8-octa``, ``9-nona``, and ``10-deca``, the ligand is so bulky that it is assumed there are no other ligands around the complex and only one isomer is generated.

            If ``force_all_isomers`` is set to True, DART will generate these types of isomers as well. However, exact 3D duplicates of isomers will still be filtered out in the post-assembly checks if ``duplicate_tolerance`` is set.

        :param float | None duplicate_tolerance: Tolerance used to identify duplicate isomers in the post-assembly checks. Increase this value to make more isomers pass the duplicate filter. Decrease (up to 0.0) to make the filter more strict. If None, no duplicate filtering is applied.
        :param float | None clashing_tolerance: Tolerance for detecting clashing ligands in the post-assembly checks. The tolerance will be added to the sum of the van der Waals radii of two atoms to determine whether they clash. Increase this value to make the clash filter more strict. If None, no clash filtering is applied.
        :param bool clashing_metal: If True, include ligand-metal and metal-metal pairs in clash checks.
        :param str complex_name_suffix: Optional suffix to append to each generated complex name, e.g. ``'ZUMUVAMI_OOH'`` if suffix is set to ``'_OOH'``. Useful for tracking complexes in downstream analysis.
        :param int | None random_seed: Random seed for reproducibility. Defaults to the batch index if None.
        :param str | None background_file: EXPERIMENTAL OPTION: Optional path to an .xyz file to use as a "background" structure for each isomer. This background structure will be combined with each generated isomer and saved as an additional structure named ``<isomer_name>_combined.xyz`` in each complex directory. The background structure can contain e.g. solvent molecules or a surface. The background structure is an experimental option that is not considered during assembly, clash checks, or duplicate checks; and it is neither present in the concatenated .xyz files nor in the summary .csv file or the complex .json files. It is only saved as an additional single .xyz file per isomer for user convenience. If None, no background structure is used.
        :param list[float, float, float] background_translation: EXPERIMENTAL OPTION: Translation applied to the background structure. If e.g. ``[0, 0, -1]``, the background structure is shifted by -1 Å in z-direction before combining with the isomer. If None, no translation is applied.
        :return: None
        :rtype: None
        """
        # Set random seed for reproducibility. Do this batch-wise so every batch is reproducible independently.
        if pd.isna(random_seed):
            random_seed = batch_idx
        random.seed(random_seed)
        np.random.seed(random_seed)

        # Handle defaults
        if isinstance(ligand_db_files, str):    # Expand a single path to a list of paths for each ligand.
            ligand_db_files = [ligand_db_files for _ in target_vectors]
        # Validate input
        if ligand_archetypes:
            if len(ligand_archetypes) != len(target_vectors):
                raise ValueError(f"When providing 'ligand_archetypes', the length of the list must match the number of ligands ({len(target_vectors)}). Provided length: {len(ligand_archetypes)}.")
            invalid_archetypes = [arch for arch in ligand_archetypes if arch not in valid_archetypes]
            if invalid_archetypes:
                raise ValueError(f"The following provided ligand archetypes are invalid: {invalid_archetypes}. Please choose from the following valid archetypes: {valid_archetypes}.")
            for arch in ligand_archetypes:
                denticity = int(arch.split('-')[0])
                n_target_vectors = len(target_vectors[ligand_archetypes.index(arch)])
                if denticity != n_target_vectors:
                    raise ValueError(f"The provided ligand archetype '{arch}' has an effective denticity of {denticity}, which does not match the number of target vectors ({n_target_vectors}) specified for that ligand position. Please ensure that the archetype and target vectors are consistent.")

        self.batch_name = name
        self.batch_idx = batch_idx
        self.ligand_db_files = ligand_db_files   # do not resolve path here to keep keywords like 'metalig', 'same_as_previous'
        self.n_max_complexes = n_max_complexes
        self.random_seed = random_seed
        self.total_ligand_charges = total_ligand_charges
        self.target_vectors = target_vectors
        self.ligand_origins = ligand_origins
        self.metal_centers = metal_centers
        self.complex_name_suffix = complex_name_suffix
        self.permutable_ligands = permutable_ligands
        self.clashing_tolerance = clashing_tolerance
        self.clashing_metal = clashing_metal
        self.duplicate_tolerance = duplicate_tolerance
        self.monoaxial_optimization = monoaxial_optimization
        self.ligand_archetypes = ligand_archetypes
        self.force_all_isomers = force_all_isomers
        self.batch_output_path = Path(self.gbl_outcontrol.batch_dir, self.batch_name)
        self.batch_outcontrol = BatchAssemblerOutput(self.batch_output_path)
        self.background_file = background_file
        self.background_translation = background_translation

        # Redirect tqdm to the logging module so that messages appear properly on two different lines
        with (logging_redirect_tqdm()):

            # Load the ligand databases and cache them for later use
            self.ligand_dbs = self._get_ligand_databases()

            # Set up an iterator for the ligand combinations
            ligand_choice = LigandChoice(
                                    ligand_dbs=self.ligand_dbs,
                                    total_ligand_charges=self.total_ligand_charges,
                                    n_max_complexes=self.n_max_complexes,
                                    )
            ligand_combinations = ligand_choice.choose_ligands()

            # Set progress bar with or without final number of assembled complexes
            total = self.n_max_complexes if self.n_max_complexes != 'all' else None
            progressbar = tqdm(desc='Assembling complexes', unit=' complexes', file=sys.stdout, total=total, disable=self.verbosity < 2)

            batch_sum_assembled_complexes = 0  # Number of assembled complexes in this batch
            while ligand_choice.if_make_more_complexes(batch_sum_assembled_complexes):

                # Choose ligands for complex
                try:
                    ligands = next(ligand_combinations)
                except StopIteration:
                    break # If all ligand combinations are exhausted, stop the batch

                complex = AssembledComplex(
                    ligands=ligands,
                    target_vectors=self.target_vectors,
                    ligand_origins=self.ligand_origins,
                    metal_centers=self.metal_centers,
                )
                complex.generate_isomers(
                                            clashing_tolerance= self.clashing_tolerance,
                                            clashing_metal= self.clashing_metal,
                                            duplicate_tolerance= self.duplicate_tolerance,
                                            permutable_ligands = self.permutable_ligands,
                                            monoaxial_optimization = self.monoaxial_optimization,
                                            force_all_isomers=self.force_all_isomers,
                                            complex_name_length=self.complex_name_length,
                                            complex_name_suffix=self.complex_name_suffix,
                                            avoid_names=self.all_tried_complex_names,  # Avoid names of already tried complexes
                )
                # Add the complex name to the set of all tried complex names to avoid duplicates in the next iteration
                self.all_tried_complex_names.add(complex.complex_name)

                self._save_assembled_isomers(complex=complex)

                # Update counters if at least one isomer was successfully assembled for this complex
                if complex.success:
                    batch_sum_assembled_complexes += 1
                    progressbar.update(1)

            progressbar.close()

        return

    def run(self, batches: list[dict]) -> None:
        """
        Execute assembly for a sequence of user-defined batches.

        Each batch dictionary describes inputs for a single assembly run. This method
        validates batch names, logs run metadata, saves the input settings and iterates
        over batches calling the internal batch runner.

        :param batches: List of batch specification dictionaries as expected by _run_batch.
        :type batches: list[dict]
        :return: None
        :rtype: None
        """
        self.batches = batches
        for idx, batch in enumerate(self.batches):
            batch['name'] = batch.get('name', f'batch_{idx}')  # Use batch name or generate a default one
        self.n_batches = len(self.batches)
        self.df_info = []
        self.all_tried_complex_names = set()    # Keep track of all tried complex names to avoid duplicates
        self.successfully_assembled_isomer_names = []

        self.batches_args = {**locals()}
        self.batches_args.pop('self')
        self.input_args = {**self.init_args, **self.batches_args}
        self._log_global_info()

        self._check_batch_settings(batches)

        # Save yml file with input arguments to output directory
        self.gbl_outcontrol.save_settings(self.input_args)

        start = datetime.datetime.now()
        for idx, batch_settings in enumerate(self.batches):
            self._log_batch_title_and_settings(batch_settings=batch_settings)
            self.run_batch(batch_idx=idx, **batch_settings)  # run the batch assembly

        self.runtime = datetime.datetime.now() - start  # keep track of the runtime to display later
        self._make_and_save_output_csv()
        self._log_summary()
        self._final_checks()

        return

    @staticmethod
    def _check_batch_settings(batches: list) -> None:
        """
        Validate that batch names are unique.

        :param batches: List of batch specification dictionaries; each must contain a 'name' key.
        :type batches: list[dict]
        :raises ValueError: If any batch names are duplicated.
        :return: None
        :rtype: None
        """
        # Check all batch names are unique
        batch_names = [batch['name'] for batch in batches]
        if not len(batch_names) == len(set(batch_names)):
            raise ValueError(f"DART batch names must be unique. The following batch names are not unique: {batch_names}")

        return

    def _get_ligand_databases(self) -> list[LigandDB]:
        """
        Load and prepare ligand databases for the current batch.

        Each entry in self.ligand_db_files is resolved to a LigandDB. Databases are cached
        to avoid repeated loading. Ligands are filtered by effective denticity and,
        optionally, by archetype.

        :return: List of LigandDB objects (one per ligand position).
        :rtype: list[ LigandDB ]
        :raises ValueError: If a required database contains no ligands matching the requested denticity/archetype.
        """
        ligand_databases = []
        for idx, (target_vectors, ligand_db_filepath) in enumerate(zip(self.target_vectors, self.ligand_db_files)):
            if ligand_db_filepath == 'same_as_previous':
                if not idx > 0:
                    raise ValueError("The first ligand database cannot be 'same_as_previous'. Please provide a valid ligand database file path.")
                ligand_databases.append('same_as_previous')
                continue

            ligand_db_filepath = get_correct_ligand_db_path_from_input(ligand_db_filepath)
            if not ligand_db_filepath in self._loaded_ligand_databases:
                self._loaded_ligand_databases[ligand_db_filepath] = LigandDB.from_json(ligand_db_filepath, n_max=self.n_max_ligands, show_progress=self.verbosity > 1)

            # Reduce to the ligands which have the correct n_eff_denticity for the specified target vectors
            n_target_vectors = len(target_vectors)
            database = {name: ligand for name, ligand in self._loaded_ligand_databases[ligand_db_filepath].db.items() if ligand.n_eff_denticities == n_target_vectors}

            # If required, reduce the database to the ligands with the correct archetype
            if self.ligand_archetypes is not None:
                archetype = self.ligand_archetypes[idx]
                database = {name: ligand for name, ligand in database.items() if ligand.archetype == archetype}

            if not database:
                with_archetypes = f' and ligand archetype `{archetype}`' if self.ligand_archetypes is not None else ''
                raise ValueError(f"The provided ligand database contains no ligands with `n_eff_denticities={n_target_vectors}`{with_archetypes}. Please check your input ligand database `{Path(ligand_db_filepath).resolve()}`.")

            database = LigandDB(database)  # Convert to LigandDB object
            ligand_databases.append(database)

        return ligand_databases

    def _save_assembled_isomers(self, complex):
        """
        Persist assembled complex isomers and per-isomer XYZ files to disk.

        For successful complexes a JSON with all isomer metadata is written. Every isomer's
        XYZ string is saved to the complex folder and appended to batch-level concatenated XYZs.
        If a background file is provided, an additional combined XYZ is written.

        :param complex: AssembledComplex instance with generated isomers and metadata.
        :type complex: AssembledComplex
        :return: None
        :rtype: None
        """
        # Save complex json file with all isomer data
        if complex.success:
            complex_dir = Path(self.batch_outcontrol.complex_dir, complex.complex_name)
            complex_json_filepath = complex_dir / f'{complex.complex_name}.json'
            data = complex.to_dict()
            save_json(db=data, path=complex_json_filepath, mkdir=True, indent=4)

        for isomer in complex.isomers:
            success = True if isomer.warning == '' else False
            isomer_name = isomer.isomer_name

            # Add a comment to the xyz file with complex and ligand names for easier identification
            xyz_comment = f'isomer_name: {isomer_name}, warning: {isomer.warning}, ligand_unique_names: ({", ".join(isomer.ligand_info["unique_names"])})'
            xyz_string = isomer.get_xyz_string(comment=xyz_comment)

            # Save xyz of isomer to complex directory
            if success:
                isomer_xyz_filepath = complex_dir / f'{isomer_name}.xyz'
                with open(str(isomer_xyz_filepath), 'w') as xyz_file:
                    xyz_file.write(xyz_string)
                self.successfully_assembled_isomer_names.append(isomer_name)
                if self.background_file is not None:
                    # Combine the DART generated xyz with the extra structure xyz into one new xyz file and save it in the same location

                    extra_structure = ase.io.read(self.background_file, format="xyz")

                    # Ensure we always use a NumPy array for the translation vector
                    translation = np.array(self.background_translation or [0.0, 0.0, 0.0], dtype=float)
                    extra_structure.positions += translation

                    # Combine structures (order doesn’t matter unless you care about atom order)
                    combined_structure = extra_structure + isomer.atoms

                    # combined_structure.set_positions(np.round(combined_structure.get_positions(), decimals=6))

                    # Write combined XYZ
                    combined_xyz_filepath = complex_dir / f"{isomer_name}_combined.xyz"
                    ase.io.write(str(combined_xyz_filepath), combined_structure, format="xyz")


            # Save to concatenated xyz files of this batch
            self.batch_outcontrol.save_xyz(xyz_string, success=success, append=True)    # passed/failed xyz files are created automatically
            self.batch_outcontrol.save_file(xyz_string, self.batch_outcontrol.all_xyz_path, append=True)        # save to all_xyz file

            # Save data for csv file.
            isomer_idx = (len(self.successfully_assembled_isomer_names) - 1) if success else None
            self._add_batch_info(isomer=isomer, success=success, isomer_idx=isomer_idx, complex_name=complex.complex_name)

        return

    def _add_batch_info(self, isomer: AssembledIsomer, success, isomer_idx: int, complex_name: str) -> None:
        """
        Append a row of per-isomer metadata to the internal batch info list.

        The resulting entries are later assembled into a pandas DataFrame and saved as CSV.

        :param AssembledIsomer isomer: The isomer instance to extract metadata from.
        :param bool success: True if the isomer passed post-filters.
        :param int isomer_idx: Index within successful isomers list or None for failures.
        :param str complex_name: Name of the parent complex.
        :return: None
        :rtype: None
        """
        elements = isomer.get_metal_symbols()
        metal_stoi = get_standardized_stoichiometry_from_atoms_list(elements)
        data = {
            'success': success,
            'isomer_idx': isomer_idx,
            'isomer_name': isomer.isomer_name,
            'complex_name': complex_name,
            'stoichiometry': isomer.stoichiometry,
            'graph_hash': isomer.graph_hash,
            'warning': isomer.warning,
            'ligand_unique_names': isomer.ligand_info['unique_names'],
            'ligand_archetypes': isomer.ligand_info['archetypes'],
            'ligand_stoichiometries': isomer.ligand_info['stoichiometries'],
            'ligand_charges': isomer.ligand_info['charges'],
            'ligand_donors': isomer.ligand_info['donors'],
            'batch_idx': self.batch_idx,
            'batch_name': self.batch_name,
            'metal_centers': metal_stoi,
            'total_ligand_charges': self.total_ligand_charges,
            'random_seed': self.random_seed,
        }
        self.df_info.append(data)

        return

    @classmethod
    def run_from_yaml(cls, input: Union[str, Path, None], n_max_ligands:Optional[int]=None) -> 'Assembler':
        """
        Instantiate and run an Assembler using a YAML configuration file.

        If input is None, the project's default assembler YAML template is used.
        The YAML must contain top-level options accepted by Assembler.__init__ and a 'batches' list.

        :param input: Path to YAML configuration file or None to use the default template.
        :type input: Union[str, Path, None]
        :param n_max_ligands: Optional override for the maximum number of ligands to load from each ligand database. Takes precedence over the value in the YAML file if specified.
        :type n_max_ligands: int | None
        :return: Assembler instance after executing the specified batches.
        :rtype: Assembler
        """
        if input is None:
            input = default_assembler_yml_path

        options = read_yaml(input)
        if n_max_ligands is not None:
            options['n_max_ligands'] = n_max_ligands
        batches = options.pop('batches')

        assembler = Assembler(**options)
        assembler.run(batches=batches)

        return assembler

    @classmethod
    def run_from_cli(cls, input: Union[str, Path, None], n_max_ligands: Optional[int]=None) -> 'Assembler':
        """
        Run the Assembler from a command-line context with pre/post hooks.

        Wraps run_from_yaml and integrates BaseModule CLI logging hooks.

        :param input: Path to YAML configuration file or None to use the default template.
        :type input: Union[str, Path, None]
        :param n_max_ligands: Optional override for the maximum number of ligands to load from each ligand database. Takes precedence over the value in the YAML file if specified.
        :type n_max_ligands: int | None
        :return: Assembler instance after run completion.
        :rtype: Assembler
        """
        super()._before_run_from_cli()
        super()._print_cli_input(input=input)
        assembler = cls.run_from_yaml(input=input, n_max_ligands=n_max_ligands)
        super()._after_run_from_cli()

        return assembler

    def _make_and_save_output_csv(self) -> None:
        """
        Assemble batch info into a CSV and save it via the global output controller.

        The internal df_info list is converted to a pandas DataFrame, formatted (lists -> strings)
        and written to the run info table produced by AssemblerOutput.

        :return: None
        :rtype: None
        """
        self.df_info = pd.DataFrame(self.df_info)
        self.df_info['attempt'] = self.df_info.index
        self.df_info = self.df_info[['attempt'] + [col for col in self.df_info.columns if col != 'attempt']]  # Move attempt column to front

        outdf = self.df_info.copy()
        # Make lists in the dataframe to strings for saving to csv
        for col in outdf.columns:
            if isinstance(outdf[col].iloc[0], list):
                outdf[col] = outdf[col].apply(lambda x: f'({", ".join(str(el) for el in x)})' if isinstance(x, list) else x)
        self.gbl_outcontrol.save_run_info_table(outdf)

        return

    def _final_checks(self) -> None:
        """
        Perform final consistency checks on assembled outputs.

        Verifies that no duplicate isomer names exist within each successful batch and raises
        an AssertionError if duplicates are detected.

        :return: None
        :rtype: None
        :raises AssertionError: If duplicate isomer names are found within any batch.
        """
        df_test_success = self.df_info[self.df_info['success']]
        batches = df_test_success['batch_idx'].unique()
        for batch in batches:
            df_batch = df_test_success[df_test_success['batch_idx'] == batch]
            # Check for duplicate complex names in the batch
            duplicate_names = df_batch['isomer_name'][df_batch['isomer_name'].duplicated()].values
            assert len(duplicate_names) == 0, f"Duplicate isomer names in batch {batch}: {duplicate_names}. Please report this issue to our GitHub page."

        return

    def _log_summary(self) -> None:
        """
        Log a concise per-batch and overall summary of the assembly run.

        The summary includes counts of attempted isomers, successful isomers, failed
        filters and the output directory location. Runtime is printed when verbosity > 1.

        :return: None
        :rtype: None
        """
        batch_summary_title = '  Summary per batch  '
        logging.info(f'{batch_summary_title:=^80}')
        for batch_idx, batch in enumerate(self.batches):
            batch_name = batch['name']
            df = self.df_info[self.df_info['batch_idx'] == batch_idx]
            logging.info(f"{batch_name}:")
            self._log_success_rate(df)

        # Print total summary of run
        total_summary_title = '  Total summary of DART Assembler run  '
        logging.info(f'{total_summary_title:=^80}')
        self._log_success_rate(self.df_info)
        n_isomers = self.df_info['success'].sum()
        n_complexes = self.df_info[self.df_info['success']]['graph_hash'].nunique()
        logging.info(f"DART Assembler output files saved to directory `{self.output_directory.name}`.")

        # The runtime is printed but not logged, so that slight differences in the runtime do not cause the integration tests to fail.
        if self.verbosity > 1:
            print(f"Total runtime for assembling {n_isomers} isomers (from {n_complexes} complexes): {self.runtime}")

        return

    def _log_global_info(self) -> None:
        """
        Log initial run metadata and user-defined global settings.

        Writes the configured output directory, number of batches and initialization arguments
        to the logging facility for reproducibility and debugging.

        :return: None
        :rtype: None
        """
        logging.info('Starting DART Assembler Module.')
        logging.info(f'Output directory: {self.output_directory.name}')
        plural = 'es' if self.n_batches > 1 else ''  # print plural or singular in next line
        logging.info(f"Running {self.n_batches} batch{plural}...")
        logging.info(f"User-defined global settings:")
        for key, value in self.init_args.items():
            logging.info(f"    {key: <30}{value}")

        return

    @staticmethod
    def _log_success_rate(df):
        """
        Log success/failure statistics and most common failure modes for a DataFrame of attempts.

        :param df: DataFrame produced by _make_and_save_output_csv containing per-attempt metadata.
        :type df: pd.DataFrame
        :return: None
        :rtype: None
        """
        n_total = len(df)
        n_isomers = df['success'].sum()
        n_complexes = df[df['success']]['complex_name'].nunique()

        # Output statistics how many isomers failed each filter
        post_filters = df['warning'].value_counts().to_dict()
        # Merge all warnings that start with 'duplicate' into one
        n_duplicates = sum(count for note, count in post_filters.items() if note.startswith('duplicate'))
        post_filters = {note: count for note, count in post_filters.items() if not note.startswith('duplicate')}
        if n_duplicates > 0:
            post_filters['duplicate'] = n_duplicates
        # Sort the post-filters by the number of occurrences
        post_filters = dict(sorted(post_filters.items(), key=lambda item: item[1], reverse=True))
        post_filter_notes = '\n'.join([f'    - {filter}: {n}' for filter, n in post_filters.items() if not filter == ''])

        logging.info(f"  - {n_total} isomers tried, {n_isomers} isomers (from {n_complexes} complexes) successfully assembled.")
        if post_filter_notes != '':
            logging.info(f"  - {n_total - n_isomers} isomers failed because of filters:")
            logging.info(post_filter_notes)

        return

    @staticmethod
    def _log_batch_title_and_settings(batch_settings: Dict[Any, Any]) -> None:
        """
        Log batch title and user-specified settings for that batch.

        :param batch_settings: Dictionary containing the batch configuration.
        :type batch_settings: dict
        :return: None
        :rtype: None
        """
        batch_title = f'  {batch_settings["name"]}  '
        logging.info(f'{batch_title:=^80}')
        logging.info(f"User-defined settings for {batch_settings['name']}:")
        for key, value in batch_settings.items():
            logging.info(f"    {key: <30}{value}")

        return