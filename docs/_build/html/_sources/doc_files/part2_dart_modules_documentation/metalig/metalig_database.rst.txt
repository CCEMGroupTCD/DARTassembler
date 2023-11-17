.. _metalig:

MetaLig Ligand Database
==========================

The MetaLig database is a collection of 41,018 ligands with assigned formal charges, meticulously extracted from the Cambridge Structural Database (CSD). It's a treasure trove for chemists looking to explore and design metal complexes, offering a wide range of ligand structures essential for catalysis and material design. Beyond being a rich source of information, MetaLig is closely integrated with the DART assembler module, helping users to select and assemble ligands for their specific research needs. The diversity within MetaLig is showcased in the :ref:`metalig_ligand_statistics` section.

.. figure:: /_static/part2/metalig/metalig_fig.png
   :width: 100%
   :align: center

   Figure 1: A snapshot of the MetaLig database

For chemists who prefer a hands-on approach, the database is designed to be user-friendly, with the ability to access 3D coordinates and molecular graphs through a Python interface. A concise summary of key properties is also available in a downloadable csv file, as indicated in :ref:`metalig_ligand_properties`.

**Extrinsic information:**

Unique to the MetaLig database, alongside the usual chemical information (intrinsic properties), you'll find valuable extrinsic data extracted from the CSD's extensive records. This includes statistics like how often a ligand appears in the :confval:`CSD Occurrences` or which metals it typically bonds with in the :confval:`CSD Metal Count` — information that can guide your selection of ligands for specific metal centers. Moreover, if you're looking to replicate or adapt a known synthesis, the :confval:`CSD Complex IDs` associated with each ligand point you to past literature where these ligands have been used, connecting the dots between computational design and practical laboratory work.


.. _metalig_ligand_properties:

Ligand Properties
-----------------

A suite of precomputed properties is also available for each ligand. Most properties relevant for the ligand filter module are covered here. Running the command ``dart dbinfo --path metalig`` will generate a .csv file listing the following properties:

**Intrinsic properties :**
    - **Ligand ID :** The ID of the ligand in the database.
    - **Stoichiometry :** Chemical formula of the ligand.
    - **Denticity :** The denticity of the ligand.
    - **Formal Charge :** Charge of the ligand, determined as per the DART paper.
    - **Donors :** Elements constituting the donor atoms.
    - **Number of Atoms :** Total atom count within the ligand.
    - **Molecular Weight :** Mass of the ligand in g/mol.
    - **Ligand Planarity :** Degree of planarity of all ligand atoms between 0 and 1, where 1.0 represents a perfect plane.
    - **Haptic :** Potentially haptic interactions are identified by checking if any of the coordinating atoms are bound to each other in the molecular graph.
    - **Beta-Hydrogen :** If the ligand has a hydrogen in beta position to the metal center.
    - **Max. Interatomic Distance :** Largest distance between any two atoms in the ligand, which is a measure for the size of the ligand.
    - **Avg. M-D Bond Length :** Mean of all bond lengths from metal to donor atoms.
    - **Graph ID :** The ID of the molecular graph in the database, unique for each unique ligand.

**Extrinsic properties :**
    - **CSD Occurrences :** The number of occurrences of the ligand in the CSD.
    - **CSD Complex IDs :** The IDs of the complexes in the CSD that contain the ligand.
    - **CSD Metal Count :** All metals that the ligand is coordinated to in the CSD, along with their counts.

.. _metalig_ligand_statistics:

Ligand Statistics
-----------------

.. figure:: /_static/part2/metalig/hist_donors.png
   :width: 100%
   :align: center

   Figure 2: Histogram of donor atoms in the MetaLig. For instance, there are nearly 8,000 N-N donor ligands present.

.. figure:: /_static/part2/metalig/hist_metal_center.png
   :width: 100%
   :align: center

   Figure 3: Histogram showing the prevalence of ligands coordinating to specific metals, such as over 7,000 instances of ligands which were found in the CSD coordinating to Ni.