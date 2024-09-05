.. _metalig:

MetaLig Ligand Database
==========================

.. contents:: :local:

The organo\ **Meta**\ llic **Lig**\ and database (MetaLig) contains 41,018 ligands extracted from the Cambridge Structural Database (CSD). It contains 3D coordinates, formal charge, molecular graph and a variety of physical properties. Each ligand also includes statistical data about its occurrences in the CSD, such as which metals it typically coordinates to.

The MetaLig ligand database can be used in a variety of applications:

   - **DART Assembler :** As a source of ligands for the DART Assembler module.
   - **DART Ligand Filters :** Filter ligands based on their properties to target specific chemical spaces in the Assembler module.
   - **Ligand Analysis :** To analyze and explore ligands across the CSD.
   - **Ligand Property Prediction :** As a dataset for training machine learning models.

.. figure:: /_static/part2/metalig/metalig_fig.png
   :width: 100%
   :align: center



.. _metalig_ligand_properties:

Explore Ligand Structures and Properties
------------------------------------------

To explore the ligands in the MetaLig, run the command ``DARTassembler dbinfo --path metalig``. This will generate two files, an .xyz file and a .csv file:

The .xyz file contains the 3D structures of all ligands concatenated. To view and browse through the ligands with ase, you can use the command ``ase gui concat_MetaLigDB_v1.0.0.xyz``. Each ligand is coordinated to a Cu metal center for visualization purposes.

The .csv file lists most of the physical and statistical properties of each ligand:

**Physical properties :**
    - **Ligand ID**
    - **Denticity**
    - **Donors**
    - **Stoichiometry**
    - **Number of Atoms**
    - **Formal Charge** - Determined as per the DART paper.
    - **Molecular Weight** - in g/mol.
    - **Ligand Planarity** - Degree of planarity of all ligand atoms between 0 and 1, where 1.0 represents a perfect plane.
    - **Haptic** - If the molecular graph has any neighboring donor atoms.
    - **Beta-Hydrogen** - If the ligand has a hydrogen in beta position to the metal center.
    - **Max. Interatomic Distance** - Largest distance between any two atoms in the ligand
    - **Avg. M-D Bond Length** - Mean of all bond lengths from metal to donor atoms.
    - **Graph ID** - The ID of the molecular graph in the database, unique for each unique ligand.

**Statistical CSD properties :**
    - **CSD Occurrences** - The number of occurrences of the ligand in the CSD.
    - **CSD Complex IDs** - The IDs of the complexes in the CSD that contain the ligand.
    - **CSD Metal Count** - All metals that the ligand is coordinated to in the CSD, along with their counts.

By the way, the command ``DARTassembler dbinfo --path LIGAND_DB_PATH.jsonlines`` can be used to display the same information for any
DART ligand database file in .jsonlines format. Users can generate ligand database files with the :ref:`Ligand Filters Module <ligandfilters>`.

.. _metalig_python_filtering:

Explore and Filter the MetaLig in Python
----------------------------------------------
For many users, the DART Ligand Filters module will be enough to filter ligands which exactly defined properties. For complete freedom in filtering and exploring, the MetaLig database can be accessed with the DART Python API. First, read in the MetaLig. To speed things up, let's only load the first 1000 ligands:

.. code-block:: python

    from DARTassembler.src.ligand_extraction.DataBase import LigandDB

    # Load the first 1000 out of 41,018 ligands in the MetaLig database.
    metalig = LigandDB.load_from_json(path='metalig', n_max=1000)

Now, you can filter the MetaLig database based on your requirements. For example, let's filter the MetaLig so that we retain only ligands with a formal charge of -1, with denticity of 2 and with a maximum of 50 atoms:

.. code-block:: python

    # Set some criteria to filter ligands
    keep_denticity = 2
    keep_charge = -1
    max_n_atoms = 50

    ligands_to_keep = []
    for ligand_name, ligand in metalig.db.items():
        correct_denticity = ligand.denticity == keep_denticity
        correct_charge = ligand.pred_charge == keep_charge
        correct_n_atoms = ligand.n_atoms <= max_n_atoms
        if correct_denticity and correct_charge and correct_n_atoms:
            ligands_to_keep.append(ligand_name)

    # Reduce MetaLig database to only keep ligands which adhere to the above criteria
    filtered_metalig_dict = {ligand_name: ligand for ligand_name, ligand in metalig.db.items() if ligand_name in ligands_to_keep}
    filtered_metalig = LigandDB(filtered_metalig_dict)

Now, we can save the filtered MetaLig database to a .jsonlines file.

.. code-block:: python

    filtered_metalig.save_to_file('filtered_metalig.jsonlines')

This .jsonlines file can be used in the DART Assembler module as source for ligands. Since we made a database of bidentate ligands with a formal charge of -1, we could use it to assemble e.g. neutral Fe(II) square planar complexes with geometry 2-2.
We can also save an overview table of the filtered ligand database as .csv file:

.. code-block:: python

    filtered_metalig.save_reduced_csv('filtered_metalig.csv')

This table will display 136 bidentate ligands with a formal charge of -1 and a maximum of 50 atoms.


.. _metalig_ligand_statistics:

Ligand Statistics
-----------------

.. figure:: /_static/part2/metalig/hist_donors.png
   :width: 100%
   :align: center

   Figure 2: Bar chart of donor atoms in the MetaLig. For instance, there are nearly 8,000 N-N donor ligands present.

.. figure:: /_static/part2/metalig/hist_metal_center.png
   :width: 100%
   :align: center

   Figure 3: Bar chart showing the prevalence of ligands coordinating to specific metals, such as over 7,000 instances of ligands which were found in the CSD coordinating to Ni.