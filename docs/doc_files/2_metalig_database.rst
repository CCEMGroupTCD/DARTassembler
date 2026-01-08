.. _metalig:

MetaLig Ligand Database
==========================

.. contents:: :local:

The organo\ **Meta**\ llic **Lig**\ and database (MetaLig) contains 41,018 ligands extracted from 107,185 complexes from the Cambridge Structural Database (CSD). It contains 3D coordinates, formal charge, molecular graph and a variety of physical properties. Each ligand also includes statistical data about its occurrences in the CSD, such as which metals it typically coordinates to.

The MetaLig ligand database can be used in a variety of applications:

   - **DART Assembler:** As a source of ligands for the DART Assembler module.
   - **DART LigandFilters:** Filter ligands based on their properties to target specific chemical spaces in the Assembler module.
   - **Ligand Analysis:** To analyze and explore ligands across the CSD.
   - **Ligand Property Prediction:** As a dataset for training machine learning models.

.. figure:: /_static/part2/metalig/metalig_fig.png
   :width: 100%
   :align: center

To explore the ligands in the MetaLig, use the terminal to run the command

.. code-block:: bash

    DARTassembler dbinfo --db metalig

This will generate two files, an .xyz file and a .csv file:

The .xyz file contains the 3D structures of all ligands. To view and browse through the ligands with ase, you can use the command ``ase gui concat_MetaLigDB_v1.1.0.xyz``. Each ligand is coordinated to a Cu metal center for visualization purposes. The Cu metal center is not part of the ligands in the MetaLig, it is only added to the .xyz file to display the coordination of each ligand.

The ``MetaLigDB_v1.1.0.csv`` file displays a tabular overview of all ligands and their properties (see below). You can open this file with a program like Excel to sort and filter the ligands based on their properties.

.. _metalig_properties:

MetaLig Ligand Properties
-------------------------------

The MetaLig contains a wide range of properties, such as the 3D geometry and the molecular graph of each ligand. You can find all properties under the :class:`Ligand <DARTassembler.src.metalig.mol.Ligand>` class documentation. The MetaLig also contains 38 tabular properties for each ligand. 29 of these are useful for filtering ligands in the DART :ref:`LigandFilters <ligandfilters>` module:

.. csv-table::
   :file: ../../DARTassembler/data/docs/ligand_filter_properties.csv
   :header-rows: 1
   :widths: 15, 5, 5, 40, 35
   :align: center

Another 9 properties are mostly useful for inspection and analysis:

.. csv-table::
   :file: ../../DARTassembler/data/docs/other_properties.csv
   :header-rows: 1
   :widths: 25, 5, 40, 30
   :align: center

.. _metalig_python_filtering:


Filter the MetaLig in Python
----------------------------------------------
For many users, the DART Ligand Filters module will be enough to filter ligands with exactly defined properties. For complete freedom in filtering and exploring, the MetaLig database can be accessed via the DARTassembler Python API, specifically the :class:`LigandDB <DARTassembler.src.metalig.db.LigandDB>` and :class:`Ligand <DARTassembler.src.metalig.mol.Ligand>` classes. This allows you to write your own custom filtering scripts in Python to target ligands with exactly the properties you need.

As an example, let us extract Cp-like ligands from the MetaLig database. First, read in the MetaLig. To speed things up in this example, let's only load the first 5000 ligands.

.. code-block:: python

    from DARTassembler import LigandDB

    # Load the first 1000 out of 41,018 ligands in the MetaLig database.
    metalig = LigandDB.from_json(path='metalig', n_max=5000)

Now, you can filter the MetaLig database based on your requirements. For example, let's filter the MetaLig so that we retain only Cp-like ligands with an archetype of ``1-mono``, a charge of -1 and 5 C donor atoms. For more information how to use the Ligand objects from the MetaLig see :ref:`its documentation <ligands>`.

.. code-block:: python

    # Set some criteria to filter Cp-like ligands
    archetype = '1-mono'
    charge = -1
    donor_elements = ['C', 'C', 'C', 'C', 'C']

    # Filter ligands and keep only those which adhere to all the above criteria
    ligands_to_keep = []
    for ligand_name, ligand in metalig.db.items():
        correct_denticity = ligand.archetype == archetype
        correct_charge = ligand.charge == charge
        correct_donor_elements = ligand.donor_elements == donor_elements
        if correct_denticity and correct_charge and correct_donor_elements:
            ligands_to_keep.append(ligand_name)

    # Reduce MetaLig database to only keep ligands which adhere to the above criteria
    filtered_metalig = metalig.get_sub_db(ligand_names=ligands_to_keep)
    print(f'Number of ligands after filtering: {len(filtered_metalig.db)}')

Now, we can save the filtered MetaLig database to a .jsonlines file and a concatenated .xyz file.

.. code-block:: python

    filtered_metalig.save_to_file('filtered_metalig.jsonlines')
    filtered_metalig.save_to_concat_xyz('filtered_metalig.xyz')
    filtered_metalig.save_to_csv('filtered_metalig.csv')


This .jsonlines file can be used in the DART Assembler module as source for ligands. By opening the .csv file with a program like Excel, you will see that this table displays 7 ligands Cp-like ligand with a formal charge of -1. You can also inspect the ligand structures in the concatenated .xyz file using ``ase gui filtered_metalig.xyz``. In this way, you can use Python to filter the MetaLig database to your exact requirements and then save the filtered database to a .jsonlines file for use in the DART Assembler module.

.. _metalig_ligand_statistics:

Ligand Statistics
-----------------

.. figure:: /_static/part2/metalig/hist_donors.png
   :width: 100%
   :align: center

   Bar chart of donor atoms in the MetaLig. For instance, there are nearly 8,000 N-N donor ligands present.

.. figure:: /_static/part2/metalig/hist_metal_center.png
   :width: 100%
   :align: center

   Bar chart showing the prevalence of ligands coordinating to specific metals, such as over 8,000 instances of ligands which were found in the CSD coordinating to Cu.
