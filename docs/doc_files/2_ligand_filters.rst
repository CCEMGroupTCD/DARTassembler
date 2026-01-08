.. _ligandfilters:

LigandFilters Module
========================

.. contents:: :local:

LigandFilters Input
-----------------------------------

The LigandFilters Module enables users to obtain a set of ligands with well-defined properties from the entire :ref:`MetaLig Database <metalig>`. These filters are invaluable for assembling complexes targeted to a user-defined chemical space.

There are four types of filters, which can be applied to any of the :ref:`MetaLig ligand properties <metalig_properties>`:

- :ref:`property <property_filter>` : filter by a simple named property such as :filter:`charge`, :filter:`archetype` or :filter:`n_haptic_groups`.
- :ref:`composition <composition_filter>` : filter by composition of donor atoms or the entire ligand.
- :ref:`smarts <smarts_filter>` : filter by a SMARTS pattern to exclude/include ligands with specific substructures, e.g. ``'[N]=[N]'`` for azo groups.
- :ref:`parents <parents_filter>` : filter by known parent metal centers of the ligand, e.g. ``['Pt2+', 'Pd']``.

The LigandFilters module is run in the terminal by providing a single configuration file:

.. code-block:: bash

    DARTassembler ligandfilters --input ligandfilters.yml

**Copy-Paste Template:**

.. literalinclude:: ../../DARTassembler/data/default/ligandfilters.yml
   :language: yaml
   :linenos:

Users can download this template into their current directory as ``ligandfilters.yml`` by running:

.. code-block:: bash

    DARTassembler configs --outdir .

LigandFilters Output
-----------------------------------

The output of the LigandFilters module is a file called e.g. ``filtered_ligand_db.jsonlines``, which contains all ligands that passed the specified filters. These files can be used as input ``ligand_db_files`` for the :ref:`Assembler Module <assembler>`. Thus, users can assemble complexes with ligands that have exactly the desired properties.

If ``dbinfo`` is ``True``, an additional directory ``info_filtered_ligand_db/`` is created, which contains summary information about the filtered ligands, including:

    - ``ligands_overview.csv``: A summary table of the filters. The column ``filter`` indicates which filter each ligand failed. If a ligand passed all filters, this column contains ``'Passed'``. The passed ligands are at the very end of the table.
    - ``filters.txt``: This is the main log file of the filtering process, recording all messages, warnings, and errors during filtering.
    - ``concat_xyz/``: This folder contains concatenated .xyz files of all ligands that failed a certain filter. Each file is named according to the filter it corresponds to. All the ligands which passed all filters are stored in ``concat_Passed.xyz``. These concatenated .xyz files can easily be browsed using the ``ase gui`` command from the ASE package:

    .. code-block:: bash

        ase gui concat_Passed.xyz

LigandFilters Options
--------------------------------

The provided filters are applied in the order they are listed in the ``filters`` list. Each filter is a dictionary with a key ``'filter'`` specifying the filter type and other filter-specific keys. Each filter is then passed to the corresponding method of the :class:`DARTassembler.src.metalig.mol.Ligand` class. Please refer to the docstrings of these 4 methods for a detailed description of all available options. You will see that the options match perfectly with the options in the .yml configuration template above.

**Usage:** To use these filters on a database, it is recommended to use the terminal command ``DARTassembler ligandfilters --input ligandfilters.yml`` and edit the provided configuration file as needed. The python API below is mainly to document the available options in a single place and for users who want to write their own filtering scripts in Python.

.. _property_filter:
.. autofunction:: DARTassembler.src.metalig.mol.Ligand.property_filter
   :no-index:

.. _composition_filter:
.. autofunction:: DARTassembler.src.metalig.mol.Ligand.composition_filter
   :no-index:

.. _smarts_filter:
.. autofunction:: DARTassembler.src.metalig.mol.Ligand.smarts_filter
   :no-index:

.. _parents_filter:
.. autofunction:: DARTassembler.src.metalig.mol.Ligand.parents_filter
   :no-index:
