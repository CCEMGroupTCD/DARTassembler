.. _ligandfilters:

Ligand Filters Module
========================

.. contents:: :local:

The Ligand Filters Module enables users to obtain a set of ligands with well-defined properties from the extensive :ref:`MetaLig Database <metalig>`. These filters are invaluable for assembling complexes targeted to a user-defined chemical space.

Users can apply a large range of predefined filters. For those requiring precise control over the structures, the :confval:`smarts` filter allows for the application of powerful SMARTS patterns to filter ligands based on their 2D chemical structure. Furthermore, the :confval:`graph_IDs` filter enables the selection of individual ligands. Alternatively to using the Ligand Filters Module with pre-defined filters, users can also :ref:`use Python to explore the MetaLig and create custom filters <metalig_python_filtering>`.

The ligand filters module is run in the terminal by providing a single configuration file:

.. code-block:: bash

    DARTassembler ligandfilters --path ligandfilters_input.yml

The following filters are currently implemented:

    Physical Property Filters:
        - :confval:`denticities`
        - :confval:`ligand_charges`
        - :confval:`ligand_composition`
        - :confval:`coordinating_atoms_composition`
        - :confval:`metal_donor_bond_lengths`
        - :confval:`number_of_atoms`
        - :confval:`molecular_weight`
        - :confval:`interatomic_distances`
        - :confval:`planarity`
    Molecular Graph Filters:
        - :confval:`remove_ligands_with_adjacent_coordinating_atoms`
        - :confval:`remove_ligands_with_beta_hydrogens`
        - :confval:`remove_ligands_with_missing_bond_orders`
        - :confval:`atomic_neighbors`
        - :confval:`smarts`
    Statistical CSD Filters:
        - :confval:`occurrences`
        - :confval:`metal_ligand_binding_history`
    Special Filters:
        - :confval:`graph_IDs`

Input File
~~~~~~~~~~

Users interact with the Ligand Filters Module by providing an input file in YAML format. In this file, users can specify parameters for each filter, repeat the same filter with different parameters, or omit filters they don't need. The order of filters doesn't matter.

This template specifies all available filters and examples of their parameters:

**Copy-Paste Template:**

.. code-block:: yaml

    ################## Settings for the DART ligand filters module. ##################
    # Everything after '#' is ignored by the program and only there for the user.

    input_db_file: test_metalig                     # path, 'metalig' or 'test_metalig'. Default: 'metalig'
    output_db_file: filtered_ligand_db.jsonlines    # path. Default: 'filtered_ligand_db.jsonlines'
    output_ligands_info: true                       # true or false. If true, an overview of the filtered and passed ligands will be saved. Default: true

    filters:

      ####### Physical Property Filters #######

      - filter: denticities
        denticities: [2, 3, 4]                      # Only keep ligands with these denticities

      - filter: ligand_charges
        ligand_charges: [-1, 0, 1]                  # Only keep ligands with these charges
        apply_to_denticities:                       # List of denticities to apply this filter to. If empty, applies to all denticities.

      - filter: ligand_composition                  # Filters ligands by their stoichiometry
        elements: CHN                               # Stoichiometry/list of elements to apply this filter to
        instruction: must_only_contain_in_any_amount    # Instruction for how to apply this filter. Options: 'must_contain_and_only_contain', 'must_at_least_contain', 'must_exclude', 'must_only_contain_in_any_amount'
        apply_to_denticities:                       # List of denticities to apply this filter to. If empty, applies to all denticities.

      - filter: coordinating_atoms_composition      # Filters ligands by their donor atoms
        elements: CN                                # Stoichiometry/list of elements to apply this filter to
        instruction: must_contain_and_only_contain  # Instruction for how to apply this filter. Options: 'must_contain_and_only_contain', 'must_at_least_contain', 'must_exclude', 'must_only_contain_in_any_amount'
        apply_to_denticities:                       # List of denticities to apply this filter to. If empty, applies to all denticities.

      - filter: metal_donor_bond_lengths            # Filters ligands by the bond lengths between the metal and the donor atoms (in Angstrom).
        min: 1.3                                    # If empty, defaults to 0.
        max: 2.0                                    # If empty, defaults to infinity.
        apply_to_denticities: [2]                   # List of denticities to apply this filter to. If empty, applies to all denticities.

      - filter: number_of_atoms                     # Filters ligands by their total atom count.
        min: 10                                     # If empty, defaults to 0.
        max: 100                                    # If empty, defaults to infinity.
        apply_to_denticities: [1]                   # List of denticities to apply this filter to. If empty, applies to all denticities.

      - filter: molecular_weight                    # Filters ligands by their molecular weight (in g/mol).
        min:                                        # If empty, defaults to 0.
        max: 200                                    # If empty, defaults to infinity.
        apply_to_denticities:                       # List of denticities to apply this filter to. If empty, applies to all denticities.

      - filter: interatomic_distances               # Filters ligands by interatomic distances (in Angstrom), but not only bonds.
        min: 0.6                                    # If empty, defaults to 0.
        max:                                        # If empty, defaults to infinity.
        apply_to_denticities:                       # List of denticities to apply this filter to. If empty, applies to all denticities.

      - filter: planarity                           # The 'planarity score' is a number between 0 and 1. 1 means all ligand atoms are perfectly planar.
        min: 0.9                                    # If empty, defaults to 0.
        max: 1.0                                    # If empty, defaults to 1.0.
        apply_to_denticities:                       # List of denticities to apply this filter to. If empty, applies to all denticities.

      ####### Molecular Graph Filters #######

      - filter: remove_ligands_with_adjacent_coordinating_atoms     # Filter out ligands with neighboring coordinating atoms
        remove_ligands_with_adjacent_coordinating_atoms: true       # true or false. If false, filter will have no effect.
        apply_to_denticities:                                       # List of denticities to apply this filter to. If empty, applies to all denticities.

      - filter: remove_ligands_with_beta_hydrogens                  # Filter out ligands with beta hydrogens
        remove_ligands_with_beta_hydrogens: true                    # true or false. If false, filter will have no effect.
        apply_to_denticities:                                       # List of denticities to apply this filter to. If empty, applies to all denticities.

      - filter: remove_ligands_with_missing_bond_orders             # Filter out ligands with missing bond orders
        remove_ligands_with_missing_bond_orders: true               # true or false. If false, filter will be ignored.
        apply_to_denticities:                                       # List of denticities to apply this filter to. If empty, applies to all denticities.

      - filter: atomic_neighbors                    # Filters out ligands in which a chemical element is connected to the specified neighbors
        atom: C                                     # Chemical element of the central atom
        neighbors: H2                               # List of chemical elements/stoichiometry of the neighbors
        apply_to_denticities:                       # List of denticities to apply this filter to. If empty, applies to all denticities.

      - filter: smarts                              # Filter ligands using SMARTS patterns. Recommended to be used with filter:remove_ligands_with_missing_bond_orders
        smarts: '[C&H2]'                            # SMARTS pattern to match. Important: use single quotes around the SMARTS pattern.
        should_contain: false                       # If true, the ligand must contain the SMARTS pattern to pass the filter. If false, the ligand must not contain the SMARTS pattern to pass.
        include_metal: false                        # If true, the ligand structure will contain a 'Cu' metal center connected to the coordinating atoms when matching the SMARTS pattern.
        apply_to_denticities:                       # List of denticities to apply this filter to. If empty, applies to all denticities.

      ####### Statistical CSD Filters #######

      - filter: occurrences                         # Filter out ligands based on the number of times they have been observed in the CSD
        min: 20                                     # If empty, defaults to 0.
        max:                                        # If empty, defaults to infinity.
        apply_to_denticities:                       # List of denticities to apply this filter to. If empty, applies to all denticities.

      - filter: metal_ligand_binding_history        # Only keep ligands which have been observed to coordinate to these metals
        metal_ligand_binding_history: [Pd, Ni]      # List of metals to keep
        apply_to_denticities:                       # List of denticities to apply this filter to. If empty, applies to all denticities.

      ####### Special Filters #######

      - filter: graph_IDs                           # Only keep ligands with specified graph IDs
        graph_IDs: [a2b7bbb6ca4ce36dc3147760335e7374, 53b7a3d91a1be6e167a3975bb7921206]     # List of graph IDs to keep



You can also download this template into your current directory by running:

.. code-block:: bash

    DARTassembler configs --path .

.. Note::

    Every filter, except :confval:`denticities` and :confval:`graph_IDs` filter, has an optional parameter **apply_to_denticities**. This parameter allows users to apply the respective filter only to ligands with the specified denticities, which can be very handy. If this parameter is empty or omitted from the file, the filter will be applied to all ligands.

Global Options
~~~~~~~~~~~~~~~~~~~~

The following options specify global settings for the Ligand Filters Module. If a setting is missing, the default value is used.

.. confval:: input_db_file

    :type: `filepath`, ``metalig``, ``test_metalig``
    :default: ``metalig``

    Path to the input ligand database. If empty, the entire :ref:`MetaLig ligand database<metalig>` will be used as input.

.. confval:: output_db_file

    :type: `filepath`
    :default: ``filtered_ligand_db.jsonlines``

    Path to where the filtered ligand database will be saved.

.. confval:: output_ligands_info

    :type: ``true``, ``false``
    :default: ``true``

    If ``false``, only the ligand database file will be saved. If ``true``, a directory with info files about the database and the filtering process will be saved.

Physical Property Filters
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _filter_denticities:

.. confval:: denticities

    Keeps only ligands with denticities specified in the list.

    :options:

        denticities :
            List of denticities to keep.

    :example: This example will keep only ligands with denticity 2, 3 and 5.

        .. code-block:: yaml

            - filter: denticities
              denticities: [2, 3, 5]

.. _filter_ligand_charges:

.. confval:: ligand_charges

    Keep only ligands with formal charges specified in the list.

    :options:

        ligand_charges :
            List of formal charges to keep.

        apply_to_denticities :
            Denticity or list of denticities. This filter will be applied only to ligands with the specified denticities. If empty or omitted, will apply to all ligands.

    :example: For ligands with denticity of 2 or 3, this example will keep only ligands which have a formal charge of -1, 0 or 1. Ligands with denticities other than 2 or 3 will always pass.

        .. code-block:: yaml

            - filter: ligand_charges
              ligand_charges: [-1, 0, 1]
              apply_to_denticities: [2, 3]

.. _filter_ligand_composition:

.. confval:: ligand_composition

    Filter ligands based on their chemical composition, e.g. C\ :sub:`6`\H\ :sub:`5` for phenyl. The filter has four different modes: depending on the value of **instruction**, the specified **elements** are used to check a different condition. This filter works exactly like the :confval:`coordinating_atoms_composition` filter, except that it applies to all atoms of the ligand instead of only the set of coordinating atoms.

    :options:

        **elements :**

            Stoichiometry or list of chemical elements to apply this filter to. For example, specifying ``CH2N`` is equivalent to ``[C, H, H, N]``. For most instructions, the atom count is irrelevant and only the specified elements are used by the filter.

        **instruction :**

            Instruction for how to apply this filter. The following instructions are available:

            - ``must_contain_and_only_contain``
                Ligands must consist of exactly these atoms in exactly this count. Use this to filter for exact stoichiometry.
            - ``must_at_least_contain``
                Ligands must contain all specified elements but can also contain other elements. Atom count is ignored, only elements are important.
            - ``must_exclude``
                Ligands must not contain any of the specified elements. Atom count is ignored, only elements are important.
            - ``must_only_contain_in_any_amount``
                Ligands must contain no other elements than the specified elements, but may contain not all of the specified elements. Atom count is ignored, only elements are important.

        **apply_to_denticities :**

            Denticity or list of denticities. This filter will be applied only to ligands with the specified denticities. If empty or omitted, will apply to all ligands.

    :example: This will keep only ligands with exact stoichiometry of C\ :sub:`2`\H\ :sub:`6`\N.

        .. code-block:: yaml

            - filter: ligand_composition
              elements: C2H6N
              instruction: must_contain_and_only_contain
              apply_to_denticities:

    :example: This will keep only ligands which contain at least the elements C, H, N and may contain other elements.

        .. code-block:: yaml

            - filter: ligand_composition
              elements: CHN
              instruction: must_at_least_contain
              apply_to_denticities:

    :example: This will keep only ligands which do not contain any C, H or N atoms.

        .. code-block:: yaml

            - filter: ligand_composition
              elements: CHN
              instruction: must_exclude
              apply_to_denticities:

    :example: This will keep only ligands which contain C, H, N or subsets of these elements (e.g. C, H or only H).
    
        .. code-block:: yaml
    
            - filter: ligand_composition
              elements: CHN
              instruction: must_only_contain_in_any_amount
              apply_to_denticities:


.. _filter_coordinating_atoms_composition:

.. confval:: coordinating_atoms_composition

    Filter ligands based on their donor atoms. The filter has four different modes: depending on the value of **instruction**, the specified **elements** are used to check a different condition. This filter works exactly like the :confval:`ligand_composition` filter, except that it applies only to the set of donor atoms instead of all atoms in the ligand.

    :options:

        **elements :**

            Stoichiometry or list of chemical elements to apply this filter to. For example, specifying ``N2`` is equivalent to ``[N, N]``. For most instructions, the atom count is irrelevant and only the specified elements are used by the filter.

        **instruction :**

            Instruction for how to apply this filter. The following instructions are available:

            - ``must_contain_and_only_contain``
                Donor atoms must consist of exactly these atoms in exactly this count. Use this to filter for an exact list of donor atoms, e.g. N-N ligands.
            - ``must_at_least_contain``
                Donor atoms must contain all specified elements but can also contain other elements. Atom count is ignored, only elements are important.
            - ``must_exclude``
                Donor atoms must not contain any of the specified elements. Atom count is ignored, only elements are important.
            - ``must_only_contain_in_any_amount``
                Donor atoms must contain no other elements than the specified elements, but may contain not all of the specified elements. Atom count is ignored, only elements are important.

        **apply_to_denticities :**

            Denticity or list of denticities. This filter will be applied only to ligands with the specified denticities. If empty or omitted, will apply to all ligands.

    :example: This will keep only bidentate N-N donors.

        .. code-block:: yaml

            - filter: coordinating_atoms_composition
              elements: N2
              instruction: must_contain_and_only_contain
              apply_to_denticities:

    :example: This will keep only ligands which coordinate via at least one C and one N atom, such as C-N or C-N-H donors.

        .. code-block:: yaml

            - filter: coordinating_atoms_composition
              elements: CN
              instruction: must_at_least_contain
              apply_to_denticities:

    :example: This will keep only ligands which do not coordinate via any C or N atoms, such as O-O donors.

        .. code-block:: yaml

            - filter: coordinating_atoms_composition
              elements: CN
              instruction: must_exclude
              apply_to_denticities:

    :example: This will keep only ligands which coordinate only via C and N atoms or subsets of these atoms, such as C-N-N or N-N donors.

        .. code-block:: yaml

            - filter: coordinating_atoms_composition
              elements: CN
              instruction: must_only_contain_in_any_amount
              apply_to_denticities:

.. tip::

    The :confval:`ligand_composition` and :confval:`coordinating_atoms_composition` filters have four different modes depending on the **instruction** parameter. On first glance, these modes might seem too general to make a useful filter, but by combining the same filter multiple times with different instructions, users can achieve very specific filters.

.. _filter_metal_donor_bond_lengths:

.. confval:: metal_donor_bond_lengths

    Only keeps ligands where all metal-donor bond lengths are within the specified range.

    :options:

        min :
            Minimum bond length in Angstrom. If empty, will be set to 0.

        max :
            Maximum bond length in Angstrom. If empty, will be treated as infinity.

        apply_to_denticities :

            Denticity or list of denticities. This filter will be applied only to ligands with the specified denticities. If empty or omitted, will apply to all ligands.

    :example: This filter would remove a bidentate ligand with metal-donor bond lengths of (1.4, 2.2) Angstrom, but keep another bidentate ligand with metal-donor bond lengths of (1.6, 1.8) Angstrom.

        .. code-block:: yaml

            - filter: metal_donor_bond_lengths
              min: 1.3
              max: 2.0
              apply_to_denticities: [2]

.. _filter_number_of_atoms:

.. confval:: number_of_atoms

    Removes ligands with number of atoms outside of the specified range.

    :options:

        min :
            Minimum number of atoms. If empty, will be set to 0.

        max :
            Maximum number of atoms. If empty, will be treated as infinity.

        apply_to_denticities :
            Denticity or list of denticities. This filter will be applied only to ligands with the specified denticities. If empty or omitted, will apply to all ligands.

    :example: This example will remove all monodentate ligands with less than 10 atoms or more than 100 atoms. Ligands with denticities other than 1 will always pass.

        .. code-block:: yaml

            - filter: number_of_atoms
              min: 10
              max: 100
              apply_to_denticities: [1]

.. _filter_molecular_weight:

.. confval::  molecular_weight

    Only keeps ligands with molecular weight within the specified range.

    :options:

        min :
            Minimum molecular weight in g/mol. If empty, will be set to 0.

        max :
            Maximum molecular weight in g/mol. If empty, will be treated as infinity.

        apply_to_denticities :
            Denticity or list of denticities. This filter will be applied only to ligands with the specified denticities. If empty or omitted, will apply to all ligands.

    :example: This example will keep only ligands with a maximum molecular weight of 200 g/mol.

        .. code-block:: yaml

            - filter: molecular_weight
              min:
              max: 200
              apply_to_denticities:

.. _filter_interatomic_distances:

.. confval:: interatomic_distances

    Only keeps ligands in which all interatomic distances are within the specified range. The calculated interatomic distances are not only between atoms with a bond, but between all atoms in the ligand, even far-away ones! The maximum interatomic distance can be used as a measure for the size of a ligand, while the minimum interatomic distance can be used as a measure for how close atoms are in the ligand. Therefore, this filter is basically a 2-in-1 filter which can be used to remove ligands which are either too big or have atoms which are too close to each other.

    :options:

        min :
            Minimum interatomic distance in Angstrom. If empty, will be set to 0.

        max :
            Maximum interatomic distance in Angstrom. If empty, will be treated as infinity.

        apply_to_denticities :
            Denticity or list of denticities. This filter will be applied only to ligands with the specified denticities. If empty or omitted, will apply to all ligands.

    :example: This filter will remove ligands if any two atoms in the ligand are closer than 0.6 Angstrom.

        .. code-block:: yaml

            - filter: interatomic_distances
              min: 0.6
              max:
              apply_to_denticities:

    :example: This filter will remove "big" ligands which are more than 30 Angstroms long in any direction, without considering bulkiness.

        .. code-block:: yaml

            - filter: interatomic_distances
              min:
              max: 30
              apply_to_denticities:

.. _filter_planarity:

.. confval:: planarity

    This filter uses a 'planarity score' to filter ligands based on how planar all their atoms are. Very planar ligands are ones in which all atoms lie in one plane, while very non-planar ligands are ones which are sphere-like. The planarity score is a number between 0 and 1, where 0 is not planar (a perfect sphere) and 1 is perfectly planar. Because this planarity score has no physical intuition behind it, it is recommended to try different values and see what works best for your application.

    :options:

        min :
            Minimum planarity score. If empty, will be set to 0.

        max :
            Maximum planarity score. If empty, will be set to 1.

        apply_to_denticities :
            Denticity or list of denticities. This filter will be applied only to ligands with the specified denticities. If empty or omitted, will apply to all ligands.

    :example: This filter will keep only relatively planar ligands in which most atoms lie mostly in the same plane.

        .. code-block:: yaml

            - filter: planarity
              min: 0.9
              max: 1
              apply_to_denticities:

.. tip::

    There are four filters which can be used as a measure for the size and bulkiness of a ligand: :confval:`number_of_atoms`, :confval:`molecular_weight`,  :confval:`interatomic_distances` and :confval:`planarity`. They all measure different aspects and can be used in combination to define the dimension of your ligands.

Molecular Graph Filters
~~~~~~~~~~~~~~~~~~~~~~~~

.. _filter_remove_ligands_with_adjacent_coordinating_atoms:

.. confval:: remove_ligands_with_adjacent_coordinating_atoms

    Removes ligands that have a donor atom bonding to another donor atom, which often correlates with haptic interactions. It is recommended to always apply this filter because DART in its current version cannot assemble these ligands yet and they are filtered out during the assembly anyway.

    :options:

        remove_ligands_with_adjacent_coordinating_atoms :
            If ``true``, apply this filter. If ``false``, this filter has no effect.

        apply_to_denticities :
            Denticity or list of denticities. This filter will be applied only to ligands with the specified denticities. If empty or omitted, will apply to all ligands.

    :example: This example will remove all ligands with neighboring coordinating atoms.

        .. code-block:: yaml

              - filter: remove_ligands_with_adjacent_coordinating_atoms
                remove_ligands_with_adjacent_coordinating_atoms: true
                apply_to_denticities:

.. _filter_remove_ligands_with_beta_hydrogens:

.. confval:: remove_ligands_with_beta_hydrogens

    Removes ligands with beta hydrogen atoms, i.e. hydrogen atoms bound to donor atoms.

    :options:

        remove_ligands_with_beta_hydrogens :
            If ``true``, apply this filter. If ``false``, this filter has no effect.

        apply_to_denticities :
            Denticity or list of denticities. This filter will be applied only to ligands with the specified denticities. If empty or omitted, will apply to all ligands.

    :example: This example will remove all ligands with beta hydrogen atoms.

        .. code-block:: yaml

              - filter: remove_ligands_with_beta_hydrogens
                remove_ligands_with_beta_hydrogens: true
                apply_to_denticities:

.. _filter_remove_ligands_with_missing_bond_orders:

.. confval:: remove_ligands_with_missing_bond_orders

    Removes ligands with missing bond orders (~4% of ligands in the MetaLig). Most helpful in concert with the filter :confval:`smarts`, since that filter will automatically pass ligands with unknown bond orders. If you want to be sure that all passed ligands obey the SMARTS filter, it is recommended to apply this filter together with the SMARTS filter.

    :options:

        remove_ligands_with_missing_bond_orders :
            If ``true``, apply this filter. If ``false``, this filter has no effect.

        apply_to_denticities :
            Denticity or list of denticities. This filter will be applied only to ligands with the specified denticities. If empty or omitted, will apply to all ligands.

    :example: This example will remove all ligands with missing bond orders.

        .. code-block:: yaml

              - filter: remove_ligands_with_missing_bond_orders
                remove_ligands_with_missing_bond_orders: true
                apply_to_denticities:

.. _filter_atomic_neighbors:

.. confval:: atomic_neighbors

        This filter removes all ligands in which a chemical element :confval:`atom` is connected to the atoms specified in :confval:`neighbors`. Importantly, this filter only checks if the specified atom has at least the specified neighbors, but there might be more neighbors than specified and the ligand will still be removed. For more control, use the :confval:`smarts` filter.

        :options:

            **atom :**

                Chemical element of the central atom.

            **neighbors :**

                List of chemical elements or stoichiometry. The ligand will be removed if the :confval:`atom` is connected to at least the specified neighbors.

            **apply_to_denticities :**

                Denticity or list of denticities. This filter will be applied only to ligands with the specified denticities. If empty or omitted, will apply to all ligands.

        :example: This example removes all ligands in which a C is connected to 2 H atoms, plus potentially other neighbors.

            .. code-block:: yaml

                - filter: atomic_neighbors
                  atom: C
                  neighbors: H2
                  apply_to_denticities:

.. _filter_smarts:

.. confval:: smarts

        This filter is a very powerful tool to filter ligands based on their 2D chemical structure, including bond orders. `SMARTS <https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html>`_ is a language to describe and match chemical patterns and motifs in molecules. It can be thought of as a way to search chemical motifs in SMILES strings.

        The smarts filter works by first computing the SMILES string of the ligand (with or without 'Cu' metal center depending on the parameter **include_metal**) and then matching the specified SMARTS pattern to the SMILES string using rdkit.

        .. warning::
            If a ligand has unknown bond orders (~4% of ligands in the MetaLig), it will automatically pass this filter. If you want to be sure that all passed ligands obey the SMARTS filter, it is recommended to apply this filter together with the filter :confval:`remove_ligands_with_missing_bond_orders`.

        .. note::
            SMARTS patterns are very expressive, but can be difficult to come up with. We recommended to use tools like `SMARTSviewer <https://smartsview.zbh.uni-hamburg.de/>`_ to design your SMARTS pattern. We have also made very good experiences with using Large Language Models like ChatGPT. Either way, always make sure your SMARTS pattern works as intended by checking the passed and failed output ligands of the filter.

        :options:

            **smarts :**

                `SMARTS <https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html>`_ pattern to match. Please note that the SMARTS pattern must be enclosed in single or double quotes, e.g. '[C&H2]'. Otherwise it is likely that the YAML parser will throw an error.

            **should_contain :**

                If ``true``, the ligand `must contain` the SMARTS pattern to pass. If ``false``, the ligand `must not contain` the SMARTS pattern to pass.

            **include_metal :**

                If ``true``, the ligand's coordinating atoms will be connected to a Cu metal center. The bonds between Cu and the coordinating atoms are defined as single bonds. This allows to target coordinating atoms in the SMARTS pattern in contrast to other atoms. If ``false``, the ligand will be treated as just the ligand structure without a metal center.

            **apply_to_denticities :**

                Denticity or list of denticities. This filter will be applied only to ligands with the specified denticities. If empty or omitted, will apply to all ligands.

        :example: This example will remove all ligands in which any C atom bonds to exactly 2 H atoms.

            .. code-block:: yaml

                - filter: smarts
                  smarts: '[C&H2]'
                  should_contain: false
                  include_metal: false
                  apply_to_denticities:


Statistical CSD Filters
~~~~~~~~~~~~~~~~~~~~~~~~

.. _filter_occurrences:

.. confval:: occurrences

    Filters ligands based on how often they were observed in the Cambridge Structural Database (CSD).

    :options:

        min :
            Minimum number of occurrences. If empty, will be set to 0.

        max :
            Maximum number of occurrences. If empty, will be treated as infinity.

        apply_to_denticities :
            Denticity or list of denticities. This filter will be applied only to ligands with the specified denticities. If empty or omitted, will apply to all ligands.

    :example: This example will keep only ligands which have been observed in the CSD at least 20 times. This might be helpful to avoid exotic ligands and help with synthetic feasibility.

        .. code-block:: yaml

            - filter: occurrences
              min: 20
              max:
              apply_to_denticities:


.. _filter_metal_ligand_binding_history:

.. confval:: metal_ligand_binding_history

    Keep only ligands which have been observed in the Cambridge Structural Database to coordinate to specific metals. If a ligand has never been observed coordinating to any of the specified metals it will be filtered out.

    :options:

        metal_ligand_binding_history :
            List of metals, e.g. [Pd, Ni]. Any metal from the d- or f-block can be specified.

        apply_to_denticities :
            Denticity or list of denticities. This filter will be applied only to ligands with the specified denticities. If empty or omitted, will apply to all ligands.

    :example: This filter will keep only ligands which have been observed to coordinate to Pd or Ni.

        .. code-block:: yaml

            - filter: metal_ligand_binding_history
              metal_ligand_binding_history: [Pd, Ni]
              apply_to_denticities:


Special Filters
~~~~~~~~~~~~~~~~~~~~~~~~

.. _filter_graph_IDs:

.. confval:: graph_IDs

    A filter to keep only individually specified ligands. Graph IDs are unique IDs for each ligand which can be found in all ligand .csv files, generated e.g. by the :ref:`dbinfo module <module_overview>`. Together with :ref:`writing custom filters using python <metalig_python_filtering>`, this filter is very useful for special requirements.

    :options:

        graph_IDs :
            List of graph IDs of the ligands to keep.

    :example: This example will keep only the 2 ligands with the graph IDs `a2b7bbb6ca4ce36dc3147760335e7374` and `53b7a3d91a1be6e167a3975bb7921206`.

        .. code-block:: yaml

            - filter: graph_IDs
              graph_IDs: [a2b7bbb6ca4ce36dc3147760335e7374, 53b7a3d91a1be6e167a3975bb7921206]