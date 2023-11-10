Ligand Filter Module
====================


Introduction
------------

The Ligand Filter module allows users to filter down the entire MetaLig database with 41,018 ligands to a subset of ligands based on various criteria. Given the vast and diverse chemical space covered in the MetaLig database, this is helpful when using DART to assemble complexes with a certain chemistry, for example square-planar Pd complexes with P-N donors as shown in one of the examples.

The module is designed to provide a list of already implemented filters which can be used by specifying various options. The input to the module is a YAML file which allows users to specify various filters to apply to the ligand database. The module will then apply these filters and save the filtered ligand database to a new file. The module can be used to filter the entire MetaLig database (default) or a user-provided ligand database.

In case none of the usual filters are sufficient, the module has one extra filter (`graph_IDs`) which allows users to specify a list of ligand IDs to have exact control over which ligands make it through the filter. This allows users to select ligands by doing their Excel magic on a ligand_overview.csv file, extract the ligand IDs and input them as list into the graph_IDs filter.


Input YAML File
---------------

The Ligand Filter module uses a YAML file to allow the user to supply various options for each filter. At the beginning of each ligand filter input file, the user has to specify the input and output ligand database paths or can let them empty for default values. Then, the user can specify a list of filters after the keyword `'filters:'`. Not every filter has to be specified, and the order of the filters does not matter. Some filters which are simple yes/no switches can also be ignored by setting them to False.

Below we provide an example YAML file which can be used as a template. The template contains all available filters and examples how to use them.

Copy-Paste Ligand Filter Input File Template
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block::

    input_ligand_db_path:           # path or empty. If empty, the entire MetaLig ligand database will be used as input
    output_ligand_db_path: data/filtered_ligand_db.json   # path or empty. If empty defaults to 'filtered_ligand_db.json' in the current directory.

    filters:
      - filter: denticities         # filter out ligands with denticities not in this list
        denticities: [1, 2, 3, 4, 5]  # list of denticities to keep

      - filter: ligand_charges                            # filter out ligands with charges not provided in this list
        ligand_charges: [-1, 0]                             # only keep ligands with these charges
        apply_to_denticities: [2, 3]                                 # List of denticities to apply this filter to. If empty, will be ignored.

      - filter: ligand_composition                      # filter ligands based on their composition
        elements: [ C, H, N, P]          # elements to apply this filter to
        instruction: must_only_contain_in_any_amount    # instruction for how to apply this filter
        apply_to_denticities:                                      # List of denticities to apply this filter to. If empty, will be ignored.

      - filter: coordinating_atoms_composition            # filter ligands based on their coordinating atoms
        elements: [P, N]                       # TODO
        instruction: must_contain_and_only_contain        #
        apply_to_denticities: [2]

      - filter: number_of_atoms                  # filter out ligands with atom count outside of this range
        min: 10                               # Minimum number of atoms. If empty, will be ignored.
        max: 100                              # Maximum number of atoms. If empty, will be ignored.
        apply_to_denticities: [3]                      # List of denticities to apply this filter to. If empty, will be ignored.

      - filter: molecular_weight            # filter out ligands with molecular weight outside of this range.
        min: 30.1                             # Minimum molecular weight (todo unit). If empty, will be ignored.
        max: 300.2                            # Maximum molecular weight (todo unit). If empty, will be ignored.
        apply_to_denticities: [3]                      # List of denticities to apply this filter to. If empty, will be ignored.

      - filter: metal_donor_bond_lengths    # filter out ligands based on the bond lengths between the metal and the donor atoms
        min: 1.0                              # Minimum bond length (in Angstrom). If empty, will be ignored.
        max: 2.5                              # Maximum bond length (in Angstrom). If empty, will be ignored.
        apply_to_denticities: [2]                      # List of denticities to apply this filter to. If empty, will be ignored.

      - filter: interatomic_distances       # filter out ligands with interatomic distances outside of this range
        min: 0.5                              # Minimum interatomic distance (in Angstrom). If empty, will be ignored.
        max: 40.0                             # Maximum interatomic distance (in Angstrom). If empty, will be ignored.
        apply_to_denticities: [3]                      # List of denticities to apply this filter to. If empty, will be ignored.

      - filter: planarity                   # filter out ligands based on their 'planarity score', which is a number between 0 and 1. 0 is not planar, 1 is perfectly planar.
        min: 0.2                              # Minimum planarity score. If empty, will be ignored.
        max: 1.0                              # Maximum planarity score. If empty, will be ignored.
        apply_to_denticities: [1]                      # List of denticities to apply this filter to. If empty, will be ignored.

      - filter: occurrences                 # filter out ligands based on the number of times they have been observed in the CSD
        min: 3                                # Minimum number of occurrences. If empty, will be ignored.
        max:                                  # Maximum number of occurrences. If empty, will be ignored.
        apply_to_denticities: [4]                      # List of denticities to apply this filter to. If empty, will be ignored.

      - filter: metal_ligand_binding_history                                                  # only keep ligands which have been observed to coordinate to these metals
        metal_ligand_binding_history: [Pd, Ni]                        # list of metals to keep
        apply_to_denticities:                                                                  # List of denticities to apply this filter to. If empty, will be ignored.

      - filter: remove_ligands_with_adjacent_coordinating_atoms  # filter out ligands with neighboring coordinating atoms
        remove_ligands_with_adjacent_coordinating_atoms: True      # True or False. If False, will be ignored.

      - filter: remove_ligands_with_beta_hydrogens        # filter out ligands with beta hydrogens
        remove_ligands_with_beta_hydrogens: True            # True or False. If False, will be ignored.

      - filter: graph_IDs                   # only keep ligands with the following graph IDs
        graph_IDs: [a2b7bbb6ca4ce36dc3147760335e7374, 53b7a3d91a1be6e167a3975bb7921206]     # list of graph IDs to keep


Input/Output Options
~~~~~~~~~~~~~~~~~~~~

The following two options have to be specified at the beginning of each ligand filter input file. They can be let empty but they have to be specified.

input_ligand_db_path : [str, empty]
    Path to the input ligand database. If empty, the entire MetaLig ligand database will be used as input.
output_ligand_db_path : [str, empty]
    Path where the filtered ligand database will be saved. If empty, will default to 'filtered_ligand_db.json' in the current directory.

Filter Options
~~~~~~~~~~~~~~

denticities
^^^^^^^^^^^

Keeps only ligands with denticities specified in the list.

**Options**:
    denticities : [list[int]]
        List of denticities to keep.

**Example**:

This example will keep only ligands with denticity 2, 3 or 5.

.. code-block::

    - filter: denticities
        denticities: (2, 3, 5)


ligand_charges
^^^^^^^^^^^^^^

Keep only ligands with formal charges which are specified in the list.

**Options**:
    ligand_charges : [list[int]]
        List of charges to keep.
    denticities : [list[int], empty]
        A list of denticities. This filter will be applied only to ligands with a denticity in this list. If empty, will apply to all ligands.

**Example**:

For ligands with denticity of 2 or 3, this example will keep only ligands which have a formal charge of -1, 0 or 1. Ligands with denticities other than 2 or 3 will always pass.

.. code-block::

    - filter: ligand_charges
        ligand_charges: (-1, 0, 1)
        apply_to_denticities: (2, 3)


ligand_composition
^^^^^^^^^^^^^^^^^^

Filter ligands based on their chemical composition, i.e. the atoms in their chemical formula. The 'elements' parameter specifies the elements to apply this filter to. The 'instruction' parameter specifies exactly how to apply this filter. This filter works exactly like the 'coordinating_atoms_composition' filter, except that it applies to all atoms instead of only the coordinating atoms.

**Options**:
    elements : [list[str]]
        List of chemical elements to apply this filter to. Depending on the instruction, duplicate elements in this list may or may not be ignored.
    instruction : [str]
        Instruction for how to apply this filter. Note that instructions must be lowercase and exactly the correct string, so the best is to copy-paste it from the documentation. The following instructions are available:

        - **must_contain_and_only_contain**
            Ligands must consist of exactly these atoms in exactly this count. For example, if the 'elements' are '(C, C, H, N)', then a ligand must consist of exactly two Carbon, one Hydrogen and one Nitrogen atom to pass this filter.
        - **must_at_least_contain**
            Ligands must contain all specified elements but can also contain other elements. Duplicate elements are ignored. For example, if the 'elements' are '(C, C, H, N)', then a ligand must contain at least one Carbon, one Hydrogen and one Nitrogen atom to pass this filter.
        - **must_exclude**
            Ligands must not contain any of the specified elements. Duplicate elements are ignored. For example, if the 'elements' are '(C, C, H, N)', then a ligand must not contain any Carbon, Hydrogen or Nitrogen atoms to pass this filter.
        - **must_only_contain_in_any_amount**
            Ligands must only contain the specified elements, but the amount of each element is not important and can even be zero. Duplicate elements are ignored. For example, if the 'elements' are '(C, C, H, N)', then any ligand that contains no other elements than Carbon, Hydrogen and Nitrogen will pass this filter, and even ligands containing subsets such as ligands containing only Carbon.
    denticities : [list[int], empty]
        A list of denticities. This filter will be applied only to ligands with a denticity in this list. If empty, will apply to all ligands.

**Example**:

This example will keep only ligands with denticity 3 which consist of only Carbon, Hydrogen, Nitrogen and Phosphorus atoms or a subset of these elements. Ligands with denticities other than 3 will always pass.

.. code-block::

    - filter: ligand_composition
        elements: (C, H, N, P)
        instruction: must_only_contain_in_any_amount
        apply_to_denticities: (3)



coordinating_atoms_composition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Filter ligands based on their coordinating atoms, i.e. the atoms bound to the metal center. The 'elements' parameter specifies the elements to apply this filter to. The 'instruction' parameter specifies exactly how to apply this filter. This filter works exactly like the 'ligand_composition' filter, except that it only applies to the coordinating atoms of the ligand.

**Options**:
    elements : [list[str]]
        List of chemical elements to apply this filter to. Depending on the instruction, duplicate elements in this list may or may not be ignored.
    instruction : [str]
        Instruction for how to apply this filter. Note that instructions must be lowercase and exactly the correct string, so the best is to copy-paste it from the documentation. The following instructions are available:

        - **must_contain_and_only_contain**
            The ligand must have exactly these coordinating atoms in exactly this count. For example, if the 'elements' are '(C, C, N)', the ligand must have exactly two Carbon and one Nitrogen atom coordinating to the metal.
        - **must_at_least_contain**
            The coordinating atoms of the ligand must contain all specified elements but can also contain other elements. Duplicate elements are ignored. For example, if the 'elements' are '(C, C, N)', then the list of coordinating atoms must contain at least one Carbon and one Nitrogen atom to pass this filter.
        - **must_exclude**
            The coordinating atoms of the ligand must not contain any of the specified elements. Duplicate elements are ignored. For example, if the 'elements' are '(C, C, N)', then the list of coordinating atoms must not contain any Carbon or Nitrogen atoms to pass this filter.
        - **must_only_contain_in_any_amount**
            The coordinating atoms of the ligand must only contain the specified elements, but the amount of each element is not important and can even be zero. Duplicate elements are ignored. For example, if the 'elements' are '(C, C, N)', then any ligand with coordinating atoms which contain no other elements than Carbon and Nitrogen will pass this filter, and even ligands containing subsets such as ligands containing only Carbon.
    denticities : [list[int], empty]
        A list of denticities. This filter will be applied only to ligands with a denticity in this list. If empty, will apply to all ligands.

**Example**:

This example will keep only ligands with denticity of 3 which have exactly one Carbon, one Nitrogen and one Oxygen coordinating to the metal center. Ligands with denticities other than 3 will be removed automatically, since these will always have more or less coordinating atoms.

.. code-block::

    - filter: coordinating_atoms_composition
        elements: (C, N, O)
        instruction: must_contain_and_only_contain
        apply_to_denticities:


number_of_atoms
^^^^^^^^^^^^^^^

Removes ligands with number of atoms outside of the specified range. The 'min' and 'max parameters specify the minimum and maximum number of atoms, respectively.

**Options**:
    min : [float, empty]
        Minimum number of atoms. If empty, will be set to 0.
    max : [float, empty]
        Maximum number of atoms. If empty, will be set to infinity.
    denticities : [list[int], empty]
        A list of denticities. This filter will be applied only to ligands with a denticity in this list. If empty, will apply to all ligands.

**Example**:

This example will remove all ligands with a denticity of 1 or 2 with less than 10 atoms or more than 100 atoms. Ligands with denticities other than 1 or 2 will always pass.

.. code-block::

    - filter: number_of_atoms
        min: 10
        max: 100
        apply_to_denticities: (1, 2)


molecular_weight
^^^^^^^^^^^^^^^^

Only keeps ligands with molecular weight within the specified range. The 'min' and 'max' parameters specify the minimum and maximum molecular weight, respectively. For example, setting 'min' to 30 and 'max' to 300 will remove all ligands with molecular weight less than 30g/mol or more than 300g/mol.

**Options**:
    min : [float, empty]
        Minimum molecular weight in g/mol. If empty, will be set to 0.
    max : [float, empty]
        Maximum molecular weight in g/mol. If empty, will be set to infinity.
    denticities : [list[int], empty]
        A list of denticities. This filter will be applied only to ligands with a denticity in this list. If empty, will apply to all ligands.

**Example**:

This example will keep only ligands with a molecular weight between 10g/mol and 300g/mol. Because the denticities list is empty, this filter will be applied to every ligand.

.. code-block::

    - filter: molecular_weight
        min: 30
        max: 300
        apply_to_denticities:


metal_donor_bond_lengths
^^^^^^^^^^^^^^^^^^^^^^^^

Only keeps ligands with metal-donor bond lengths within the specified range. All bond lengths between the metal and the donor atoms are considered. The 'min' and 'max' parameters specify the minimum and maximum allowed bond length for at least one bond.

**Options**:
    min : [float, empty]
        Minimum bond length in Angstrom. If empty, will be set to 0.
    max : [float, empty]
        Maximum bond length in Angstrom. If empty, will be set to infinity.
    denticities : [list[int], empty]
        A list of denticities. This filter will be applied only to ligands with a denticity in this list. If empty, will apply to all ligands.

**Example**:

For ligands with a denticity of 2 or 3, this example will only keep ligands which have a metal-donor bond length between 1.0 Angstrom and 2.5 Angstrom. Ligands with denticities other than 2 or 3 will always pass.

.. code-block::

    - filter: metal_donor_bond_lengths
        min: 1.0
        max: 2.5
        apply_to_denticities: (2, 3)


interatomic_distances
^^^^^^^^^^^^^^^^^^^^^

Only keeps ligands with interatomic distances within the specified range. The calculated interatomic distances are not only between atoms with a bond, but between all atoms in the ligand. The maximum interatomic distance is a measure for the total size of the ligand, while the minimum interatomic distance is a measure for the smallest bond length. Therefore, this filter is basically a 2-in-1 filter which can be used to remove either too big ligands or ligands with too small bond lengths.

**Options**:
    min : [float, empty]
        Minimum interatomic distance in Angstrom. If empty, will be set to 0.
    max : [float, empty]
        Maximum interatomic distance in Angstrom. If empty, will be set to infinity.
    denticities : [list[int], empty]
        A list of denticities. This filter will be applied only to ligands with a denticity in this list. If empty, will apply to all ligands.

**Example**:

For ligands with a denticity of 3 or 4, this example will only keep ligands which have an interatomic distance between 0.5 Angstrom and 40 Angstrom. Ligands with denticities other than 3 or 4 will always pass.

.. code-block::

    - filter: interatomic_distances
        min: 0.5
        max: 40
        apply_to_denticities: (3, 4)

planarity
^^^^^^^^^

This filter uses a 'planarity score' to filter ligands based on how planar all their atoms are. Very planar ligands are ones in which all atoms lie in one plane, while very non-planar ligands are ones which are sphere-like. The planarity score is a number between 0 and 1, where 0 is not planar (a perfect sphere) and 1 is perfectly planar. Because this planarity score has no physical intuition behind it, it is recommended to try different values for the 'min' and 'max' parameters to see what works best for your application.

**Options**:
    min : [float]
        Minimum planarity score. If empty, will be set to 0.
    max : [float]
        Maximum planarity score. If empty, will be set to 1.
    denticities : [list[int], empty]
        A list of denticities. This filter will be applied only to ligands with a denticity in this list. If empty, will apply to all ligands.

**Example**:

This example will keep only ligands with a denticity of 1 which have a planarity score between 0.9 and 1.0, i.e. very planar ligands. Ligands with denticities other than 1 will always pass.

.. code-block::

    - filter: planarity
        min: 0.9
        max: 1
        apply_to_denticities: (1)

occurrences
^^^^^^^^^^^

Filters ligands based on how often they were observed in the Cambridge Structural Database (CSD). The 'min' and 'max' parameters specify the minimum and maximum number of occurrences, respectively.

**Options**:
    min : [int]
        Minimum number of occurrences. If empty, will be set to 0.
    max : [int]
        Maximum number of occurrences. If empty, will be set to infinity.
    denticities : [list[int], empty]
        A list of denticities. This filter will be applied only to ligands with a denticity in this list. If empty, will apply to all ligands.

**Example**:

For ligands with denticities of 3 or 4, this example will keep only ligands which have been observed in the CSD at least 3 times. Ligands with denticities other than 3 or 4 will always pass.

.. code-block::

    - filter: occurrences
        min: 3
        max:
        apply_to_denticities: (3, 4)


metal_ligand_binding_history
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Keep only ligands which have been observed in the Cambridge Structural Database to coordinate to the metals specified in the 'metal_ligand_binding_history' list. If a ligand has never been observed coordinating to any of the metals in the 'metal_ligand_binding_history' list, it will be filtered out.

**Options**:
    metal_ligand_binding_history : [list[str]]
        List of metals, e.g. (Pd, Ni). Any metal from the d- or f-block can be specified.
    denticities : [list[int], empty]
        A list of denticities. This filter will be applied only to ligands with a denticity in this list. If empty, will apply to all ligands.

**Example**:

For ligands with denticity of 2 or 3, this example will keep only ligands which have been observed to coordinate to Pd or Ni. Ligands with denticities other than 2 or 3 will always pass.

.. code-block::

    - filter: metal_ligand_binding_history
        metal_ligand_binding_history: (Pd, Ni)
        apply_to_denticities: (2, 3)



remove_ligands_with_adjacent_coordinating_atoms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Removes ligands that have coordinating atoms with a bond between them, i.e. coordinating atoms which are neighbors.

**Options**:
    remove_ligands_with_adjacent_coordinating_atoms : [True or False]
        If True, apply this filter. If False, will be ignored.

**Example**:

This example will remove all ligands with neighboring coordinating atoms.

.. code-block::

      - filter: remove_ligands_with_adjacent_coordinating_atoms
            remove_ligands_with_adjacent_coordinating_atoms: True

remove_ligands_with_beta_hydrogens
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Removes ligands with beta Hydrogen atoms, i.e. Hydrogen atoms bound to coordinating atoms.

**Options**:
    remove_ligands_with_beta_hydrogens : [True or False]
        If True, apply this filter. If False, will be ignored.

**Example**:

This example will remove all ligands with beta Hydrogen atoms.

.. code-block::

      - filter: remove_ligands_with_beta_hydrogens
            remove_ligands_with_beta_hydrogens: True


graph_IDs
^^^^^^^^^

A filter to keep only the exactly specified ligands. Graph IDs are unique IDs for each ligand which can be taken from the ligand overview csv. This filter will remove all other ligands except for the ones specified.

**Options**:
    graph_IDs : [list[str]]
        List of graph IDs to keep.

**Example**:

This example will keep only the 2 ligands with the graph IDs 'a2b7bbb6ca4ce36dc3147760335e7374' and '53b7a3d91a1be6e167a3975bb7921206'.

.. code-block::

    - filter: graph_IDs
        graph_IDs: [a2b7bbb6ca4ce36dc3147760335e7374, 53b7a3d91a1be6e167a3975bb7921206]


Table: overview of all available filters
----------------------------------------

.. csv-table:: Overview of all available filters
   :file: ligand_filters_overview.csv
   :widths: 20, 20, 10, 50
   :header-rows: 1
