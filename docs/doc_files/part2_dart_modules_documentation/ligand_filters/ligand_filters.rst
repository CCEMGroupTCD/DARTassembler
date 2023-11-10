.. _ligandfilters:

Ligand Filter Module
====================

The Ligand Filter module streamlines the selection of ligands from the extensive :ref:`MetaLig Database <metalig>`, which boasts 41,018 entries. This tool is invaluable for assembling complexes targeted to a certain chemical space.

Users can apply a range of predefined filters via an input file to refine the ligand pool. The module supports both the entire MetaLig database and user-supplied subsets. For those requiring precise control, the `graph_IDs` filter enables selection based on specific ligand IDs.

**Currently implemented filters :**

- :ref:`filter_denticities`
- :ref:`filter_ligand_charges`
- :ref:`filter_ligand_composition`
- :ref:`filter_coordinating_atoms_composition`
- :ref:`filter_number_of_atoms`
- :ref:`filter_molecular_weight`
- :ref:`filter_metal_donor_bond_lengths`
- :ref:`filter_interatomic_distances`
- :ref:`filter_planarity`
- :ref:`filter_occurrences`
- :ref:`filter_metal_ligand_binding_history`
- :ref:`filter_remove_ligands_with_adjacent_coordinating_atoms`
- :ref:`filter_remove_ligands_with_beta_hydrogens`
- :ref:`filter_graph_IDs`


Input File
---------------

The module's input is a YAML file where users define the input and output paths for the ligand databases, followed by the desired filters under the `filters:` keyword. Filters are customizable; users can omit unwanted filters and the order of filters doesn't matter.

Below is an example YAML file, serving as a comprehensive template showcasing all available filters and their usage:


.. code-block::

    input_ligand_db_path:                           # Path or empty. If empty, the entire MetaLig ligand database will be used as input
    output_ligand_db_path: filtered_ligand_db.json  # Path or empty. If empty defaults to 'filtered_ligand_db.json' in the current directory.

    filters:
      - filter: denticities
        denticities: [1, 2, 3, 4, 5]                # Only keep ligands with these denticities

      - filter: ligand_charges
        ligand_charges: [-1, 0]                     # Only keep ligands with these charges
        apply_to_denticities: [2, 3]                # List of denticities to apply this filter to. If empty, applies to all denticities.

      - filter: ligand_composition
        elements: [C, H, N, P]                      # Elements to apply this filter to
        instruction: must_only_contain_in_any_amount    # Instruction for how to apply this filter
        apply_to_denticities:                       # List of denticities to apply this filter to. If empty, applies to all denticities.

      - filter: coordinating_atoms_composition
        elements: [P, N]                            # Elements to apply this filter to
        instruction: must_contain_and_only_contain  # Instruction for how to apply this filter
        apply_to_denticities: [2]                   # List of denticities to apply this filter to. If empty, applies to all denticities.

      - filter: number_of_atoms                     # Filters ligands by their total atom count.
        min: 10                                     # Minimum number of atoms. If empty, defaults to 0.
        max: 100                                    # Maximum number of atoms. If empty, defaults to infinity.
        apply_to_denticities: [3]                   # List of denticities to apply this filter to. If empty, applies to all denticities.

      - filter: molecular_weight
        min: 30.1                                   # Minimum molecular weight (in g/mol). If empty, defaults to 0.
        max: 300.2                                  # Maximum molecular weight (in g/mol). If empty, defaults to infinity.
        apply_to_denticities: [3]                   # List of denticities to apply this filter to. If empty, applies to all denticities.

      - filter: metal_donor_bond_lengths
        min: 1.0                                    # Minimum bond length (in Angstrom). If empty, defaults to 0.
        max: 2.5                                    # Maximum bond length (in Angstrom). If empty, defaults to infinity.
        apply_to_denticities: [2]                   # List of denticities to apply this filter to. If empty, applies to all denticities.

      - filter: interatomic_distances
        min: 0.5                                    # Minimum interatomic distance (in Angstrom). If empty, defaults to 0.
        max: 40.0                                   # Maximum interatomic distance (in Angstrom). If empty, defaults to infinity.
        apply_to_denticities: [3]                   # List of denticities to apply this filter to. If empty, applies to all denticities.

      - filter: planarity                           # The 'planarity score' is a number between 0 and 1. 0 is not planar, 1 is perfectly planar.
        min: 0.2                                    # Minimum planarity score. If empty, defaults to 0.
        max: 1.0                                    # Maximum planarity score. If empty, defaults to 1.0.
        apply_to_denticities: [1]                   # List of denticities to apply this filter to. If empty, applies to all denticities.

      - filter: occurrences                         # Filter out ligands based on the number of times they have been observed in the CSD
        min: 3                                      # Minimum number of occurrences. If empty, defaults to 0.
        max:                                        # Maximum number of occurrences. If empty, defaults to infinity.
        apply_to_denticities: [4]                   # List of denticities to apply this filter to. If empty, applies to all denticities.

      - filter: metal_ligand_binding_history        # Only keep ligands which have been observed to coordinate to these metals
        metal_ligand_binding_history: [Pd, Ni]      # List of metals to keep
        apply_to_denticities:                       # List of denticities to apply this filter to. If empty, applies to all denticities.

      - filter: remove_ligands_with_adjacent_coordinating_atoms     # Filter out ligands with neighboring coordinating atoms
        remove_ligands_with_adjacent_coordinating_atoms: true       # true or false. If false, filter will be ignored. Recommended to set to true.

      - filter: remove_ligands_with_beta_hydrogens                  # Filter out ligands with beta hydrogens
        remove_ligands_with_beta_hydrogens: true                    # true or false. If false, filter will be ignored.

      - filter: graph_IDs                           # Only keep ligands with the following graph IDs
        graph_IDs: [a2b7bbb6ca4ce36dc3147760335e7374, 53b7a3d91a1be6e167a3975bb7921206]     # List of graph IDs to keep


Input/Output Options
~~~~~~~~~~~~~~~~~~~~

The following two options have to be specified at the beginning of each ligand filter input file. They can be let empty but they have to be specified.

input_ligand_db_path :
    Path to the input ligand database. If empty, the entire MetaLig ligand database will be used as input.

output_ligand_db_path :
    Path to where the filtered ligand database will be saved. If empty, will default to 'filtered_ligand_db.json' in the current directory.

Filters
~~~~~~~~~~~~~~

.. _filter_denticities:

denticities
^^^^^^^^^^^

Keeps only ligands with denticities specified in the list.


denticities :
    List of denticities to keep.

**Example** :    This example will keep only ligands with denticity 2, 3 or 5.

.. code-block::

    - filter: denticities
        denticities: [2, 3, 5]

.. _filter_ligand_charges:

ligand_charges
^^^^^^^^^^^^^^

Keep only ligands with formal charges which are specified in the list.


ligand_charges :
    List of charges to keep.

denticities :
    A list of denticities. This filter will be applied only to ligands with a denticity in this list. If empty, will apply to all ligands.

**Example** :    For ligands with denticity of 2 or 3, this example will keep only ligands which have a formal charge of -1, 0 or 1. Ligands with denticities other than 2 or 3 will always pass.

.. code-block::

    - filter: ligand_charges
        ligand_charges: [-1, 0, 1]
        apply_to_denticities: [2, 3]

.. _filter_ligand_composition:

ligand_composition
^^^^^^^^^^^^^^^^^^

Filter ligands based on their chemical composition, i.e. the atoms in their chemical formula. The 'elements' parameter specifies the elements to apply this filter to. The 'instruction' parameter specifies exactly how to apply this filter. This filter works exactly like the 'coordinating_atoms_composition' filter, except that it applies to all atoms instead of only the coordinating atoms.


elements :
    List of chemical elements to apply this filter to. Depending on the instruction, duplicate elements in this list may or may not be ignored.

instruction :
    Options: any of the words below

    Instruction for how to apply this filter. Note that instructions must be lowercase and exactly the correct string, so the best is to copy-paste it from the documentation. The following instructions are available:

    - **must_contain_and_only_contain**
        Ligands must consist of exactly these atoms in exactly this count. For example, if the 'elements' are '[C, C, H, N]', then a ligand must consist of exactly two Carbon, one Hydrogen and one Nitrogen atom to pass this filter.
    - **must_at_least_contain**
        Ligands must contain all specified elements but can also contain other elements. Duplicate elements are ignored. For example, if the 'elements' are '[C, C, H, N]', then a ligand must contain at least one Carbon, one Hydrogen and one Nitrogen atom to pass this filter.
    - **must_exclude**
        Ligands must not contain any of the specified elements. Duplicate elements are ignored. For example, if the 'elements' are '[C, C, H, N]', then a ligand must not contain any Carbon, Hydrogen or Nitrogen atoms to pass this filter.
    - **must_only_contain_in_any_amount**
        Ligands must only contain the specified elements, but the amount of each element is not important and can even be zero. Duplicate elements are ignored. For example, if the 'elements' are '[C, C, H, N]', then any ligand that contains no other elements than Carbon, Hydrogen and Nitrogen will pass this filter, and even ligands containing subsets such as ligands containing only Carbon.

denticities :
    A list of denticities. This filter will be applied only to ligands with a denticity in this list. If empty, will apply to all ligands.

**Example** :    This example will keep only ligands with denticity 3 which consist of only Carbon, Hydrogen, Nitrogen and Phosphorus atoms or a subset of these elements. Ligands with denticities other than 3 will always pass.

.. code-block::

    - filter: ligand_composition
        elements: [C, H, N, P]
        instruction: must_only_contain_in_any_amount
        apply_to_denticities: [3]

.. _filter_coordinating_atoms_composition:

coordinating_atoms_composition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Filter ligands based on their coordinating atoms, i.e. the atoms bound to the metal center. The 'elements' parameter specifies the elements to apply this filter to. The 'instruction' parameter specifies exactly how to apply this filter. This filter works exactly like the 'ligand_composition' filter, except that it only applies to the coordinating atoms of the ligand.

elements :
    List of chemical elements to apply this filter to. Depending on the instruction, duplicate elements in this list may or may not be ignored.

instruction :
    Options: any of the words below

    Instruction for how to apply this filter. Note that instructions must be lowercase and exactly the correct string, so the best is to copy-paste it from the documentation. The following instructions are available:

    - **must_contain_and_only_contain**
        The ligand must have exactly these coordinating atoms in exactly this count. For example, if the 'elements' are '[C, C, N]', the ligand must have exactly two Carbon and one Nitrogen atom coordinating to the metal.
    - **must_at_least_contain**
        The coordinating atoms of the ligand must contain all specified elements but can also contain other elements. Duplicate elements are ignored. For example, if the 'elements' are '[C, C, N]', then the list of coordinating atoms must contain at least one Carbon and one Nitrogen atom to pass this filter.
    - **must_exclude**
        The coordinating atoms of the ligand must not contain any of the specified elements. Duplicate elements are ignored. For example, if the 'elements' are '[C, C, N]', then the list of coordinating atoms must not contain any Carbon or Nitrogen atoms to pass this filter.
    - **must_only_contain_in_any_amount**
        The coordinating atoms of the ligand must only contain the specified elements, but the amount of each element is not important and can even be zero. Duplicate elements are ignored. For example, if the 'elements' are '[C, C, N]', then any ligand with coordinating atoms which contain no other elements than Carbon and Nitrogen will pass this filter, and even ligands containing subsets such as ligands containing only Carbon.

denticities :
    A list of denticities or empty. This filter will be applied only to ligands with a denticity in this list. If empty, will apply to all ligands.

**Example** : This example will keep only ligands with denticity of 3 which have exactly one Carbon, one Nitrogen and one Oxygen coordinating to the metal center. Ligands with denticities other than 3 will be removed automatically, since these will always have more or less coordinating atoms.

.. code-block::

    - filter: coordinating_atoms_composition
        elements: [C, N, O]
        instruction: must_contain_and_only_contain
        apply_to_denticities:

.. _filter_number_of_atoms:

number_of_atoms
^^^^^^^^^^^^^^^

Removes ligands with number of atoms outside of the specified range. The 'min' and 'max parameters specify the minimum and maximum number of atoms, respectively.


min :
    Minimum number of atoms. If empty, will be set to 0.

max :
    Maximum number of atoms. If empty, will be set to infinity.

denticities :
    A list of denticities or empty. This filter will be applied only to ligands with a denticity in this list. If empty, will apply to all ligands.

**Example** : This example will remove all ligands with a denticity of 1 or 2 with less than 10 atoms or more than 100 atoms. Ligands with denticities other than 1 or 2 will always pass.

.. code-block::

    - filter: number_of_atoms
        min: 10
        max: 100
        apply_to_denticities: [1, 2]

.. _filter_molecular_weight:

molecular_weight
^^^^^^^^^^^^^^^^

Only keeps ligands with molecular weight within the specified range. The 'min' and 'max' parameters specify the minimum and maximum molecular weight, respectively. For example, setting 'min' to 30 and 'max' to 300 will remove all ligands with molecular weight less than 30g/mol or more than 300g/mol.

min :
    Minimum molecular weight in g/mol. If empty, will be set to 0.

max :
    Maximum molecular weight in g/mol. If empty, will be set to infinity.

denticities :
    A list of denticities or empty. This filter will be applied only to ligands with a denticity in this list. If empty, will apply to all ligands.

**Example** : This example will keep only ligands with a molecular weight between 10g/mol and 300g/mol. Because the denticities list is empty, this filter will be applied to every ligand.

.. code-block::

    - filter: molecular_weight
        min: 30
        max: 300
        apply_to_denticities:

.. _filter_metal_donor_bond_lengths:

metal_donor_bond_lengths
^^^^^^^^^^^^^^^^^^^^^^^^

Only keeps ligands with metal-donor bond lengths within the specified range. All bond lengths between the metal and the donor atoms are considered. The 'min' and 'max' parameters specify the minimum and maximum allowed bond length for at least one bond.

min :
    Minimum bond length in Angstrom. If empty, will be set to 0.

max :
    Maximum bond length in Angstrom. If empty, will be set to infinity.

denticities :
    A list of denticities or empty. This filter will be applied only to ligands with a denticity in this list. If empty, will apply to all ligands.

**Example** : For ligands with a denticity of 2 or 3, this example will only keep ligands which have a metal-donor bond length between 1.0 Angstrom and 2.5 Angstrom. Ligands with denticities other than 2 or 3 will always pass.

.. code-block::

    - filter: metal_donor_bond_lengths
        min: 1.0
        max: 2.5
        apply_to_denticities: [2, 3]

.. _filter_interatomic_distances:

interatomic_distances
^^^^^^^^^^^^^^^^^^^^^

Only keeps ligands with interatomic distances within the specified range. The calculated interatomic distances are not only between atoms with a bond, but between all atoms in the ligand. The maximum interatomic distance is a measure for the total size of the ligand, while the minimum interatomic distance is a measure for the smallest bond length. Therefore, this filter is basically a 2-in-1 filter which can be used to remove either too big ligands or ligands with too small bond lengths.

min :
    Minimum interatomic distance in Angstrom. If empty, will be set to 0.

max :
    Maximum interatomic distance in Angstrom. If empty, will be set to infinity.

denticities :
    A list of denticities or empty. This filter will be applied only to ligands with a denticity in this list. If empty, will apply to all ligands.

**Example** : For ligands with a denticity of 3 or 4, this example will only keep ligands which have an interatomic distance between 0.5 Angstrom and 40 Angstrom. Ligands with denticities other than 3 or 4 will always pass.

.. code-block::

    - filter: interatomic_distances
        min: 0.5
        max: 40
        apply_to_denticities: [3, 4]

.. _filter_planarity:

planarity
^^^^^^^^^

This filter uses a 'planarity score' to filter ligands based on how planar all their atoms are. Very planar ligands are ones in which all atoms lie in one plane, while very non-planar ligands are ones which are sphere-like. The planarity score is a number between 0 and 1, where 0 is not planar (a perfect sphere) and 1 is perfectly planar. Because this planarity score has no physical intuition behind it, it is recommended to try different values for the 'min' and 'max' parameters to see what works best for your application.


min :
    Minimum planarity score. If empty, will be set to 0.

max :
    Maximum planarity score. If empty, will be set to 1.

denticities :
    A list of denticities or empty. This filter will be applied only to ligands with a denticity in this list. If empty, will apply to all ligands.

**Example** : This example will keep only ligands with a denticity of 1 which have a planarity score between 0.9 and 1.0, i.e. very planar ligands. Ligands with denticities other than 1 will always pass.

.. code-block::

    - filter: planarity
        min: 0.9
        max: 1
        apply_to_denticities: [1]

.. _filter_occurrences:

occurrences
^^^^^^^^^^^

Filters ligands based on how often they were observed in the Cambridge Structural Database (CSD).

min :
    Minimum number of occurrences. If empty, will be set to 0.

max :
    Maximum number of occurrences. If empty, will be set to infinity.

denticities :
    A list of denticities or empty. This filter will be applied only to ligands with a denticity in this list. If empty, will apply to all ligands.

**Example** : For ligands with denticities of 3 or 4, this example will keep only ligands which have been observed in the CSD at least 3 times. Ligands with denticities other than 3 or 4 will always pass.

.. code-block::

    - filter: occurrences
        min: 3
        max:
        apply_to_denticities: [3, 4]


.. _filter_metal_ligand_binding_history:

metal_ligand_binding_history
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Keep only ligands which have been observed in the Cambridge Structural Database to coordinate to the metals specified in the 'metal_ligand_binding_history' list. If a ligand has never been observed coordinating to any of the metals in the 'metal_ligand_binding_history' list, it will be filtered out.

metal_ligand_binding_history :
    List of metals, e.g. [Pd, Ni]. Any metal from the d- or f-block can be specified.

denticities :
    A list of denticities or empty. This filter will be applied only to ligands with a denticity in this list. If empty, will apply to all ligands.

**Example** :   For ligands with denticity of 2 or 3, this example will keep only ligands which have been observed to coordinate to Pd or Ni. Ligands with denticities other than 2 or 3 will always pass.

.. code-block::

    - filter: metal_ligand_binding_history
        metal_ligand_binding_history: [Pd, Ni]
        apply_to_denticities: [2, 3]

.. _filter_remove_ligands_with_adjacent_coordinating_atoms:

remove_ligands_with_adjacent_coordinating_atoms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Removes ligands that have coordinating atoms with a bond between them, i.e. coordinating atoms which are neighbors. It is recommended to apply this filter, since it filters out ligands with haptic interactions, which are difficult to assemble and might not be stable.

remove_ligands_with_adjacent_coordinating_atoms :
    If true, apply this filter. If false, will be ignored.

**Example** : This example will remove all ligands with neighboring coordinating atoms.

.. code-block::

      - filter: remove_ligands_with_adjacent_coordinating_atoms
            remove_ligands_with_adjacent_coordinating_atoms: true

.. _filter_remove_ligands_with_beta_hydrogens:

remove_ligands_with_beta_hydrogens
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Removes ligands with beta Hydrogen atoms, i.e. Hydrogen atoms bound to coordinating atoms.


remove_ligands_with_beta_hydrogens :
    If true, apply this filter. If false, will be ignored.

**Example** : This example will remove all ligands with beta Hydrogen atoms.

.. code-block::

      - filter: remove_ligands_with_beta_hydrogens
            remove_ligands_with_beta_hydrogens: true

.. _filter_graph_IDs:

graph_IDs
^^^^^^^^^

A filter to keep only the exactly specified ligands. Graph IDs are unique IDs for each ligand which can be taken from the ligand overview csv. This filter will remove all other ligands except for the ones specified.
This filter allows users to select ligands by doing their Excel magic on a ligand_overview.csv file generated by ``dart dbinfo``, extract the ligand IDs and input them as list into the graph_IDs filter.

graph_IDs :
    List of graph IDs to keep.

**Example** : This example will keep only the 2 ligands with the graph IDs `a2b7bbb6ca4ce36dc3147760335e7374` and `53b7a3d91a1be6e167a3975bb7921206`.

.. code-block::

    - filter: graph_IDs
        graph_IDs: [a2b7bbb6ca4ce36dc3147760335e7374, 53b7a3d91a1be6e167a3975bb7921206]



