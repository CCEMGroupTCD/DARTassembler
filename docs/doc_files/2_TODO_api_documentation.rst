API Documentation
=====================

Ligands
^^^^^^^^^

This section gives an overview of the ligands in the MetaLig database. Each ligand is extracted from a complex in the Cambridge Structural Database as explained in our DART paper.


Ligand properties
-----------------

Each ligand in the MetaLig database has the following properties:

`atomic_props`:
    A dictionary containing the atomic properties of the ligand. The keys are the atomic symbols and the values are lists of the atomic properties. The order of the atomic properties is the same for all atoms. The atomic properties are:

    - `x`: The x-coordinate of the atom in the complex.
    - `y`: The y-coordinate of the atom in the complex.
    - `z`: The z-coordinate of the atom in the complex.
    - `atoms`: The atomic symbol of the atom.
    - `original_complex_indices`: The index of the atom in the complex.

`global_props`:
    A dictionary containing the global properties of the ligand. The keys are the property names and the values are the property values. The global properties are:

    - `CSD_code`: The CSD code of the complex.
    - `n_atoms`: The number of atoms in the ligand.
    - `n_elements`: The number of different elements in the ligand.
    - `molecular_weight`: The molecular weight of the ligand.
    - `n_C_H_bonds`: The number of C-H bonds in the ligand.
    - `LCS_pred_charge`: The predicted charge of the ligand.
    - `LCS_pred_charge_confidence`: The confidence of the predicted charge.
    - `LCS_pred_charge_is_confident`: Whether the predicted charge is confident.
    - `LCS_pred_charge_exact`: The exact predicted charge.

`mol`:
    The ASE Atoms object of the ligand.

`n_atoms`:
    The number of atoms in the ligand.

`n_hydrogens`:
    The number of hydrogen atoms in the ligand.

`n_protons`:
    The number of protons in the ligand.

`graph`:
    The molecular graph of the ligand, i.e. which atoms are connected to which atoms.

`graph_hash`:
    The Weisfeiler-Lehman hash of the molecular graph of the ligand.

`n_bonds`:
    The number of bonds in the ligand.

`has_bond_order_attribute`:
    Whether the bonds in the ligand have a bond order attribute.

`has_unknown_bond_orders`:
    Whether the ligand has bonds with unknown bond orders.

`has_good_bond_orders`:
    Whether the ligand has bond orders and all of them are known.

`heavy_atoms_graph_hash`:
    The Weisfeiler-Lehman hash of the molecular graph of the ligand without the hydrogen atoms.

`bond_order_graph_hash`:
    The Weisfeiler-Lehman hash of the molecular graph of the ligand with the bond orders.

`stoichiometry`:
    The stoichiometry of the ligand. If an element has a count of 1, this 1 is explicitly written.

`original_complex_id`:
    The CSD code of the complex.

`local_elements`:
    The chemical symbols of the coordinating atoms.

`was_connected_to_metal`:
    Whether the ligand was connected to a metal.

`original_metal`:
    The atomic number of the original metal.

`original_metal_position`:
    The position of the metal in the complex.

`original_metal_symbol`:
    The chemical symbol of the metal.

`original_metal_os`:
    The oxidation state of the metal.

`is_centrosymmetric`:
    Whether the ligand is centrosymmetric.

`centrosymmetry_ang_dev`:
    The angular deviation from centrosymmetry.

`graph_hash_with_metal`:
    The Weisfeiler-Lehman hash of the molecular graph of the ligand with the metal. This is the property used to find identical ligands.

`heavy_atoms_graph_hash_with_metal`:
    The Weisfeiler-Lehman hash of the molecular graph of the ligand without the hydrogen atoms with the metal.

`has_betaH`:
    Whether the ligand has a beta hydrogen.

`has_neighboring_coordinating_atoms`:
    Whether the ligand has neighboring coordinating atoms, i.e. haptic interactions.

`stats`:
    A dictionary containing the statistics of the ligand. The keys are the property names and the values are the property values. The statistics are:

    - `min_distance_to_metal`: The minimum distance of an atom to the metal.
    - `max_distance_to_metal`: The maximum distance of an atom to the metal.
    - `coordinating_atom_distances_to_metal`: The distances of the coordinating atoms to the metal.
    - `min_atomic_distance`: The minimum distance between two atoms in the ligand.
    - `max_atomic_distance`: The maximum distance between two atoms in the ligand.
    - `max_dist_per_atoms`: The maximum distance between two atoms in the ligand divided by the number of atoms in the ligand.

`unique_name`:
    The unique name of the ligand.

`pred_charge`:
    The predicted formal charge of the ligand.

`pred_charge_is_confident`:
    Whether the formal charge prediction was confident or not.

`all_ligand_names`:
    A list containing all names of

`identical_ligand_info`:
    A dictionary containing the information of identical ligands. The keys are the property names and the values are the property values. The identical ligand information is:

    - `name`: The name of the ligand.
    - `original_metal_symbol`: The chemical symbol of the metal.
    - `original_metal_os`: The oxidation state of the metal.
    - `original_complex_charge`: The charge of the complex.
    - `original_complex_id`: The CSD code of the complex.

`occurrences`:
    The number of occurrences of the ligand in the CSD.

`same_graph_denticities`:
























{'warnings': [], 'node_label': 'node_label', 'atomic_props': {'x': [-4.7867999999999995, 0.38070000000000004, -1.7165, 0.20350000000000001, -2.6942, -2.6641, -3.7560000000000002, -4.4297, -3.7766, -2.7913, -2.8202, -1.762, -0.7388999999999999, -0.6641999999999999, -1.2865, 0.3246000000000002, 1.2610999999999999, 1.9398999999999997, 1.1818, 1.8239999999999998], 'y': [0.125, 5.1037, -0.01180000000000092, 1.8072999999999997, -0.8981000000000003, -1.5598, -0.8884000000000007, -1.5281000000000002, 0.08589999999999876, 1.030899999999999, 1.6974999999999998, 0.9795999999999996, 1.9879999999999995, 3.1014999999999997, 3.2309, 4.0098, 3.826699999999999, 4.449299999999999, 2.7010000000000005, 2.554499999999999], 'z': [-3.9460000000000006, -3.4881, -1.2637999999999998, -1.1483000000000008, -1.222900000000001, -0.5707000000000004, -2.1006, -2.0541, -3.024700000000001, -3.0923, -3.7389, -2.1949000000000005, -2.1061999999999994, -2.9162, -3.5952, -2.7133000000000003, -1.7553999999999998, -1.6243999999999996, -0.9848999999999997, -0.32840000000000025], 'atoms': ['F', 'F', 'N', 'N', 'C', 'H', 'C', 'H', 'C', 'C', 'H', 'C', 'C', 'C', 'H', 'C', 'C', 'H', 'C', 'H'], 'original_complex_indices': [1, 2, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26]}, 'global_props': {'CSD_code': 'OZIYON', 'n_atoms': 20, 'n_elements': 4, 'molecular_weight': 192.16880632599998, 'n_C_H_bonds': 6, 'LCS_pred_charge': 0.0, 'LCS_pred_charge_confidence': 1.0, 'LCS_pred_charge_is_confident': True, 'LCS_pred_charge_exact': 0.0}, 'mol': Atoms(symbols='F2N2CHCHC2HC3HC2HCH', pbc=False), 'n_atoms': 20, 'n_hydrogens': 6, 'n_protons': 98, 'graph': <networkx.classes.graph.Graph object at 0x181d01a00>, 'graph_hash': 'd810e651de6b310aabd5ca7060829beb', 'n_bonds': 21, 'has_bond_order_attribute': True, 'has_unknown_bond_orders': False, 'has_good_bond_orders': True, 'heavy_atoms_graph_hash': 'f76078eb3bbe68614cc779c42ff70282', 'bond_order_graph_hash': '59194cf8052a23ab8b2e41804249930e', 'stoichiometry': 'C10H6F2N2', 'original_complex_id': 'OZIYON', 'local_elements': ['N', 'N'], 'was_connected_to_metal': True, 'original_metal': 77, 'original_metal_position': [0.0, 0.0, 0.0], 'original_metal_symbol': 'Ir', 'original_metal_os': 3.0, 'is_centrosymmetric': False, 'centrosymmetry_ang_dev': nan, 'graph_hash_with_metal': '9cfe1644c35cf7f9ef3b747b268cd586', 'heavy_atoms_graph_hash_with_metal': '8d31df32a8d11ecf0b01db06d7cba93f', 'has_betaH': True, 'has_neighboring_coordinating_atoms': False, 'stats': {'min_distance_to_metal': 2.1315960991707597, 'max_distance_to_metal': 6.204836439423686, 'coordinating_atom_distances_to_metal': [2.1315960991707597, 2.150892472905143], 'min_atomic_distance': 0.9291298778965196, 'max_atomic_distance': 8.745613586821682, 'max_dist_per_atoms': 0.43728067934108406}, 'unique_name': 'unq_CSD-OZIYON-02-a', 'pred_charge': 0.0, 'pred_charge_is_confident': True, 'all_ligand_names': ['CSD-OZIYON-02-a'], 'identical_ligand_info': {'name': ['CSD-OZIYON-02-a'], 'original_metal_symbol': ['Ir'], 'original_metal_os': [3.0], 'original_complex_charge': [0], 'original_complex_id': ['OZIYON']}, 'occurrences': 1, 'same_graph_denticities': [2], 'count_metals': {'Ir': 1}, 'n_same_graph_denticities': 1, 'n_metals': 1, 'n_same_graphs': 1, 'has_unconnected_ligands': False, 'all_ligands_metals': ['Ir'], 'same_graph_charges': {'0.0': 1}, 'n_pred_charges': 1, 'common_graph_with_diff_n_hydrogens': False, 'n_electrons': 98.0, 'odd_n_electron_count': False, 'has_warnings': False, 'hash': 246019352813142257304139174990534806063, 'coordinates': {0: ['F', [-4.7867999999999995, 0.125, -3.9460000000000006]], 1: ['F', [0.38070000000000004, 5.1037, -3.4881]], 2: ['N', [-1.7165, -0.01180000000000092, -1.2637999999999998]], 3: ['N', [0.20350000000000001, 1.8072999999999997, -1.1483000000000008]], 4: ['C', [-2.6942, -0.8981000000000003, -1.222900000000001]], 5: ['H', [-2.6641, -1.5598, -0.5707000000000004]], 6: ['C', [-3.7560000000000002, -0.8884000000000007, -2.1006]], 7: ['H', [-4.4297, -1.5281000000000002, -2.0541]], 8: ['C', [-3.7766, 0.08589999999999876, -3.024700000000001]], 9: ['C', [-2.7913, 1.030899999999999, -3.0923]], 10: ['H', [-2.8202, 1.6974999999999998, -3.7389]], 11: ['C', [-1.762, 0.9795999999999996, -2.1949000000000005]], 12: ['C', [-0.7388999999999999, 1.9879999999999995, -2.1061999999999994]], 13: ['C', [-0.6641999999999999, 3.1014999999999997, -2.9162]], 14: ['H', [-1.2865, 3.2309, -3.5952]], 15: ['C', [0.3246000000000002, 4.0098, -2.7133000000000003]], 16: ['C', [1.2610999999999999, 3.826699999999999, -1.7553999999999998]], 17: ['H', [1.9398999999999997, 4.449299999999999, -1.6243999999999996]], 18: ['C', [1.1818, 2.7010000000000005, -0.9848999999999997]], 19: ['H', [1.8239999999999998, 2.554499999999999, -0.32840000000000025]]}, 'atomic_index_to_graph_index': {0: 1, 1: 2, 2: 7, 3: 8, 4: 11, 5: 12, 6: 13, 7: 14, 8: 15, 9: 16, 10: 17, 11: 18, 12: 19, 13: 20, 14: 21, 15: 22, 16: 23, 17: 24, 18: 25, 19: 26}, 'graph_index_to_atomic_index': {1: 0, 2: 1, 7: 2, 8: 3, 11: 4, 12: 5, 13: 6, 14: 7, 15: 8, 16: 9, 17: 10, 18: 11, 19: 12, 20: 13, 21: 14, 22: 15, 23: 16, 24: 17, 25: 18, 26: 19}, 'denticity': 2, 'name': 'CSD-OZIYON-02-a', 'ligand_to_metal': [2, 3], 'csd_code': 'unq_CSD-OZIYON-02-a'}