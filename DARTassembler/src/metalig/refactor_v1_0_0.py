from copy import deepcopy

from DARTassembler.src.ligand_extraction.utilities import check_equal

refactor_ligand_dict = {
    "atomic_props": "atomic_props",
    "global_props": "global_props",
    "graph_dict": "graph",
    "ligand_to_metal": "donor_idc",
    "identical_ligand_info": "ligand_instances",
    "original_metal_position": "parent_metal_position",
}
refactor_global_props_dict = {
    "unique_name": "unique_name",
    "pred_charge": "charge",
    "molecular_weight": "molecular_weight",
    "denticity": "n_donors",
    "n_atoms": "n_atoms",
    "n_elements": "n_elements",
    "n_bonds": "n_bonds",
    "n_electrons": "n_electrons",
    "n_protons": "n_protons",
    "occurrences": "n_ligand_instances",
    "has_good_bond_orders": "has_all_bond_orders_valid",
    'pred_charge_is_confident': 'has_confident_charge',
    "original_complex_id": "parent_complex_id",
    "original_metal_symbol": "parent_metal",
    "original_metal_os": "parent_metal_os",
    # graph hashes
    "graph_hash": "graph_hash",
    "graph_hash_with_metal": "graph_hash_with_metal",
    "heavy_atoms_graph_hash": "heavy_atoms_graph_hash",
    "heavy_atoms_graph_hash_with_metal": "heavy_atoms_graph_hash_with_metal",
    "bond_order_graph_hash": "bond_order_graph_hash",
}
refactor_other_ligand_instances_dict = {
    "name": "ligand_name",
    "original_complex_id": "parent_complex_id",
    "original_complex_charge": "parent_complex_charge",
    "original_metal_symbol": "parent_metal",
    "original_metal_os": "parent_metal_os",
}
remove_properties = [
    'all_ligand_names',
    "n_C_H_bonds",
    "n_hydrogens",
    "has_neighboring_coordinating_atoms",
    "has_betaH",
    'stoichiometry',
    "has_bond_order_attribute",
    "has_unknown_bond_orders",
    'n_metals',
    'odd_n_electron_count',
    'is_centrosymmetric',
    'centrosymmetry_ang_dev',
    'original_metal',
    'local_elements',
    'all_ligands_metals',
    'count_metals',
    'was_connected_to_metal',
    'warnings',
    'has_warnings',
    'CSD_code',
    'LCS_pred_charge',
    'LCS_pred_charge_is_confident',
    'LCS_pred_charge_confidence',
    'LCS_pred_charge_exact',
    'same_graph_charges',
    'n_pred_charges',
    'same_graph_denticities',
    'n_same_graph_denticities',
    'n_same_graphs',
    'common_graph_with_diff_n_hydrogens',
    'has_unconnected_ligands',
    # from stats
    'min_distance_to_metal',
    'max_distance_to_metal',
    'coordinating_atom_distances_to_metal',
    'max_dist_per_atoms',
    "min_atomic_distance",
    "max_atomic_distance",
]
make_integer = ['n_electrons', 'pred_charge']

def reindex_graph_dict(graph_dict: dict[dict]) -> dict[dict]:
    """
    Re-index the graph dictionary such that the indices are integers starting from 0.
    :param graph_dict: A dictionary with the graph representation of a ligand, containing two subdictionaries: 'graph' and 'node_attributes'.
    :return: The re-indexed graph dictionary in the same format as the input.
    """

    old_idx_to_new_idx = {old_idx: new_idx for new_idx, old_idx in enumerate(graph_dict['graph'].keys())}

    graph_bonds_new = {}
    node_attributes_new = {}
    for old_idx1, bonds in graph_dict['graph'].items():
        # Update node attributes
        new_idx1 = old_idx_to_new_idx[old_idx1]
        node_attributes_new[new_idx1] = graph_dict['node_attributes'][old_idx1]
        # Update bonds
        bonds = {old_idx_to_new_idx[old_idx2]: bond for old_idx2, bond in bonds.items()}
        graph_bonds_new[new_idx1] = bonds

    graph_dict_new = {
        'graph': graph_bonds_new,
        'node_attributes': node_attributes_new,
    }

    return graph_dict_new

def refactor_metalig_entry_from_v1_0_0_to_v1_1_0(ligand: dict) -> dict:
    """
    Refactor the MetaLig ligand database from version 1.0.0 to 1.1.0. This includes renaming, moving and removing some properties.
    """
    ligand = deepcopy(ligand)

    # Change the order of the properties in identical_ligand_info so that the current ligand is the first entry
    idx = ligand['identical_ligand_info']['name'].index(ligand['name'])
    for prop in ligand['identical_ligand_info']:
        ligand['identical_ligand_info'][prop] = [ligand['identical_ligand_info'][prop][idx]] + [ligand['identical_ligand_info'][prop][i] for i in range(len(ligand['identical_ligand_info'][prop])) if i != idx]
    # Assert ["original_complex_id", "original_metal_symbol", "original_metal_os"] are also in other instances
    assert check_equal(ligand['original_complex_id'], ligand['identical_ligand_info']['original_complex_id'][0])
    assert check_equal(ligand['original_metal_symbol'], ligand['identical_ligand_info']['original_metal_symbol'][0])
    assert check_equal(ligand['original_metal_os'], ligand['identical_ligand_info']['original_metal_os'][0])

    # Expand the properties from stats into global_props
    for prop in ligand['stats']:
        ligand['global_props'][prop] = ligand['stats'][prop]
    del ligand['stats']

    # Get set of all properties for later checking
    all_props = set(ligand.keys()).union(set(ligand['global_props'].keys())).union(
        set(ligand['identical_ligand_info'].keys()))

    # Make some properties integer
    for prop in make_integer:
        if prop in ligand:
            ligand[prop] = int(ligand[prop])
        if prop in ligand['global_props']:
            ligand['global_props'][prop] = int(ligand['global_props'][prop])

    # Remove some properties
    if 'name' in ligand:
        # Remove 'name' attribute from ligand separately because this label is also in identical_ligand_info
        ligand.pop('name')
    for prop in remove_properties:
        # In ligand
        if prop in ligand:
            ligand.pop(prop)
        # In global_props
        if prop in ligand['global_props']:
            ligand['global_props'].pop(prop)

    # Put some properties into global_props
    for prop in refactor_global_props_dict.keys():
        if prop in ligand:
            ligand['global_props'][prop] = ligand.pop(prop)

    # Sort global_props by order of refactor_global_props_dict
    props_before = set(ligand['global_props'].keys())
    ligand['global_props'] = {k: ligand['global_props'][k] for k in refactor_global_props_dict.keys()}
    props_lost = props_before - set(ligand['global_props'].keys())
    assert len(props_lost) == 0, f'Lost properties in global_props while sorting: {props_lost}'

    # Rename properties
    for prop, new_prop in refactor_ligand_dict.items():
        if prop in ligand:
            ligand[new_prop] = ligand.pop(prop)
    # Rename properties in identical_ligand_info
    for prop, new_prop in refactor_other_ligand_instances_dict.items():
        if prop in ligand['ligand_instances']:
            ligand['ligand_instances'][new_prop] = ligand['ligand_instances'].pop(prop)
    # Rename properties in global_props
    for prop, new_prop in refactor_global_props_dict.items():
        if prop in ligand['global_props']:
            ligand['global_props'][new_prop] = ligand['global_props'].pop(prop)

    # Check if all properties are covered
    all_props_after_expected_old = all_props - set(remove_properties)
    # Convert all old prop names to new prop names
    all_props_after_expected_new = []
    for prop in all_props_after_expected_old:
        if prop in refactor_ligand_dict:
            all_props_after_expected_new.append(refactor_ligand_dict[prop])
        elif prop in refactor_global_props_dict:
            all_props_after_expected_new.append(refactor_global_props_dict[prop])
        elif prop in refactor_other_ligand_instances_dict:
            all_props_after_expected_new.append(refactor_other_ligand_instances_dict[prop])
        else:
            raise ValueError(f'Property "{prop}" is not covered by the refactoring.')
    all_props_after_expected_new = set(all_props_after_expected_new)
    all_real_props = set(ligand.keys()).union(set(ligand['global_props'].keys())).union(
        set(ligand['ligand_instances'].keys()))
    assert all_props_after_expected_new == all_real_props, f'Not all properties are covered by the refactoring. Uncovered: {all_props_after_expected_new.symmetric_difference(all_real_props)}'

    # Remove some more properties which would throw an error in the last assertion because they have the same name as properties in the ligand_instances.
    ligand['global_props'].pop('parent_complex_id')
    ligand['global_props'].pop('parent_metal')
    ligand['global_props'].pop('parent_metal_os')
    ligand.pop('parent_metal_position')

    # Rename some special properties
    ligand['atomic_props']['parent_complex_idc'] = ligand['atomic_props'].pop('original_complex_indices')

    # Reindex graph and make indices integers
    ligand['graph'] = reindex_graph_dict(ligand['graph'])

    # Add new properties
    optional_props = {'hapdent_idc': None, 'geometric_isomers_hapdent_idc': None}
    for prop, default in optional_props.items():
        if not prop in ligand:
            ligand[prop] = default
    ligand['ligand_instances']['parent_metal_position'] = [[0., 0., 0.] for _ in range(len(ligand['ligand_instances']['ligand_name']))]

    return ligand