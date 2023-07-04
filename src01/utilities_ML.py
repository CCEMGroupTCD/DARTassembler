from pymatgen.core.periodic_table import Element as Pymatgen_Element


def get_element_descriptors(el: str) -> dict:
    el = Pymatgen_Element(el)
    descriptors = {
        'Z': el.Z,
        'row': el.row,
        'group': el.group,
        'atomic_mass': el.atomic_mass,
        'atomic_radius': el.atomic_radius,
        'electron_affinity': el.electron_affinity,
        'electronegativity': el.X,
        'min_oxidation_state': el.min_oxidation_state,
        'max_oxidation_state': el.max_oxidation_state,
        'ionization_energy': el.ionization_energy
            }

    return descriptors

