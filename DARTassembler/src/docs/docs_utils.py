from DARTassembler.src.assembly.TransitionMetalComplex import TransitionMetalComplex
from DARTassembler.src.constants.Paths import docs_run_dirs
from DARTassembler.src.constants.Periodic_Table import DART_Element
import moldoc.molecule as molecule
from pathlib import Path

def get_moldoc_molecule(run: str, complex_name: str) -> molecule.Molecule:
    """
    Get a moldoc molecule object from a complex json file.
    :param run: The run name, e.g. 'quickstart'.
    :param complex_name: The complex name, e.g. 'YOZEPECO'.
    :return: A moldoc molecule object.
    """
    rundir = docs_run_dirs[run]
    jsonpath = Path(rundir, complex_name, f'{complex_name}_data.json')
    tmc = TransitionMetalComplex.from_json(jsonpath)
    atm = tmc.atomic_props

    # Extract atoms and bonds from the TransitionMetalComplex object
    atoms = []
    for atom, x, y, z in zip(atm['atoms'], atm['x'], atm['y'], atm['z']):
        atomic_number = DART_Element(atom).atomic_number
        position = (x, y, z)
        atoms.append(molecule.Atom(atomic_number=atomic_number, position=position))
    bonds = [molecule.Bond(atom1_id=id1, atom2_id=id2, order=1,) for id1, id2 in tmc.graph.edges]

    # Create a moldoc molecule object with the specified atoms and bonds
    moldoc_display_molecule = molecule.Molecule(atoms=atoms,bonds=bonds)

    return moldoc_display_molecule


if __name__ == '__main__':
    run = 'quickstart'
    complex_name = 'YOZEPECO'
    moldoc_display_molecule = get_moldoc_molecule(run, complex_name)
