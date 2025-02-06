
import ase
import itertools
AssemblyLigand = 2
ligand_isomers = 2

@classmethod
def from_ligand_combination_return_isomers(cls,
                                           ligands: list[AssemblyLigand],
                                           geometry: str,
                                           metal: str
                                           ):
    # Get all isomers of ligands...
    # Get all complex isomers
    isomers = []
    ligand_isomer_combinations = list(itertools.combinations(ligand_isomers, len(ligand_isomers)))
    for ligands in ligand_isomer_combinations:
        isomer = cls.from_ligands(ligands=ligands, metal=metal)
        isomers.append(isomer)
    # Do simple 'isomer search' with ff-sp to find the best orientation of the monodentate ligand and maybe others. Also do other refactoring for tridentate, trans-bidentate etc. ...
    # ...
    return isomers
@classmethod
def from_ligands(cls,
                  ligands: list[AssemblyLigand],
                  metal: str
                  ):
    complex = ase.Atoms()
    complex.append(ase.Atom(symbol=metal, position=[0, 0, 0]))
    for ligand in ligands:
        complex.extend(ligand)
    return cls(atoms=complex, ligands=ligands)