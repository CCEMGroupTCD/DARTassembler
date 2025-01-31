from DARTassembler.src.ligand_extraction.DataBase import MoleculeDB


class LigRotNorm:
    """
    This class is responsible for rotating ligands such that they are normal to cartesian coordinate systems
    """
    def __init__(self, ligand) -> None:     # ligand: MoleculeDB (I removed the hint because it was causing an error)
        self.ligand = ligand    # The ligand object we are working on
        self.process_ligand()   # Process the ligand object


    def rotate_monodenate(self):
        """
        Monodentate ligand is
        :return:
        """
        print("Rotating monodentate ligand ...")

    def rotate_bidentate(self):
        print("Rotating bidentate ligand ...")

    def rotate_bidentate_hapt(self):
        print("Rotating bidentate hapt ligand ...")

    def rotate_tridentate_planar(self):
        print("Rotating tridentate planar ligand ...")

    def rotate_tridentate_non_planar(self):
        print("Rotating tridentate non-planar ligand ...")

    def rotate_tetradentate_planar(self):
        print("Rotating tetradentate planar ligand ...")

    def rotate_tetradentate_non_planar(self):
        print("Rotating tetradentate non-planar ligand ...")

    def rotate_pentadentate(self):
        print("Rotating pentadentate ligand ...")

    def rotate_pentadentate_hapt(self):
        print("Rotating pentadentate hapt ligand ...")

    def process_ligand(self):
        """
        This function determines which ligand rotation helper function to call based on the ligand denticity
        :return:
        """
        if self.ligand.denticity == 1:
            self.rotate_monodenate()
        elif self.ligand.denticity == 2 and not self.ligand.has_neighboring_coordinating_atoms:
            self.rotate_bidentate()
        elif self.ligand.denticity == 2 and self.ligand.has_neighboring_coordinating_atoms:
            self.rotate_bidentate_hapt()
        elif self.ligand.denticity == 3 and not self.ligand.has_neighboring_coordinating_atoms:
            self.rotate_tridentate_planar()
        elif self.ligand.denticity == 3 and self.ligand.non_planar: # todo this attribute may not be correct
            self.rotate_tridentate_non_planar()
        elif self.ligand.denticity == 4:
            self.rotate_tetradentate_planar()
        elif self.ligand.denticity == 4 and self.ligand.non_planar:
            self.rotate_tetradentate_non_planar()
        elif self.ligand.denticity == 5:
            self.rotate_pentadentate()
        elif self.ligand.denticity == 5 and self.ligand.has_neighboring_coordinating_atoms:
            self.rotate_pentadentate_hapt()
        else:
            raise ValueError("Ligand denticity not supported")
        pass
