import json
from pathlib import Path
import os
from typing import Union, Tuple, List

from DARTassembler.src.assembly.TransitionMetalComplex import TransitionMetalComplex
from DARTassembler.src.ligand_extraction.io_custom import load_json


class AssembledComplexOutput(object):

    def __init__(self, dirpath: Union[str, Path]):

        if not self.is_complex_directory(dirpath):
            raise ValueError(f'"{dirpath}" is not a valid complex directory!')

        self.dirpath = Path(dirpath).resolve()
        self.name = self.dirpath.name

        self.dartstructure_file, self.data_file, self.ligandinfo_file, self.xtbstructure_file, self.xtbdata_file = self.get_filepaths()
        self.dartcomplex = TransitionMetalComplex.from_json(self.data_file)

        self.has_xtb = self.xtbstructure_file.exists()
        self.xtbcomplex = TransitionMetalComplex.from_json(self.data_file, xyz=self.xtbstructure_file) if self.has_xtb else None
        # self.xtbdata = load_json(self.xtbdata_file) if self.has_xtb else None




    def get_filepaths(self):
        return self.get_filepaths_from_dirpath(self.dirpath)

    @staticmethod
    def get_filepaths_from_dirpath(dirpath):
        dirpath = Path(dirpath).resolve()
        name = dirpath.name

        structure_file = Path(dirpath, f'{name}_structure.xyz')
        data_file = Path(dirpath, f'{name}_data.json')
        ligandinfo_file = Path(dirpath, f'{name}_ligandinfo.csv')
        xtbstructure_file = Path(dirpath, f'{name}_xtbstructure.xyz')
        xtbdata_file = Path(dirpath, f'{name}_xtbdata.json')

        return (structure_file, data_file, ligandinfo_file, xtbstructure_file, xtbdata_file)

    @staticmethod
    def is_complex_directory(dir: Union[str, Path]) -> bool:
        important_files = ['_structure.xyz', '_data.json']
        name = Path(dir).name
        return all([Path(dir, f'{name}{file}').exists() for file in important_files])

    @staticmethod
    def get_all_complex_directories(dir: Union[str, Path]) -> List[str]:
        dir = Path(dir).resolve()
        if not dir.exists():
            raise ValueError(f'Directory "{dir}" does not exist!')

        return [str(dirpath) for dirpath, _, _ in os.walk(dir) if AssembledComplexOutput.is_complex_directory(dirpath)]



if __name__ == '__main__':
    pass

    # Copy data once
    # import shutil
    # all_xtb_dir = '/Users/timosommer/PhD/projects/RCA/projects/DART/examples/Pd_Ni_Cross_Coupling/dev/xtb_calculations/231013_relaxations/complexes'
    # all_data_dir = '/Users/timosommer/PhD/projects/RCA/projects/DART/examples/Pd_Ni_Cross_Coupling/dev/xtb_calculations/231013_relaxations_and_original_dirs'
    # # copy data
    # dirs = AssembledComplexOutput.get_all_complex_directories(all_data_dir)
    # for complex_dir in dirs:
    #     name = Path(complex_dir).name
    #     xtb_dir = Path(all_xtb_dir, name)
    #     # Copy files from xtb dir to complex dir
    #     for file in os.listdir(xtb_dir):
    #         src = Path(xtb_dir, file)
    #         dst = Path(complex_dir, file)
    #         assert src.exists(), f'File "{src}" does not exist!'
    #         assert not dst.exists(), f'File "{dst}" already exists!'
    #         print(f'From: {src}, To: {dst}')
    #         shutil.copy(src, dst)
