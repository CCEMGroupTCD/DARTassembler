from pathlib import Path
from typing import Union, List
import os

def get_filepaths_from_dirpath(dirpath):
    dirpath = Path(dirpath).resolve()
    name = dirpath.name

    structure_file = Path(dirpath, f'{name}_structure.xyz')
    data_file = Path(dirpath, f'{name}_data.json')
    ligandinfo_file = Path(dirpath, f'{name}_ligandinfo.csv')
    xtbstructure_file = Path(dirpath, f'{name}_xtbstructure.xyz')
    xtbdata_file = Path(dirpath, f'{name}_xtbdata.json')
    gaussiandata_file = Path(dirpath, f'{name}_gaussian.fchk')

    return (structure_file, data_file, ligandinfo_file, xtbstructure_file, xtbdata_file, gaussiandata_file)


def is_complex_directory(dir: Union[str, Path]) -> bool:
    important_files = ['_structure.xyz', '_data.json']
    name = Path(dir).name
    return all([Path(dir, f'{name}{file}').exists() for file in important_files])


def get_all_complex_directories(dir: Union[str, Path]) -> List[str]:
    dir = Path(dir).resolve()
    if not dir.exists():
        raise ValueError(f'Directory "{dir}" does not exist!')

    return [str(dirpath) for dirpath, _, _ in os.walk(dir) if is_complex_directory(dirpath)]