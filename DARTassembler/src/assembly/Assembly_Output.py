import shutil
from pathlib import Path
from typing import Union
import json
import pandas as pd
from DARTassembler.src.assembly.Assembly_Input import AssemblyInput

_gbl_optimization_movie = 'ffmovie.xyz'
_gbl_concatenated_xyz = 'INTEGRATION_TEST.xyz'
_gbl_run_info_table = 'info_table.csv'
_gbl_batch_dir = 'batches'

# Batch output files
_batch_passed_ff_movie = 'concat_passed_ffmovie.xyz'     # All xyz movies of the forcefield optimization of passed complexes
_batch_failed_ff_movie = 'concat_failed_ffmovie.xyz'     # All xyz movies of the forcefield optimization of failed complexes
_batch_passed_xyz = 'concat_passed_complexes.xyz'        # All xyz files of passed complexes
_batch_failed_xyz = 'concat_failed_complexes.xyz'        # All xyz files of failed complexes
_batch_output = 'batch_output.txt'                       # The batch output file with all stdout output
_batch_errors = 'batch_errors.txt'                       # The batch errors file with all stderr output
_batch_complex_dir = 'complexes'                         # The directory where all the complex output files are stored

# Complex output files
_complex_data = 'data.json'
_complex_structure = 'structure.xyz'
_complex_gaussian = 'gaussian.com'
_complex_info = 'info.txt'
_complex_settings = 'settings.yml'
_complex_ligandfilters = 'ligandfilters.yml'
_complex_warnings = 'warnings.txt'
_complex_ff_movie = 'ffmovie.xyz'
_complex_ligand_info = 'ligandinfo.csv'


def save_file(string: str, outpath: Union[str,Path]):
    with open(outpath, 'w') as f:
        f.write(string)

def append_file(string: str, outpath: Union[str,Path]):
    with open(outpath, 'a') as f:
        f.write(string)

def append_global_optimization_movie(xyz_string, outdir: [str, Path]):
    outpath = Path(outdir, _gbl_optimization_movie)
    append_file(xyz_string, outpath)

def append_global_concatenated_xyz(xyz_string, outdir: [str, Path]):
    outpath = Path(outdir, _gbl_concatenated_xyz)
    append_file(xyz_string, outpath)

def save_batch_optimization_movie(xyz_string, outdir: [str,Path]):
    outpath = Path(outdir, _gbl_optimization_movie)
    append_file(xyz_string, outpath)


class AssemblyOutput(object):

    def __init__(self, outdir: [str,Path], ensure_empty_output_dir: bool = False):
        self.outdir = outdir
        self.run_info_table = Path(self.outdir, _gbl_run_info_table)
        self.batch_dir = Path(self.outdir, _gbl_batch_dir)
        if ensure_empty_output_dir:
            self.ensure_output_directory_empty()

    def save_global_optimization_movie(self, xyz_string):
        append_global_optimization_movie(xyz_string, self.outdir)

    def save_global_concatenated_xyz(self, xyz_string):
        append_global_concatenated_xyz(xyz_string, self.outdir)

    def save_run_info_table(self, df_info: pd.DataFrame):
        df_info.to_csv(self.run_info_table, index=False)

    def ensure_output_directory_empty(self) -> None:
        """
        Checks if the output directory is valid.
        """
        # Delete and recreate directory
        if self.outdir.is_dir():
            shutil.rmtree(self.outdir)
        self.outdir.mkdir(parents=True)

        return



class BatchAssemblyOutput(object):

    def __init__(self, batchdir: Union[str,Path]):
        self.batchdir = Path(batchdir)
        self.batchdir.mkdir(parents=True, exist_ok=True)

        self.passed_ff_movie_path = Path(self.batchdir, _batch_passed_ff_movie)   # All xyz movies of the forcefield optimization of passed complexes
        self.failed_ff_movie_path = Path(self.batchdir, _batch_failed_ff_movie)   # All xyz movies of the forcefield optimization of failed complexes
        self.passed_xyz_path = Path(self.batchdir, _batch_passed_xyz)             # All xyz files of passed complexes
        self.failed_xyz_path = Path(self.batchdir, _batch_failed_xyz)             # All xyz files of failed complexes
        self.output_path = Path(self.batchdir, _batch_output)                     # The batch output file with all stdout output
        self.errors_path = Path(self.batchdir, _batch_errors)                     # The batch errors file with all stderr output
        self.complex_dir = Path(self.batchdir, _batch_complex_dir)                # The directory where all the complex output files are stored

    def save_passed_ff_movie(self, xyz_string: str, append=False) -> None:
        """
        Saves a movie of the successful forcefield optimizations as a concatenated .xyz file.
        @param xyz_string: The concatenated .xyz string.
        """
        self.save_file(xyz_string, self.passed_ff_movie_path, append=append)

    def save_failed_ff_movie(self, xyz_string: str, append=False) -> None:
        """
        Saves a movie of the unsuccessful forcefield optimizations as a concatenated .xyz file.
        @param xyz_string: The concatenated .xyz string.
        """
        self.save_file(xyz_string, self.failed_ff_movie_path, append=append)

    def save_passed_xyz(self, xyz_string: str, append=False) -> None:
        """
        Saves the successful complexes as a concatenated .xyz file.
        @param xyz_string: The concatenated .xyz string.
        """
        self.save_file(xyz_string, self.passed_xyz_path, append=append)

    def save_failed_xyz(self, xyz_string: str, append=False) -> None:
        """
        Saves the unsuccessful/failed complexes as a concatenated .xyz file.
        @param xyz_string: The concatenated .xyz string.
        """
        self.save_file(xyz_string, self.failed_xyz_path, append=append)

    def save_output(self, string: str, append=False) -> None:
        """
        Appends the given string to the batch output file.
        @param string: The string to save.
        """
        self.save_file(string, self.output_path, append=append)

    def save_errors(self, string: str, append=False) -> None:
        """
        Appends the given string to the batch errors file.
        @param string: The string to save.
        """
        self.save_file(string, self.errors_path, append=append)

    def save_file(self, string: str, outpath: Union[str,Path], append=False) -> None:
        """
        Saves a file to the batch directory.
        """
        if append:
            append_file(string, outpath)
        else:
            save_file(string, outpath)

class ComplexAssemblyOutput(object):

    def __init__(self,
                 complexdir: [str,Path],
                 name: Union[str,None] = None,
                 add_name: bool = True
                 ):
        self.complexdir = Path(complexdir)
        self.complexdir.mkdir(parents=True, exist_ok=True)
        self.name = name or self.complexdir.name

        start = f'{self.name}_' if add_name else ''
        self.data_path = Path(complexdir, start + _complex_data)
        self.structure_path = Path(complexdir, start + _complex_structure)
        self.gaussian_path = Path(complexdir, start + _complex_gaussian)
        self.info_path = Path(complexdir, start + _complex_info)
        self.settings_path = Path(complexdir, start + _complex_settings)
        self.ligandfilters_path = Path(complexdir, start + _complex_ligandfilters)
        self.warnings_path = Path(complexdir, start + _complex_warnings)
        self.ff_movie_path = Path(complexdir, start + _complex_ff_movie)
        self.ligand_output_path = Path(complexdir, start + _complex_ligand_info)

    def save_all_complex_data(self,
                              complex,
                              complex_idx: int,
                              xyz_structure: str,
                              ff_movie: Union[str, None] = None,
                              assembly_input_path: [str, Path, None] = None,
                              batch_idx: Union[int, None] = None,
                              ligands: Union[dict, None] = None,
                              ) -> None:

        self.save_structure(xyz_structure)  # Save the structure as xyz
        if ff_movie is not None:
            self.save_ff_movie(ff_movie)  # Save the force field movie
        if assembly_input_path is not None:
            self.save_settings(assembly_input_path)  # Save the assembly input settings
        if ligands is not None:
            self.save_ligand_info(ligands)  # Save the ligand info
        self.save_data_json(
                            complex=complex,
                            complex_idx=complex_idx,
                            xyz_structure=xyz_structure,
                            ff_movie=ff_movie,
                            assembly_input_path=assembly_input_path,
                            batch_idx=batch_idx,
                            ligands=ligands
                            )

    def get_ligand_info_dict(self, ligands: dict) -> dict:
        return {f'Ligand {i}': lig.get_ligand_output_info() for i, lig in enumerate(ligands.values())}
    def save_ligand_info(self, ligands: dict) -> None:
        ligand_infos = self.get_ligand_info_dict(ligands)
        df = pd.DataFrame.from_dict(ligand_infos, orient='index')
        df.to_csv(self.ligand_output_path)

    def save_data_json(self,
                        complex,
                        complex_idx: int,
                        xyz_structure: str,
                        ff_movie: Union[str,None] = None,
                        assembly_input_path: [str,Path,None] = None,
                        batch_idx: Union[int,None] = None,
                        ligands: Union[dict,None] = None
                        ) -> None:
        """
        Saves all data in a contained json file.
        """
        ligand_info = self.get_ligand_info_dict(ligands) if ligands is not None else None
        data = {
                'complex': complex.to_data_dict(),
                'complex_idx': complex_idx,
                'xyz_structure': xyz_structure,
                'assemby_input_path': str(assembly_input_path),
                'assembly_input_settings': AssemblyInput.get_settings_from_input_file(assembly_input_path),
                'batch_idx': batch_idx,
                'ff_movie': ff_movie,
                'ligand_info': ligand_info,
                }
        self.save_file(data, self.data_path)

    def save_structure(self, xyz_string: str) -> None:
        """
        Saves the structure file.
        @param xyz_string: The .xyz string of the molecule
        """
        self.save_file(xyz_string, self.structure_path)

    def save_ff_movie(self, xyz_string: str) -> None:
        """
        Saves the forcefield optimization movie.
        @param xyz_string: The concatenated .xyz string.
        """
        self.save_file(xyz_string, self.ff_movie_path)

    def save_settings(self, assembly_input_filepath: Union[str,Path]) -> None:
        shutil.copy(str(assembly_input_filepath), str(self.settings_path))

    def save_file(self, file: [str,dict], outpath: Union[str,Path], append=False) -> None:
        """
        Saves a file to the complex directory.
        """
        if isinstance(file, (dict,list)):
            file = json.dumps(file, indent=4)

        if append:
            append_file(file, outpath)
        else:
            save_file(file, outpath)



