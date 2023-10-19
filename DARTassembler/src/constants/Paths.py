from pathlib import Path
import os
#todo rename path or Dart or something

def split_paths_on_unix_and_windows(relpaths: str) -> tuple:
    return tuple(p for paths in relpaths.split('\\') for p in paths.split('/'))

class project_path:
    def __init__(self):
        # This is the file or folder that marks the origin of the project. It is used to find the origin of the project.
        self.files_marking_origin_of_project = ['.git', 'DARTassembler', 'data']
        self.project_path = self.get_project_origin()

    def extend(self, *relpaths: object) -> Path:
        if len(relpaths) == 1:
            relpaths = split_paths_on_unix_and_windows(relpaths[0])
        return Path(self.project_path, *relpaths)

    def get_project_origin(self) -> Path:
        # Check if the current working directory is the origin of the project
        cwd = Path.cwd()
        if self.check_if_origin_of_project(cwd):
            return cwd

        # Check if any parent directory is the origin of the project
        for parent_path in cwd.parents:
            if self.check_if_origin_of_project(parent_path):
                return parent_path  # Found the origin of the project
        raise FileNotFoundError(f'Could not find the origin of the project. Please make sure that the project is located in a folder that contains the files or folders `{self.files_marking_origin_of_project}`.')

    def check_if_origin_of_project(self, path: Path) -> bool:
        """
        Check if the path is the origin of the project. This is done by checking if the path contains the file or folder `self.file_marking_origin_of_project`.
        :param path: Path to check
        :return: True if the path is the origin of the project, False otherwise
        """
        path = Path(path)
        all_files_on_same_level = set(os.listdir(path))
        origin_files_in_dir = all_files_on_same_level.intersection(self.files_marking_origin_of_project)
        all_files_exist = len(origin_files_in_dir) == len(self.files_marking_origin_of_project)
        return all_files_exist


# Default paths, to use if not user specified
csd_graph_path = project_path().extend(*'data_input/CSD/graphs'.split('/'))
default_ligand_db_path = project_path().extend(*'data/final_db_versions/MetaLigDB_v1.0.0.jsonlines'.split('/'))
element_data_path = project_path().extend(*'DARTassembler/src/constants/element_data_table.csv'.split('/'))
charge_benchmark_dir = project_path().extend('dev', 'test', 'debug', 'databases', 'charge_benchmark')


