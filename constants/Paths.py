from pathlib import Path
import os
#todo rename path or Dart or something
class project_path:
    def __init__(self):
        # This is the file or folder that marks the origin of the project. It is used to find the origin of the project.
        self.files_marking_origin_of_project = ['.git', 'constants', 'data']
        self.project_path = self.get_project_origin()

    def extend(self, *relpaths: object) -> object:
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


serverpath = Path(project_path().extend('ccem_server'))
