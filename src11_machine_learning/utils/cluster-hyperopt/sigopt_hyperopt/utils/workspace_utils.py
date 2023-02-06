import logging
import os.path
import os
import shutil
import subprocess
from pathlib import Path
from typing import Optional

from sigopt_hyperopt.utils.const import LOCAL_EXPERIMENT_FOLDER_NAME
from sigopt_hyperopt.utils.logger_utils import setup_logger

logger = setup_logger(__name__)


class WorkspaceError(Exception):
    """Raised when sth. fails on the workspace domain."""


class Workspace:
    """This class includes logic for creating and managing workspace for each experiment"""
    def __init__(self, experiment_name: str, use_local_workspace: False, duration: int = 30):
        self.experiment_name = experiment_name
        self.use_local_workspace = use_local_workspace
        self.duration = duration
        self.local_workspace_base_path: str = os.path.join(str(Path(__file__).parent.parent.parent.resolve()),
                                                           LOCAL_EXPERIMENT_FOLDER_NAME)

    @property
    def workspace_path(self) -> str:
        if self.use_local_workspace:
            path = os.path.join(self.local_workspace_base_path, self.experiment_name)
        else:
            path = subprocess.run(['ws_find', self.experiment_name], stdout=subprocess.PIPE)
            path = path.stdout.decode('utf-8').split('\n')[0]
        return path

    def create_workspace(self):
        if not self.check_workspace_exist():
            logger.info(f"Create workspace for experiment {self.experiment_name}!")
            if self.use_local_workspace:
                path = os.path.join(self.local_workspace_base_path, self.experiment_name)
                try:
                    os.mkdir(path)
                except OSError as error:
                    raise WorkspaceError(f"Could not create experiment folder in {path}, see stacktrace: {error}")
            else:
                result = subprocess.run(
                    ['ws_allocate',
                     f'{self.experiment_name}',
                     f'{self.duration}'
                     ],
                    stdout=subprocess.PIPE
                )
                path = result.stdout.decode('utf-8').split('\n')[0]
            logger.info(f"Created workspace in path {path}")
        else:
            raise WorkspaceError(f"There already exist a workspace for the experiment {self.experiment_name} in "
                                 f"{self.workspace_path}! Please delete the workspace, rename your experiment in your "
                                 f"config or use the continue command to continue your experiment.")

    def check_workspace_exist(self) -> bool:
        ws_exists = os.path.isdir(self.workspace_path)
        return ws_exists

    def delete_workspace(self):
        if self.use_local_workspace:
            try:
                shutil.rmtree(self.workspace_path)
            except FileNotFoundError:
                raise WorkspaceError(f"The workspace for the experiment {self.experiment_name} "
                                     f"could not be automatically deleted.")
        else:
            command = f"ws_release {self.experiment_name}"
            try:
                subprocess.run(command, stdout=subprocess.PIPE, encoding='utf-8', shell=True, check=True)
            except subprocess.CalledProcessError:
                raise WorkspaceError(f"The workspace for the experiment {self.experiment_name} could not be "
                                     f"automatically deleted wth command: {command}")
