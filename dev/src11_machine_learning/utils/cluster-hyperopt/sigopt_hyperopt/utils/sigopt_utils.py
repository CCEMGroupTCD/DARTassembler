import os
import subprocess
from typing import Optional

import sigopt
from sigopt import exception

from sigopt_hyperopt.utils.const import SIGOPT_TOKEN_FILE_ENV_NAME


class SigoptError(Exception):
    """Raised when sth. fails on the sigopt side."""


def load_sigopt_token(env_path: Optional[str] = None, use_dev_token: bool = True) -> str:
    """Loads a client token from a env_path argument.
    Returns:
        token: (str) The token as a string.
    """
    assert os.path.isfile(env_path), SigoptError(f"No Sigopt env file found in {env_path}")
    with open(env_path, 'r') as f:
        token = None
        for line in f.readlines():
            if not use_dev_token and line.startswith("SIGOPT_TOKEN"):
                token = line.split("=")[1]
            elif use_dev_token and line.startswith("SIGOPT_DEV_TOKEN"):
                token = line.split("=")[1]
        assert token is not None, SigoptError(f"A token could not be resolved in the config file {env_path}. Please "
                                              f"make sure that the file provides SIGOPT_TOKEN, SIGOPT_DEV_TOKEN or "
                                              f"both as keys!")
        return token.strip("\n")


def check_experiment_exists(connection: sigopt.interface.Connection, client_id: int, experiment_name: str,
                            project_name: str):
    experiments = connection.clients(client_id).projects(project_name).experiments().fetch()
    for experiment in experiments.iterate_pages():
        if experiment.name == experiment_name:
            raise SigoptError(f"The experiment {experiment_name} already exists in the project {project_name}. "
                              f"Please delete the experiment in sigopt, rename your experiment in your config or use"
                              f"the continue command to continue your experiment.")


def delete_sigopt_experiment(experiment_id: str, connection: sigopt.interface.Connection):
    try:
        connection.experiments(experiment_id).delete()
    except exception.ApiException:
        raise SigoptError(f"Could not delete the sigopt experiment with the id {experiment_id}. It is probably "
                          f"already deleted or not existing.")
