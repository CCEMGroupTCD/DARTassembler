import os
import subprocess
from pathlib import Path
from typing import Optional, Dict

import yaml

from sigopt_hyperopt.experiments.experiment import BWUniClusterExperiment, HorekaExperiment, ExperimentException, \
    Experiment
from sigopt_hyperopt.utils.configs import Config
from sigopt_hyperopt.utils.const import SIGOPT_TOKEN_FILE_ENV_NAME, JOB_ID_FILE_NAME, LOCAL_EXPERIMENT_FOLDER_NAME, \
    EXPERIMENT_CONFIG_NAME
from sigopt_hyperopt.utils.sigopt_utils import load_sigopt_token
from sigopt_hyperopt.utils.workspace_utils import Workspace


def _read_config(config_path: str) -> Dict:
    assert os.path.isfile(config_path), f"No config file found in {config_path}"
    with open(config_path, 'r') as stream:
        data_loaded = yaml.safe_load(stream)
    return data_loaded


def _create_experiment_object(config: Config, sigopt_api_key: str, init: bool) -> Experiment:
    if config.experiment['cluster'] == "bwunicluster":
        experiment = BWUniClusterExperiment(config=config, sigopt_api_key=sigopt_api_key, init=init)
    elif config.experiment['cluster'] == "horeka":
        experiment = HorekaExperiment(config=config, sigopt_api_key=sigopt_api_key, init=init)
    else:
        raise ExperimentException("config.experiment.cluster must either be 'bwunicluster' or 'horeka'!")
    return experiment


def _check_sigopt_env() -> str:
    assert SIGOPT_TOKEN_FILE_ENV_NAME in os.environ.keys(), "Please specify the sigopt env file as an argument " \
                                                            "for the command or specfiy it in the env variable " \
                                                            f"'{SIGOPT_TOKEN_FILE_ENV_NAME}!'"
    sigopt_env_path = os.environ.get(SIGOPT_TOKEN_FILE_ENV_NAME)
    return sigopt_env_path


def _initialize_experiment_object(config: Config, init: bool = False) -> Experiment:
    sigopt_env_path = _check_sigopt_env()
    sigopt_api_key = load_sigopt_token(env_path=sigopt_env_path, use_dev_token=config.sigopt_options["dev_run"])
    experiment = _create_experiment_object(config=config, sigopt_api_key=sigopt_api_key, init=init)
    return experiment


def _search_in_slurm_workspaces(experiment_name: str) -> Optional[str]:
    try:
        path = subprocess.run(['ws_find', experiment_name], stdout=subprocess.PIPE)
        path = path.stdout.decode('utf-8').split('\n')[0]
        # if no workspace exist, the path is just an empty string
        if not path:
            path = None
    except FileNotFoundError:
        path = None
    return path


def _search_in_local_workspace(experiment_name: str) -> Optional[str]:
    path = None
    local_workspace_base_path = os.path.join(str(Path(__file__).parent.parent.parent.resolve()),
                                             LOCAL_EXPERIMENT_FOLDER_NAME)
    if os.path.exists(local_workspace_base_path):
        local_experiments = [directory for directory in os.listdir(local_workspace_base_path)]
        if experiment_name in local_experiments:
            path = os.path.join(local_workspace_base_path, experiment_name)
    return path


def _get_config_path(experiment_name: str):
    path = _search_in_slurm_workspaces(experiment_name=experiment_name)
    if path is None:
        path = _search_in_local_workspace(experiment_name=experiment_name)
    if path is not None:
        config_path = os.path.join(path, EXPERIMENT_CONFIG_NAME)
        return config_path
    else:
        raise ExperimentException(f"Could not find the experiment folder for the experiment {experiment_name}.")


def start_experiment(config_path: str):
    config = Config()
    config.load(path=config_path)
    experiment = _initialize_experiment_object(config=config, init=True)
    experiment.initialize_experiment()
    experiment.start_experiment()


def continue_experiment(experiment_name: str):
    config = Config()
    config_path = _get_config_path(experiment_name=experiment_name)
    config.load(path=config_path)
    experiment = _initialize_experiment_object(config=config)
    if not experiment.check_initialization():
        raise ExperimentException(f"The experiment {experiment.config.experiment['experiment_name']} is not yet "
                                  f"initialized. Please use the start command to initialize and start the experiment!")
    experiment.start_experiment()


def delete_experiment(experiment_name: str):
    config = Config()
    config_path = _get_config_path(experiment_name=experiment_name)
    config.load(path=config_path)
    experiment = _initialize_experiment_object(config=config)
    experiment.delete_experiment()


def kill_experiment_jobs(experiment_name: str):
    config = Config()
    config_path = _get_config_path(experiment_name=experiment_name)
    config.load(path=config_path)
    experiment = _initialize_experiment_object(config=config)
    experiment.kill_jobs()

