import argparse
import logging
import os
import subprocess
import sys
import time
from typing import Dict, List, Any, Optional, Union, Tuple, Callable

from sigopt import Connection
import yaml

from sigopt_hyperopt.utils.configs import Config
from sigopt_hyperopt.utils.const import WORKER_CHAIN_BASE_NAME, JOB_IDS_FILE_NAME, MODEL_REPO_NAME, \
    WORKER_FILE_LOGGER_NAME
from sigopt_hyperopt.utils.entry_points import _check_sigopt_env
from sigopt_hyperopt.utils.logger_utils import setup_logger, add_file_handler
from sigopt_hyperopt.utils.sigopt_utils import load_sigopt_token
from sigopt_hyperopt.utils.slurm_utils import get_jobs_time_left

logger = setup_logger(__name__)


class WorkerException(Exception):
    """Raised when sth. fails on the worker side."""
    pass


class Worker:
    """
    This class deals with the worker job that gets deployed by slurm. It handles the communication between the worker
    and Sigopt and thus, also runs the evaluation for each suggestion.
    """
    def __init__(self, config: Config, sigopt_api_key: str, chain_index: int, chain_position: int):
        self.config = config
        self.conn = Connection(client_token=sigopt_api_key)
        self.chain_index = chain_index
        self.chain_position = chain_position
        self.max_runtime = 0
        file_logger_path = os.path.join(self.config.experiment["workspace_path"],
                                        WORKER_CHAIN_BASE_NAME + f"{chain_index+1}",
                                        WORKER_FILE_LOGGER_NAME)
        add_file_handler(logger=logger, file_logger_path=file_logger_path, level=logging.ERROR)

    @property
    def remaining_observations(self):
        experiment = self.conn.experiments(self.config.experiment['experiment_id']).fetch()
        return experiment.observation_budget - experiment.progress.observation_count

    def _extract_module_import_string(self) -> str:
        if 'entry_point' not in self.config.model:
            raise WorkerException(f'The key "entry_point" has to be in the user config. Please make sure that it exists'
                                  f' in the config used for this experiment {self.config.experiment["experiment_name"]}'
                                  f' and retrigger the experiment!')
        entry_point = self.config.model['entry_point']
        if entry_point.endswith('.py'):
            entry_point = os.path.splitext(entry_point)[0]
        import_string = entry_point.replace('/', '.')
        return import_string

    def _extract_train_function(self, import_string: str) -> Callable:
        if 'function_name' not in self.config.model:
            raise WorkerException(f'The key "function_name" has to be in the user config. Please make sure that it'
                                  f' exists in the config used for this experiment'
                                  f' {self.config.experiment["experiment_name"]} and retrigger the experiment!')

        module = __import__(import_string)
        return getattr(module, self.config.model['function_name'])

    def _create_train_config(self, suggestion) -> Dict:
        train_config = dict()
        train_config["suggestion"] = suggestion.assignments
        train_config["dataset_path"] = self.config.model["dataset_path"]
        train_config["output_path"] = os.path.join(self.config.model["output_path"], suggestion.id)
        return train_config

    def _evaluate_model(self, suggestion) -> Tuple[Union[List, float], Optional[Dict]]:
        value, metadata = None, None
        train_config = self._create_train_config(suggestion=suggestion)
        import_string = self._extract_module_import_string()
        train_function = self._extract_train_function(import_string=import_string)
        os.mkdir(train_config["output_path"])
        try:
            start = time.time()
            value, metadata = train_function(config=train_config)
            end = time.time()
            elapsed_time_doubled = ((end - start) * 2) / 60  # in minutes
            if elapsed_time_doubled > self.max_runtime:
                self.max_runtime = elapsed_time_doubled
                logger.info(f"New max runtime of the model is: {elapsed_time_doubled}")
        except Exception as error:
            logger.error(f"An error was thrown during training:\n {error}")

        return value, metadata

    def _get_current_job_id(self) -> int:
        chain_job_ids_file_path = os.path.join(self.config.experiment["workspace_path"],
                                               WORKER_CHAIN_BASE_NAME + f"{self.chain_index + 1}", JOB_IDS_FILE_NAME)
        with open(chain_job_ids_file_path, "r") as file:
            jobs_dict = yaml.safe_load(file)
        return jobs_dict[self.chain_position]

    def _check_for_time_limit(self) -> bool:
        continue_worker = True
        job_id = self._get_current_job_id()
        remaining_time = get_jobs_time_left(job_id=job_id)
        logger.info(f"Remaining time left: {remaining_time}")
        if remaining_time < self.max_runtime:
            continue_worker = False
        return continue_worker

    def _report_result_to_sigopt(self, value: List, metadata: Optional[Dict], suggestion_id: int, status: bool):
        try:
            if self.config.experiment["multimetric_experiment"]:
                self.conn.experiments(self.config.experiment['experiment_id']).observations().create(
                    suggestion=suggestion_id,
                    values=value,
                    failed=status,
                    metadata=metadata
                )
            else:
                self.conn.experiments(self.config.experiment['experiment_id']).observations().create(
                    suggestion=suggestion_id,
                    value=value,
                    failed=status,
                    metadata=metadata
                )
        except Exception as error:
            logger.error(f"An error occurred during reporting the result to the sigopt server:\n {error}")
            raise WorkerException(f"An error occurred during reporting the result to the sigopt server:\n {error}")

    def run_optimization(self):
        while self.remaining_observations > 0:
            if not self._check_for_time_limit():
                logger.info("No time left for a new model evaluation. Safe exist!")
                sys.exit(0)
            suggestion = self.conn.experiments(self.config.experiment['experiment_id']).suggestions().create()
            value, metadata = self._evaluate_model(suggestion=suggestion)
            if value is not None:
                failed = False
            else:
                failed = True
            self._report_result_to_sigopt(value=value, metadata=metadata, suggestion_id=suggestion.id, status=failed)

    def change_working_directory(self):
        model_path = os.path.join(self.config.experiment["workspace_path"], MODEL_REPO_NAME)
        os.chdir(model_path)

    def check_conda_env(self):
        try:
            activated_conda_env = os.environ['CONDA_DEFAULT_ENV']
        except KeyError:
            raise WorkerException(f'Something went wrong. No conda environment is activated! '
                                  f'Please restart the script.')
        if activated_conda_env != self.config.experiment['experiment_name']:
            raise WorkerException(f'Something went wrong. The wrong conda environment is activated! '
                                  f'The activated environment is {activated_conda_env} but not '
                                  f'{self.config.experiment["experiment_name"]}.')


def _parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config_path", type=str,
                        help="The path to the experiment config.")
    parser.add_argument("--chain_index", type=int,
                        help="The chain index to which the worker belongs.")
    parser.add_argument("--chain_position", type=int,
                        help="THe position of the worker in the chain starting from 0 for the first worker.")
    args = parser.parse_args()
    return args


def _add_repository_to_sys_path():
    sys.path.append("")


def main():
    args = _parse_arguments()
    config = Config()
    config.load(path=args.config_path)
    sigopt_env_path = _check_sigopt_env()
    sigopt_api_key = load_sigopt_token(env_path=sigopt_env_path,
                                       use_dev_token=config.sigopt_options["dev_run"])
    worker = Worker(config=config,
                    sigopt_api_key=sigopt_api_key,
                    chain_index=args.chain_index,
                    chain_position=args.chain_position
                    )
    worker.check_conda_env()
    worker.change_working_directory()
    _add_repository_to_sys_path()
    worker.run_optimization()


if __name__ == '__main__':
    main()
