import os.path
import re
import subprocess
from enum import Enum
from typing import List, Optional

import yaml

from sigopt_hyperopt.utils.const import WORKER_CHAIN_BASE_NAME, JOB_IDS_FILE_NAME
from sigopt_hyperopt.utils.logger_utils import setup_logger

logger = setup_logger(__name__)


class TimeFormats(Enum):
    MINUTES = "^\d+$"
    MINUTES_SECONDS = "^\d+:\d+$"
    HOURS_MINUTES_SECONDS = "^\d+:\d+:\d+$"


class TimeFactors(Enum):
    MINUTES = 1
    MINUTES_SECONDS = (1, 1/60)
    HOURS_MINUTES_SECONDS = (60, 1, 1/60)


class SlurmError(Exception):
    """Raised when sth. fails on the slurm domain."""


def parse_slurm_options(slurm_options: dict, as_command=False) -> str:
    if as_command:
        struct = " {prepend}{key}={val}"
    else:
        struct = "#SBATCH {prepend}{key}={val}\n"

    options_str = ""
    partition = None
    for key, val in slurm_options.items():
        if key == 'p' or key == 'partition':
            partition = val
            pass
        prepend = '-' if len(key) == 1 else '--'
        options_str += struct.format(prepend=prepend, key=key, val=val)
    if partition is None:
        raise SlurmError("No partition (`p` or `partition option`) was defined!")
    return options_str


def parse_modules(module_options: list) -> str:
    struct = "module load {module}\n"
    option_str = ""
    for module in module_options:
        option_str += struct.format(module=module)
    return option_str


def get_jobs_nodelist(job_id: int) -> str:
    result = subprocess.run(['squeue', "-j", f"{job_id}", '--Format', "NodeList"],
                            stdout=subprocess.PIPE, encoding='utf-8')
    state = result.stdout.split("\n")[1].strip()
    return state


def _calc_remaining_time(time_string) -> float:
    if bool(re.match(TimeFormats.MINUTES.value, time_string)):
        remaining_time = int(time_string) * TimeFactors.MINUTES.value
    elif bool(re.match(TimeFormats.MINUTES_SECONDS.value, time_string)):
        remaining_time = sum(i*j for i, j in zip(map(int, time_string.split(':')), TimeFactors.MINUTES_SECONDS.value))
    elif bool(re.match(TimeFormats.HOURS_MINUTES_SECONDS.value, time_string)):
        remaining_time = sum(i*j for i, j in zip(map(int, time_string.split(':')),
                                                 TimeFactors.HOURS_MINUTES_SECONDS.value))
    else:
        raise SlurmError(f"Could not calculate remaining time for time string: {time_string}!")
    return remaining_time


def get_jobs_time_left(job_id: int) -> float:
    result = subprocess.run(["squeue", "--job", f"{job_id}", "--Format", "TimeLeft"], stdout=subprocess.PIPE,
                            encoding='utf-8')
    remaining_time_string = result.stdout.split("\n")[1].strip()  # days-hours:minutes:seconds
    if "-" in remaining_time_string:
        days, rest_time = remaining_time_string.split("-")
        days_in_minutes = int(days) * 24 * 60
        remaining_time = _calc_remaining_time(time_string=rest_time) + days_in_minutes

    else:
        remaining_time = _calc_remaining_time(time_string=remaining_time_string)
    return remaining_time


class SlurmWorkerChain:
    """sbatches a script so that is restarted n time as after finishing.

    Args:
        experiment_bash_script_path: (str) Path to the bash script that should be batched.
        partition: (str) Slurm queue name to which the script should be submitted.

    Kwargs:
        n: (int) How many times the script should be restarted.
    """

    def __init__(self, chain_index: int, experiment_bash_script_path: str, partition: str, workspace_path: str, number_chain_jobs: int = 4):
        self.chain_index = chain_index
        self.experiment_bash_script_path = experiment_bash_script_path
        self.workspace_path = workspace_path
        self.number_chain_jobs = number_chain_jobs
        self.partition = partition
        self.job_ids = dict()

    def _save_job_ids(self):
        chain_job_ids_file_path = os.path.join(self.workspace_path, WORKER_CHAIN_BASE_NAME + f"{self.chain_index+1}",
                                               JOB_IDS_FILE_NAME)
        logger.info(f"Saving job_ids in file {chain_job_ids_file_path}")
        with open(chain_job_ids_file_path, "w") as file:
            yaml.safe_dump(self.job_ids, file)

    def _run_chain_job(self, chain_position: int, job_id: Optional[int] = None) -> subprocess.CompletedProcess:
        if chain_position == 0:
            result = subprocess.run(
                ['sbatch',
                 '-p',
                 f'{self.partition}',
                 self.experiment_bash_script_path,
                 str(chain_position)
                 ],
                stdout=subprocess.PIPE, encoding='utf-8', stderr=subprocess.STDOUT
            )
        else:
            result = subprocess.run(
                ['sbatch',
                 '-p',
                 f'{self.partition}',
                 "-d",
                 f"afterok:{job_id}",
                 self.experiment_bash_script_path,
                 str(chain_position)
                 ],
                stdout=subprocess.PIPE, encoding='utf-8', stderr=subprocess.STDOUT
            )
        return result

    @staticmethod
    def _extract_job_id(result_output: str) -> int:
        job_id = int(re.findall(r'\d+', result_output)[0])
        return job_id

    def start(self):
        """Starts the chain.

        Returns:
            job_ids: (list) A list of job ids corresponding to the started jobs.
        """
        job_id = None
        for chain_position in range(self.number_chain_jobs):
            logger.info(f'Starting for the {self.chain_index + 1}. chain worker the chain '
                        f'{chain_position + 1}/{self.number_chain_jobs}')

            if job_id is None:
                result = self._run_chain_job(chain_position=chain_position)
            else:
                result = self._run_chain_job(chain_position=chain_position, job_id=job_id)

            if result.returncode != 0:
                raise SlurmError(f'Failed to start worker chain. \n Error: \n {result.stdout}')

            job_id = SlurmWorkerChain._extract_job_id(result_output=result.stdout)

            logger.info(f'\tStarted worker with job id {job_id}')

            self.job_ids[chain_position] = job_id

        self._save_job_ids()

    def cancel(self, job_ids: List = None):
        """Cancels the chain.

        Kwargs:
            job_ids: (list) A list of job ids to be canceled. Must be provided if the start method of the
                     object was never called.
        """
        assert job_ids or self.job_ids, "This instance has no job_ids, please provide them manually."
        job_ids = job_ids or self.job_ids

        for job_id in job_ids:
            subprocess.run(['scancel', f'{job_id}'])
