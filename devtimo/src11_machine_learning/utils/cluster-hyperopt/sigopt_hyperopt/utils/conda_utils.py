import logging
import subprocess
from re import search

from sigopt_hyperopt.utils.logger_utils import setup_logger

logger = setup_logger(__name__)


def check_conda_env_exists(env_name: str) -> bool:
    result = subprocess.run(['conda', 'env', 'list'], stdout=subprocess.PIPE)
    conda_envs = result.stdout.decode('utf-8').split('\n')
    for env in conda_envs:
        regex_result = search(rf'^{env_name}', env)
        if regex_result:
            return True
    return False


def remove_conda_env(env_name: str):
    subprocess.run(['conda', 'env', 'remove', '--name', env_name], stdout=subprocess.PIPE)



