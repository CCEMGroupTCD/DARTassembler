import argparse

from sigopt_hyperopt.utils.const import SIGOPT_TOKEN_FILE_ENV_NAME
from sigopt_hyperopt.utils.entry_points import continue_experiment, delete_experiment, start_experiment, \
    kill_experiment_jobs
from sigopt_hyperopt.utils.logger_utils import setup_logger

logger = setup_logger(__name__)


def _parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers(title="Supported Commands")

    parser_start = subparsers.add_parser("start", help="Initializes and starts the experiment with the given config.")
    parser_start.add_argument("--config_path", type=str,
                              help="The config path for hyperparamter optimization. Must include following aspects: "
                                   "model arguments, git options, sbatch arguments, model parameters, sigopt options, "
                                   "the metrics used in the experiment and experiment options.")
    parser_start.set_defaults(function=start_experiment)

    parser_init = subparsers.add_parser("continue", help="Continues a specific experiment.")
    parser_init.add_argument("--experiment_name", type=str,
                             help="The name of the experiment which should be continued.")
    parser_init.set_defaults(function=continue_experiment)

    parser_delete = subparsers.add_parser("delete", help="Deletes the experiment workspace and the sigopt experiment.")
    parser_delete.add_argument("--experiment_name", type=str,
                               help="The name of the experiment which should be deleted.")
    parser_delete.set_defaults(function=delete_experiment)

    parser_kill_jobs = subparsers.add_parser("kill_jobs", help="Kills all jobs for a specific experiment.")
    parser_kill_jobs.add_argument("--experiment_name", type=str,
                                  help="The name of the experiment for which the jobs should be deleted.")
    parser_kill_jobs.set_defaults(function=kill_experiment_jobs)

    args = parser.parse_args()
    return args


def main():
    args = _parse_arguments()
    func = args.function
    del args.function
    func(**args.__dict__)


if __name__ == "__main__":
    main()
