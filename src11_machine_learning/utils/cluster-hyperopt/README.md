# Sigopt-Hyperopt
This repository can be used to tune the hyperparameters of a model. This is achieved with help of 
[sigopt](https://sigopt.com/).

## Introduction
Sigopt-Hyperopt mainly contains of three parts for hyperparameter optimization:
- Sigopt
- Cluster Environment
- Experiment Config File

This repository automatically connects all three aspects in order to give an easy-to-use framework to optimize 
hyperparameters.

### Sigopt
[Sigopt](https://sigopt.com/) is a service provider for optimizing hyperparameters for a given model. The starting point
for any optimization is an experiment that needs to be created via the API or the UI. In general, sigopt can be used 
via the API or UI, though the API offers a more flexible and complex utilization of the sigopt features.

This experiment defines important aspects as which hyperparameters should be optimized, the metrics that should 
be used for the evaluation, or the total number of observations in the optimization.

The optimization itself is a loop that starts with Sigopt giving a suggestion for the hyperparameters of a specific 
model. An evaluation of the suggestions will be based on the performance of the corresponding model with regard to one 
or metrics. Next, the evaluated suggestions - called observations - are reported back to Sigopt, and the loop starts 
again. 

One can also divide the optimization loop in two areas of responsibilities. The first one is only
responsible for creating the suggestions which is covered by the Sigopt service. The user itself has to deal with 
evaluating the suggestions and reporting it back to sigopt.

Please read the official [documentation](https://app.sigopt.com/docs) for a more detailed overview of sigopt.

Additionally, sigopt offers an [academic plan](https://sigopt.com/solution/academia/) that provides additional features
for the user.

### Cluster Environment
All the work that needs to be dealt with by the user with regard to the optimization loop described above, is carried 
out in a cluster environment like [BW-Unicluster](https://wiki.bwhpc.de/e/Category:BwUniCluster_2.0) or 
[HoreKA](https://www.scc.kit.edu/dienste/horeka.php) - as for now, only the BW-Unicluster is supported. 
This includes the evaluation of the sigopt suggestion and the communication between sigopt and the main process.

Both HoreKA and BW-Unicluster uses [slurm](https://slurm.schedmd.com/documentation.html) for their job scheduling. 
Please refer to the official documentation for more information.

Additionally, BW-Unicluster supports different 
[software modules](https://wiki.bwhpc.de/e/BwUniCluster_2.0_Software#Display_all_available_Modules) 
that can be loaded and used in the experiments like cuda support.

The cluster environment ist also used to create a [workspace](https://wiki.bwhpc.de/e/Workspace) for 
each experiment if the config `use_local_workspace` is set to false.

### Experiment Config
The experiment config file is the main file each user needs to define for their experiments. It covers both previous 
aspects and contains all necessary information for the whole optimization process. It is divided in following sections:
- model specific information (entry point, function name, ...)
- git options for the model (git path, branch name, ...)
- experiment information that is used to create the sigopt experiment
- cluster modules (see [modules docu](https://wiki.bwhpc.de/e/BwUniCluster_2.0_Software#Display_all_available_Modules))
- sbatch options (see [sbatch options](https://wiki.bwhpc.de/e/BwUniCluster_2.0_Slurm_common_Features#Slurm_Commands_.28excerpt.29))
- hyperparameters information (needs to fulfill requirements in the [sigopt docu](https://app.sigopt.com/docs/archive/overview/create))
- metrics information (needs to fulfill requirements in the [sigopt docu](https://app.sigopt.com/docs/archive/overview/create))
- sigopt options (experiment name, client id, ...)

The example config file  `test_conf.yaml` in the directory `demo` gives a working example of all possible fields 
including an explanation for the fields.

## User Guide

The following section will give you an overview on how to use Sigopt-Hyperopt based on an example use case that tries to 
optimize the hyperparameters alpha and beta of the [dummy model](https://github.com/aimat-lab/elasticnet). 

### Preliminaries
In order to use Sigopt-Hyperopt you need access to:
- Sigopt
- Cluster environment
- Conda environment
- LSDF (optional)

After signing up to Sigopt, one can read the [API tokens](https://app.sigopt.com/tokens/info) needed for the 
communication between Sigopt and Sigopt-Hyperopt. These keys should to be copied to a file in the respective cluster:
```buildoutcfg
SIGOPT_TOKEN=your-token-belongs-here
SIGOPT_DEV_TOKEN=your-token-belons-here
```
The next step is to assign an environment variable called `SIGOPT_ENV_FILE` to the path of this newly created file.
The difference between these two tokens is that the `SIGOPT_DEV_TOKEN` only provides a random suggestion for developing 
purposes whereas `SIGOPT_TOKEN` provides serious suggestions. This option can then be activated or deactivated in the 
config file via `sigopt_options.dev_run`.

An additional step needs to be taken for your conda installation in the cluster environment which is described in this 
[wiki page](https://wiki.bwhpc.de/e/Conda). Note that it is already sufficient if you just install it in your root 
directory and not in a workspace.

In order to use LSDF, one can follow the tutorial given in the [aimat wiki](https://aimat.iti.kit.edu/wiki/index.php?title=LSDF).

### Custom Model
After doing all the preliminaries, we can start to optimize the hyperparameters of our dummy model. The repository of
the dummy model consists of three parts: 
- a config (not needed in our use case)
- an environment.yaml
- a train script 

The `environment.yaml` includes the information about the required packages for the model itself and has to be 
presented - with this exact name - in each model repository you want to optimize. This file will then be used in the 
workflow to create a conda environment.

The train script includes the desired function for the entry point of the suggestion evaluation. Sigopt-Hyperopt will 
always call this function with a config dictionary that includes the dataset path, the suggestions and the output path 
(can be used to save checkpoints of the model). Therefore, you need to process these three arguments in your script.

In the case of the dummy model, the train script will use the suggestion to train a model with cross validation and 
measures the performance with help of the mean squared error metric. At the end of the run, it returns the evaluation 
of the suggestion with help of the selected metrics.

This setup should also be present in your model repository similarly. The entry script, and the corresponding function 
can be specified in the config via `model.entry_point` and `model.function_name`.

### Run the Optimization
In order to tune the hyperparameters of you custom model the following steps need to be completed:
1. Login to the cluster environment
2. Activate a conda environment where Sigopt-Hyperopt should be installed. If the conda cli is not working even though
it was properly installed, one can run `source ~/.bashrc` or `source ~/.zshrc` depending on the used shell
3. Configure your experiment config file. Using the demo config, it should be sufficient to just change the arguments of 
    `sigopt_options`
4. Run the `python hyperopt.py start --config_path=path-to-your-config-file`

Thats it!

### CLI Commands

Sigopt-Hyperopt supports four commands that can also be viewed with `python hyperopt.py --help`.
- start: Initializes (creates sigopt experience, workspace, conda env, ...) the experiment and starts 
  the hyperparameter optimization. Argument is the `config_path` which refers to the experiment config file
- continue: The experiment will be continued with help of the experiment name via the argument `experiment_name`
- delete: The workspace and sigopt experiment will be deleted for a given experiment. 
  The argument for this command is `experiment_name`.
- kill_jobs: All jobs corresponding to an experiment will be killed form slurm for a specific experiment. 
  The argument for this command is `experiment_name`.

The arguments can also be viewed with `python hyperopt.py command --help`.

## Developer Section
 The code structure is displayed in the following chart:
```buildoutcfg
├── demo
│   └── test_conf.yaml
├── README.md
├── setup.py
└── sigopt_hyperopt
    ├── experiments
    │   ├── experiment.py
    │   └── __init__.py
    ├── hyperopt.py
    ├── __init__.py
    ├── utils
    │   ├── conda_utils.py
    │   ├── configs.py
    │   ├── const.py
    │   ├── entry_points.py
    │   ├── git_utils.py
    │   ├── __init__.py
    │   ├── logger_utils.py
    │   ├── sigopt_utils.py
    │   ├── slurm_utils.py
    │   └── workspace_utils.py
    └── workers
        ├── __init__.py
        └── worker.py

```

The main entry for this script is `sigopt_hyperopt/hyperopt.py` which will run the corresponding command specified in 
`sigopt_hyperopt/utils/entry_points.py`. 

Next, there are two important files called `sigopt_hyperopt/experiments/experiment.py` and 
`sigopt_hyperopt/workers/worker.py`. The former includes the experiment class which deals with everything experiment 
related logic like experiment initialization, deletion of the experiment (includes workspace and sigopt experiment), or 
the creation of the bash script. Moreover, it has subclasses for each cluster environment - the subclass 
`BWUniClusterExperiment` is only a dummy subclass. The latter implements the worker job that gets deployed by slurm. 
It handles the communication between the worker and Sigopt and thus, also runs the evaluation for each suggestion.

Furthermore, the folder `utils` includes different helper files for specific applications.