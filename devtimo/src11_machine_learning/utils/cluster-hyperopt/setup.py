from setuptools import setup, find_packages

setup(
        name="sigopt_hyperopt",
        version='1.0.0',
        author="marcel.hiltscher@student.kit.edu, matthias.schniewind@kit.edu",
        description="Hyperparameter Optimization on the Cluster using SigOpt as optimization server.",
        packages=find_packages(),
        package_data={"sigopt_hyperopt": ["py.typed"]},
        install_requires=[
                "sigopt",
                "GitPython",
                "yamldataclassconfig"
        ],
        python_requires=">=3.7",
)
