from src06_fullRun.full_run import Run

if __name__ == "__main__":

    # setup for the assembly
    run = Run(database="tmQM",
              testing=15000,
              calculate_charges=True,
              get_charges=True,
              batch_yml_path="example_batch.yml",
              gaussian_input_specifications="gaussian_specifications.yml"
              )

    # for debugging
    # assemble the complexes
    # assembled_complexes = run.assembled_complexes

    # c = assembled_complexes[0]

    # c.to_json("test.json")
    # from src03_Assembly_Cian.TransitionMetalComplex import TransitionMetalComplex
    # c2 = TransitionMetalComplex().from_json("test.json")

    print(2)