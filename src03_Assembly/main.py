from src03_Assembly.RandomComplexAssembler import RandomComplexAssembler
import pickle


if __name__ == "__main__":

    RCA = RandomComplexAssembler(database_path="../data/LigandDatabases/ligand_db_w_filters.pickle",
                                 store_path="../data/Assembled_Molecules")

    safe_path = "../data/Assembled_Molecules"
    successfull_generated_complexes = {}

    for i in range(10):
        random_complex = RCA.create_random_TMC(visualize_=True, optimize_=False)

        if random_complex is None:
            print("Thus no complex could be generated")
        else:
            random_complex.print_to_xyz(path=safe_path)
            successfull_generated_complexes[i] = random_complex
        input("press Enter to continue")

    with open(f"{safe_path}/00_Random_complex_dict.pickle", "wb") as handle:
        pickle.dump(successfull_generated_complexes, handle)

    print("done")
