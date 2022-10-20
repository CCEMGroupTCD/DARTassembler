from src03_Assembly.RandomComplexAssembler import RandomComplexAssembler


if __name__ == "__main__":

    RCA = RandomComplexAssembler(database_path="../data/LigandDatabases/ligand_db_w_filters.pickle",
                                 store_path="../data/Assembled_Molecules")

    for i in range(10):
        random_complex = RCA.create_random_TMC(visualize_=True, optimize_=False)

        if random_complex is None:
            print("Thus no complex could be generated")
        input("press Enter to continue")

    print("done")
