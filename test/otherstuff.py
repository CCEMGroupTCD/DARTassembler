import pickle

if __name__ == '__main__':
    #unittest.main()

    with open("test_xyz_files/DUCVIG.xyz") as file:
        lines = file.readlines()

        dict_ = {}
        for j, line in enumerate(lines):
            if j > 1:
                dict_[j-2] = [line.split()[0], [float(line.split()[i]) for i in range(1,4)]]

        new_dict = {"DUCVIG" : dict_}

    with open(f"../test/other_test_files/coordinate_dict_example.pickle", "wb") as file:
        pickle.dump(new_dict, file)

    print("done")