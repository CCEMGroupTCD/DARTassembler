###############################
# The Goal of this
# program is to provide a
# progress update on the
# calculations
###############################
import pathlib
import os
class Iterate:
    def __init__(self, global_batch_directory: str = None):

        assert global_batch_directory is not None
        self.global_path = pathlib.Path(global_batch_directory)
        assert self.global_path.is_dir()
        self.batch_dir = []
        self._get_batch_dir()
        assert self.batch_dir != []

    def _get_batch_dir(self):
        for child in self.global_path.iterdir():
            if child.is_file():
                pass
            elif child.is_dir():
                self.batch_dir.append(child)

    def _iter_calcs(self):
        # Iterate Batches
        total_calcs_counter = 0
        completed_calcs = 0
        for batch in self.batch_dir:
            # Iterate entries
            for child in batch.iterdir():
                if child.is_file():
                    pass
                if child.is_dir():
                    for entries in child.iterdir():
                        if entries.is_dir() and str(entries.name) == "complexes":
                            for calculation in entries.iterdir():
                                total_calcs_counter += 1
                                for file in calculation.iterdir():
                                    if file.is_dir():
                                        raise ValueError
                                    elif str(file.name).endswith(".log"):
                                        print(f"{str(file.name)}")
                                        if self._check_progress(file):
                                            completed_calcs += 1
                                        else:
                                            pass


                        else:
                            pass
                else:
                    print(f"!!!Fatal_Error!!! --> Unexpected directory [{child.name}] found in [{child.parent}]")
                    raise ValueError
        print(f"Total Calcs: [{total_calcs_counter}]")
        print(f"Completed Calcs: [{completed_calcs}] ")
        print(f"Percentage Completed: [{(completed_calcs/total_calcs_counter)*100}]")

    @staticmethod
    def _check_progress(file: pathlib):
        calc_completed = False
        with open(file, "r") as log_file:
            for line in log_file:
                if "-- Stationary point found" in str(line):
                    calc_completed = True
                    break
                else:
                    #Stationary point not found
                    pass
        return calc_completed







if __name__ == "__main__":
    # It is okay to use absolute paths in this case because the program will run exclusively on HPC systems
    calc_iterator = Iterate(global_batch_directory="/Users/cianclarke/Documents/PhD/Complex_Assembly/DART/src14_Assembly_Unit_Test/DART_Example_Pd_Ni_Complexes/")
    calc_iterator._iter_calcs()
    pass