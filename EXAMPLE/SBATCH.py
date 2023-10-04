###############################
# The Goal of this
# program is to iterate through
# each directory and SBATCH all
# the .sh files present.
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
        counter = 0
        for batch in self.batch_dir:
            # Iterate entries
            for child in batch.iterdir():
                if child.is_file():
                    pass
                if child.is_dir():
                    for entries in child.iterdir():
                        if entries.is_dir() and str(entries.name) == "complexes":
                            for calculation in entries.iterdir():
                                for file in calculation.iterdir():
                                    if file.is_dir():
                                        raise ValueError
                                    elif str(file.name) == "run.sh":
                                        counter += 1
                                        print(f"sbatch {str(file.parent)}/run.sh")
                                        os.chdir(str(file.parent))
                                        #os.system(f"sbatch {str(file.parent)}/run.sh")
                        else:
                            pass


                else:
                    print(f"!!!Fatal_Error!!! --> Unexpected directory [{child.name}] found in [{child.parent}]")
                    raise ValueError
        print(f"Total Calcs: [{counter}]")








if __name__ == "__main__":
    # It is okay to use absolute paths in this case because the program will run exclusively on HPC systems
    calc_iterator = Iterate(global_batch_directory="/Users/cianclarke/Documents/PhD/Complex_Assembly/DART/src14_Assembly_Unit_Test/DART_Example_Pd_Ni_Complexesqs"
                                                   "")
    calc_iterator._iter_calcs()
    pass
