from constants.Paths import project_path
from itertools import zip_longest
import numpy as np

file1name = str(project_path().extend("src14_Assembly_Unit_Test", "INTEGRATION_TEST_Benchmark.xyz"))
file2name = str(project_path().extend("src14_Assembly_Unit_Test", "INTEGRATION_TEST.xyz"))

with open(file1name) as file1, open(file2name) as file2:
    i = 1
    for line1, line2 in zip_longest(file1, file2):
        line1 = str(line1).split()
        line2 = str(line2).split()
        if len(line1) == len(line2):
            if len(line1) == 4:
                p0 = np.array([float(line1[-3]), float(line1[-2]), float(line1[-1])])
                p1 = np.array([float(line2[-3]), float(line2[-2]), float(line2[-1])])
                dist = np.linalg.norm(p0 - p1)
                if str(line1[-4]) != str(line2[-4]):
                    print(f"!!!Fatal Error!!! --> The atom types have changed [{str(line1[-4])}] != [{str(line2[-4])}] on line [{i}] --> The code is not working as expected")
                    raise AssertionError

                elif dist > 0.1:
                    print(f"!!!Fatal Error!!! --> The atom positions have changed [{p0}] != [{p1}] on line [{i}] --> The code is not working as expected")
                    raise AssertionError
                else:
                    pass
            else:
                pass
        else:
            print(f"!!!Fatal Error!!! --> The length of the following lines are not equivalent --> The code is not working as expected")
            print(line1)
            print(line2)
            raise AssertionError
        i += 1
print("!!!Success!!! --> The code still works --> you can pat yourself on the back")
