import re
check_lines = False
with open("FeTMP2-newreact-bent-NBO.log", "r") as f:
        lines = f.readlines()
        for idx, line in enumerate(lines):
                if "(2)" in line:
                        #interactions are printed in the file
                        check_lines = True
                if check_lines:
                        get_nums = re.findall(r"[-+]?\d*\.\d+|\d+",line)
                        if len(get_nums)>3:
                                #then there are numbers here, and we can condition so that we only see certain atoms interactions
                                if float(get_nums[-3])>5:# and ("Al 78" in line or "Al 79" in line):
                                        #Then the interaction is strong
                                        #if
                                        print(line)
                        #print(re.findall(r"[-+]?\d*\.\d+|\d+",line))
                if "NATURAL BOND ORBITALS (Summary):" in line:
                        #now we can stop looking at the lines
                        check_lines = False
