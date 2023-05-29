with open("/Users/cianclarke/Documents/PhD/Complex_Assembly/CreateTMC/tmp/CASINI_test_220523/AuCl2_CSD-JEKZEI-02-a/AuCl2_CSD-JEKZEI-02-a.log", "r") as f:
    lines = f.readlines()

start = "SECOND ORDER PERTURBATION THEORY ANALYSIS OF FOCK MATRIX IN NBO BASIS"
end = "NATURAL BOND ORBITALS (Summary):"
start_found = False
end_found = False
result_N_found = False
result_C_found = False

carbon_interaction_sum = 0
for idx, line in enumerate(lines):


    if start in line:
        start_found = True
    if end in line:
        end_found = True



    if start_found and not end_found:
        line_list = line.split()
        # NITROGEN
        # Here we make sure that LP comes before N which comes before Au which should come before Cl
        # if (")Au" in line) and ("-Cl" in line) and ("LP" in line) and (len(line_list) > 3) and ("N" in line) and (line.find("LP") < line.find("N")) and (line.find("N") < line.find("Au")) and (line.find("N") < line.find("Cl")):
        if (")Au" in line) and ("-Cl" in line) and ("LP" in line) and (len(line_list) > 3) and ("N" in line) and (line.find("LP") < line.find("N")) and (line.find("N") < line.find("Au")) and (line.find("N") < line.find("Cl")):
            interaction_energy_N = line_list[-3]
            print(line.strip() + f"         Interaction_N: [{interaction_energy_N}]")
            result_N_found = True
        else:
            pass

        # CARBON
        # if (")Au" in line) and ("-Au" in line) and ("BD" in line) and (len(line_list) > 3) and ("C" in line) and ("Cl" not in line) and (line.find("BD") < line.find("C")) and (line.find("C") < line.find("-Au")) and (line.find("-Au") < line.find(")Au"))
        if (")Au" in line) and ("BD" in line) and ("RY" in line) and (len(line_list) > 3) and ("- C" in line) and (") C" not in line) and ("Cl" not in line)  and ("N" not in line) and (line.find("BD") < line.find("- C")):
            interaction_energy_C = float(line_list[-3])

            if interaction_energy_C >= 5.0:
                carbon_interaction_sum += interaction_energy_C
            print(line.strip() + f"         Interaction_C: [{interaction_energy_C}]")
            result_C_found = True
        else:
            pass

if not result_N_found:
    print("!!!Fatal_Error!!!-> No value found for Au-N")
    raise ValueError

if not result_C_found:
    print("!!!Fatal_Error!!!-> No value found for Au-C")
    raise ValueError

print(carbon_interaction_sum)
