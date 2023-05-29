print("########################################################################################################")
import os
import re
if __name__ == "__main__":
    root = "/Users/cianclarke/Documents/PhD/Complex_Assembly/CreateTMC/tmp/FRANK_200423"
    list_of_directories = os.listdir(root)
    all_calculations = []
    stationary_point_found = []
    NBO_completed = []
    failed_calculation = []
    calculation_not_started =[]
    for calculation in list_of_directories:
        if calculation.startswith("AuCl2"):
            all_calculations.append(calculation)
            calc_files = os.listdir(root+"/"+calculation)
            slurm_out_present = False
            tmp_dic = {}
            for file in calc_files:
                res = len(re.findall('(?=(_cont))', file))
                tmp_dic.update({file: int(res)})
                if str(file).endswith(".out"):
                    slurm_out_present = True
                else:
                    pass

            for file in calc_files:
                if str(file).endswith(".log") and slurm_out_present and (len(re.findall('(?=(_cont))', file)) == max(tmp_dic.values())):
                    file1 = open(f'{root+"/"+calculation+"/"+file}', 'r')
                    Lines = file1.readlines()
                    opt_correct = False
                    NBO_correct = False
                    for line in Lines:
                        if "-- Stationary point found" in line:
                            #If optimization has completed succesfully
                            stationary_point_found.append(calculation)
                            opt_correct = True
                        elif "NBO analysis completed" in line:
                            #If NBO has completed successfully
                            NBO_completed.append(calculation)
                            NBO_correct = True
                    if (opt_correct is False) or (NBO_correct is False):
                        #If there is a failure
                        print(f"CALCULATION [{calculation}] has failed: Optimization completed [{opt_correct}], NBO completed [{NBO_correct}]")
                        failed_calculation.append(calculation)
                else:
                    calculation_not_started.append(calculation)
        else:
            pass
    print("\nThis Data is obtained by searching for the exit code in the *.out file")
    print(f"Total submitted calculations: {len(all_calculations)}")
    print(f"Total Optimization Completed: {len(stationary_point_found)}")
    print(f"Total NBO_completed: {len(NBO_completed)}")
    print(f"Total Failures: {len(failed_calculation)}")
    print(f"Percentage Compete: {100*round((len(NBO_completed)+len(failed_calculation))/len(all_calculations),2)}\n")
    print("########################################################################################################")