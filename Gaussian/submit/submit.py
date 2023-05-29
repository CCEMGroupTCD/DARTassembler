import os
import re
import subprocess
# os.walk('.') generates the file names in a directory tree by walking the tree
# either top-down or bottom-up. For each directory in the tree rooted at directory top
# (including top itself), it yields a 3-tuple (dirpath, dirnames, filenames).
# Iterate over each file in any subdirectory of the current directory
print("##############################################################################################################")
for root, dirs, files in os.walk('/Users/cianclarke/Documents/PhD/Complex_Assembly/CreateTMC/tmp/CASINI_test_220523'):
    # Iterate through all files in the current directory
    tmp_dic = {}
    for file in files:
        res = len(re.findall('(?=(_cont))', file))
        tmp_dic.update({file:int(res)})

    for file in files:
        # Check if the file has a ".com" extension
        if file.endswith('.com') and (len(re.findall('(?=(_cont))', file)) == max(tmp_dic.values())) and (len(re.findall('(?=(_cont))', file)) != 0):
            filename = os.path.expanduser(os.path.join('~', root, file))
            file_without_ext = os.path.splitext(filename)[0]
            name = os.path.basename(file_without_ext)
            os.chdir(os.path.dirname(file_without_ext))
            cwd = os.getcwd()
            print(f"current working directoy is {cwd}")
            print(f"subprocess.run(['submit-g16.sh', 'amd', {name}.com, '20', '72'])")
            #result = subprocess.run(['submit-g16.sh', 'amd', name, '20', '72'])
            #print(result.stdout)
            os.chdir('..')
        else:
            pass
print("##############################################################################################################")