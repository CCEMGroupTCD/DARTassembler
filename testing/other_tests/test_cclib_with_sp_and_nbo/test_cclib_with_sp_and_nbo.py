from DARTassembler.src.assembled_complex_output.gaussian import GaussianOutput
import numpy as np
import pandas as pd

if __name__ == '__main__':

    # Test the GaussianOutput class
    fchk_path = '/Users/timosommer/PhD/projects/OERdatabase/dev/timo/playground/complex_data/dftsp_output/HOQOKUME_Mn_OH/HOQOKUME_Mn_OH_structure_gfnffrot_gfn2opths_dftspls/HOQOKUME_Mn_OH_structure_gfnffrot_gfn2opths_dftspls_onlyNBO.log'  # single point calculation
    # fchk_path = '/Users/timosommer/PhD/projects/RCA/projects/DART/examples/Pd_Ni_Cross_Coupling/dev/output/gaussian_relaxed_complexes/batches/P_N_Donors_Ni_Metal_Centre/complexes/ABADEZIX_PN_Ni/ABADEZIX_PN_Ni_gaussian.fchk'  # relaxed complex

    go_nbo = GaussianOutput(fchk_path)
    nbo_data = go_nbo.get_all_raw_data()

    sp_path = fchk_path.replace('onlyNBO', 'onlySP')
    go_sp = GaussianOutput(sp_path)
    sp_data = go_sp.get_all_raw_data()

    full_path = fchk_path.replace('_onlyNBO', '')
    go_full = GaussianOutput(full_path)
    full_data = go_full.get_all_raw_data()

    # Check if the full data  is equal to the NBO or SP data
    nbo_equal = full_data == nbo_data
    sp_equal = full_data == sp_data
    print(f'Full-NBO:\t {nbo_equal}')
    print(f'Full-SP:\t {sp_equal}')


    # Check if the NBO and SP data are the same
    dict1 = sp_data
    dict2 = full_data
    for key in dict1.keys():
        if key == 'atomspins':
            keys = set(dict1[key].keys()).union(dict2[key].keys())
            for key2 in keys:
                equal = np.allclose(dict1[key][key2], dict2[key][key2])
                print(f'{key}-{key2}:\t {equal}')
        elif key == 'atomcharges':
            keys = set(dict1[key].keys()).union(dict2[key].keys())
            for key2 in keys:
                if not key2 in dict1[key].keys():
                    print(f'{key}-{key2}:\t {key2} not in NBO data')
                    continue
                if not key2 in dict2[key].keys():
                    print(f'{key}-{key2}:\t {key2} not in SP data')
                    continue
                equal = np.allclose(dict1[key][key2], dict2[key][key2])
                print(f'{key}-{key2}:\t {equal}')
        elif key == 'moments':
            for idx, (a1, a2) in enumerate(zip(dict1[key], dict2[key])):
                equal = np.allclose(a1, a2)
                print(f'{key}-{idx}:\t {equal}')
        else:
            try:
                equal = np.allclose(dict1[key], dict2[key])
            except:
                try:
                    equal = dict1[key] == dict2[key]
                except:
                    equal = 'undefined'
            print(f'{key}:\t {equal}')