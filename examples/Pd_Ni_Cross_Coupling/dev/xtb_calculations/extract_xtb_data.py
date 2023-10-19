from DARTassembler.src.assembly.assembled_complex_output import AssembledComplexOutput
import pandas as pd
from pathlib import Path

if __name__ == '__main__':
    all_structures_dir = '231013_relaxations_and_original_dirs'


    # Extract data from xtb calculations
    data = []
    dirs = AssembledComplexOutput.get_all_complex_directories(all_structures_dir)
    for dir in dirs:
        cout = AssembledComplexOutput(dir)
        complex = cout.xtbcomplex
        bond_lengths = {f'dist_{el}': complex.get_donor_metal_bond_length(el) for el in complex.donor_elements}
        data.append({
                        'complex': cout.name,
                        'P-N bite angle': complex.get_bite_angle(['P', 'N']),
                        **bond_lengths
                    })
    df = pd.DataFrame(data)
    df.to_csv(Path(all_structures_dir, 'xtb_data.csv'), index=False)
    print('Done!')

