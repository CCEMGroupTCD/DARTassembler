from pathlib import Path

import pandas as pd

from DARTassembler.src.ligand_extraction.DataBase import LigandDB
import warnings
warnings.filterwarnings("ignore", category=UserWarning)


if __name__ == '__main__':

    n_max = 1500  # None or int. Maximum number of ligands to load
    denticities = None  # None for all denticities or list of denticities
    concat_outdir = 'concat_xyz'
    remove_haptic = False
    sort_by_rssd = False
    output_all_isomers = True



    # ==========  Main script  ==========

    concat_outdir = Path(concat_outdir)
    # Load MetaLig database
    db = LigandDB.load_from_json(n_max=n_max)
    if remove_haptic:
        db.db = {name: ligand for name, ligand in db.db.items() if not ligand.has_neighboring_coordinating_atoms}
    db = db.get_db_with_only_certain_denticities(denticities=denticities)

    data = db.save_ligand_geometry_concat_xyz_files(
                                            outdir=concat_outdir,
                                            sort_by_rssd=sort_by_rssd,
                                            output_all_isomers=output_all_isomers
                                            )

    # Load data into pandas DataFrame
    long_data = []
    for geometry, info in data.items():
        for name, rssd, _, _, _, idc in info:
            long_data.append([geometry, name, rssd, idc])
    df = pd.DataFrame(long_data, columns=['geometry', 'name', 'rssd', 'idc'])

    # # Plot distribution of rssd values
    # import seaborn as sns
    # import matplotlib.pyplot as plt
    # plot_outdir = 'plots'
    # plot_outdir = Path(plot_outdir)
    # for geometry, names_rssd in data.items():
    #     plt.figure()
    #     _, rssds, _, _, _ = zip(*names_rssd)
    #     sns.histplot(rssds, label=geometry, bins=15)
    #     plt.title(geometry)
    #     plt.xlabel('RSSD')
    #     plt.ylabel('Count')
    #     outpath = Path(plot_outdir, f'{geometry}.png')
    #     plt.savefig(outpath)
    #     plt.close()

    # ==============    Doublecheck refactoring    ==================
    from dev.test.Integration_Test import IntegrationTest
    new_dir = concat_outdir
    old_dir = concat_outdir.parent / f'OLD_n={n_max}_concat_xyz'
    if old_dir.exists():
        test = IntegrationTest(new_dir=new_dir, old_dir=old_dir)
        test.compare_all()
        print('Test for ligand geometries passed!')
    else:
        print(f'ATTENTION: could not find benchmark folder "{old_dir}"!')

    print('Done!')