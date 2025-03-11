"""
Playground which reads in the unique ligand db and lets you play with it.
"""
from pathlib import Path
from DARTassembler.src.ligand_extraction.DataBase import LigandDB
from tqdm import tqdm

if __name__ == '__main__':


    outfile = Path('/Users/timosommer/PhD/projects/DARTassembler/DARTassembler/data/metalig/test1000_MetaLigDB_v1.1.0.jsonlines')
    n_max = 1000

    metalig = LigandDB.load_from_json(n_max=n_max)

    if outfile.exists():
        raise FileExistsError(f"Output path already exists. Please specify another output path.")

    # Save to .jsonlines and .csv file
    metalig.save_to_file(outfile)
    metalig.save_reduced_csv(outpath=Path(outfile).with_suffix('.csv'))

    # Read in again the new file
    new_metalig = LigandDB.load_from_json(outfile)

    #%% ==============    Doublecheck refactoring    ==================
    # from dev.test.Integration_Test import IntegrationTest
    # old_dir = Path(outfile.parent.parent, 'benchmark_data_output')
    # if old_dir.exists():
    #     test = IntegrationTest(new_dir=outfile.parent, old_dir=old_dir)
    #     test.compare_all()
    #     print('Test for assembly of complexes passed!')
    # else:
    #     print(f'ATTENTION: could not find benchmark folder "{old_dir}"!')




    print('Done')
