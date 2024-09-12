"""
Playground which reads in the unique ligand db and lets you play with it.
"""
from pathlib import Path
from DARTassembler.src.ligand_extraction.DataBase import LigandDB
from tqdm import tqdm

if __name__ == '__main__':


    new_path = Path('/Users/timosommer/Downloads/test_DART/speedup_metalig/data_output/metalig_3000.jsonlines')

    metalig = LigandDB.load_from_json(path=new_path, n_max=50)

    # Update the MetaLig
    # for uname, ligand in tqdm(metalig.db.items()):
        # ligand.has_planar_donors = ligand.planar_check()


    # Save to .jsonlines file
    outfile = new_path
    metalig.save_to_file(outfile)

    #%% ==============    Doublecheck refactoring    ==================
    from dev.test.Integration_Test import IntegrationTest
    old_dir = Path(outfile.parent.parent, 'benchmark_data_output')
    if old_dir.exists():
        test = IntegrationTest(new_dir=outfile.parent, old_dir=old_dir)
        test.compare_all()
        print('Test for assembly of complexes passed!')
    else:
        print(f'ATTENTION: could not find benchmark folder "{old_dir}"!')




    print('Done')
