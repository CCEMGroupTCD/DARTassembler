from process import run

if __name__ == '__main__':
    """
    Runs the full extraction
    TestSize: Decide if we want to run the programm on a testset (1/100) of the actual size
    Metals_of_interest: MetalAtom in TMC we want to extract the ligands from in the TMQM DB
    denticity_numbers: Denticities we are interested in
    
    """

    run(TestSize=False,
        metals_of_interest=["Fe", "Mn", "Cr", "Ru", "Co", "Ni", "Cu", "Ir", "Mo"],
        denticity_numbers=[2, 3, 4, 5],
        get_csd_Ase_dict=False,
        Duplicate_Filter=True
        )