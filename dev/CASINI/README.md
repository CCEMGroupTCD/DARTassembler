# CASINI CALCS
Here we must store all the data to do with CASINI and double check all the data to 

All original calculations can be found in:
/Users/cianclarke/Documents/PhD/Complex_Assembly/CreateTMC/Gaussian/CASINI_DATA/CASINI_230523_calcs

## Ligand Filtering

1. Filter.denticity_of_interest_filter(denticity_of_interest=[2])

2. Filter.filter_neighbouring_coordinating_atoms()

3. Filter.filter_charge_confidence(filter_for="confident")

4. Filter.filter_coordinating_group_atoms(denticity=2, atoms_of_interest=["N", "C"], instruction="must_contain_and_only_contain")

5. Filter.filter_even_odd_electron(filter_for="even")

6. Filter.filter_ligand_charges(denticity=2, charge=-1)

- Build LigandDatabase: 72082

- Before Filtering: 72082

- After Denticity Filter: 24702

- After Neighbouring Coordinating Atoms: 23584

- Charge Assignment Confident: 15140

- Charge Assignment Not Confident: 8444

- After Charge Confidence Filter: 15140

- After coordinating atom filter: 1152

- After even / odd electron filter: 1152

- After ligand charges filter: 894

- End: 894

- Filtering Done







