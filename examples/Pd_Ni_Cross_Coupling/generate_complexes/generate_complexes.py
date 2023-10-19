import DARTassembler


nmax = False
Br_filter = DARTassembler.filter_ligands('input/ligandfilters_Br.yml', nmax=100)
phenyl_filter = DARTassembler.filter_ligands('input/ligandfilters_phenyl.yml', nmax=100)
P_N_filter = DARTassembler.filter_ligands('input/ligandfilters_P_N_ligands.yml', nmax=nmax)

assembly = DARTassembler.assemble_complexes('input/Pd_Ni_assembly_input.yml')

print('Done! All complexes assembled.')

