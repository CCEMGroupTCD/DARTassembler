import DARTassembler

Br_filter = DARTassembler.ligandfilters('input/ligandfilters_Br.yml', nmax=100)
phenyl_filter = DARTassembler.ligandfilters('input/ligandfilters_phenyl.yml', nmax=100)
P_N_filter = DARTassembler.ligandfilters('input/ligandfilters_P_N_ligands.yml')

assembly = DARTassembler.assembler('input/Pd_Ni_assembly_input.yml')


