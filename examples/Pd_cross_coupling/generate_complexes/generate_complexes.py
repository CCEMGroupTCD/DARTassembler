"""
This script generates the complexes for screening the cross-coupling reaction as an example for the DART paper.
"""
import DARTassembler


nmax = False
# Br_filter = DARTassembler.filter_ligands('./ligandfilters_Br.yml', nmax=nmax)
# phenyl_filter = DARTassembler.filter_ligands('./ligandfilters_phenyl.yml', nmax=nmax)
# P_N_filter = DARTassembler.filter_ligands('./ligandfilters_P_N_ligands.yml', nmax=nmax)


assembly = DARTassembler.assemble_complexes('Pd_Ni_assembly_input.yml')

print('Done! All complexes assembled.')

