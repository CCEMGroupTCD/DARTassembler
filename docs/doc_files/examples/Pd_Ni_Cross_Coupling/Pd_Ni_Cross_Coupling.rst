How to: Assemble Pd/Ni complexes for Cross-Coupling screening (exact ligand control)
====================================================================================

In this example, we will show you how to assemble Pd and Ni complexes for cross-coupling screening. This is particularly helpful if you are interested in a particular chemistry such as cross-coupling reactions. This example was also discussed in our initial DART paper.

In this example, we will assemble square-planar Pd and Ni complexes with two monodentates and one bidentate ligands, i.e. a 2-1-1 topology: a bromine, a phenyl ligand and a bidentate P-N donor. In all complexes, the bromine and the phenyl will always be the same, but we will assemble complexes with all neutral P-N donors in the MetaLig database. This allows us to explore the effect of the P-N donor on the properties of the final complex.

This is the python code that needs to be executed to assemble the complexes:

.. code-block:: python

    import DARTassembler
    Br_filter = DARTassembler.filter_ligands('input/ligandfilters_Br.yml')
    phenyl_filter = DARTassembler.filter_ligands('input/ligandfilters_phenyl.yml')
    P_N_filter = DARTassembler.filter_ligands('input/ligandfilters_P_N_ligands.yml')

    assembly = DARTassembler.assemble_complexes('input/Pd_Ni_assembly_input.yml')




