How to: Assemble Pd/Ni complexes for Cross-Coupling screening (exact ligand control)
====================================================================================

.. _Pd_Ni_Cross_Coupling:

In this example, we will show you how to assemble Pd and Ni complexes for cross-coupling screening. This is particularly helpful if you are interested in a particular chemistry such as cross-coupling reactions. This example was also discussed in our initial DART paper.

In this example, we will assemble square-planar Pd and Ni complexes with two monodentates and one bidentate ligands, i.e. a 2-1-1 topology: a bromine, a phenyl ligand and a bidentate P-N donor. In all complexes, the bromine and the phenyl will always be the same, but we will assemble complexes with all neutral P-N donors in the MetaLig database. This allows us to explore the effect of the P-N donor on the properties of the final complex.

To run the example, lets start by making a new directory for this example and copying the files we need into it:

All the files needed to run this example are in the github repository in the folder https://github.com/CCEMGroupTCD/DART/tree/master/examples/Pd_Ni_Cross_Coupling/generate_complexes. Please download these files as zip and extract them. We will go through them one by one and explain what they do.



To execute the python code in this file, please open the terminal in the folder containing this file and make sure you have imported DART. Then, run ``python generate_complexes.py`` in the terminal.


