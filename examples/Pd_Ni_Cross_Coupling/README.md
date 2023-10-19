# How To: Targeted assembly of complexes with precise control over the ligands

This example shows you how to assemble complexes with DART and using precise control over the assembled ligands. This is particularly helpful if you are interested in a particular chemistry such as cross-coupling reactions. This example was also discussed in our initial DART paper.

In this example, we will assemble square-planar Pd and Ni complexes with two monodentates and one bidentate ligands, i.e. a (2,1,1) topology: a bromine, a phenyl ligand and a bidentate P-N donor. In all complexes, the bromine and the phenyl will always be the same, but we will assemble complexes with all neutral P-N donors in the MetaLig database. This allows us to explore the effect of the P-N donor on the properties of the final complex.

## Step 1: Create three different ligand databases
In order to have precise control over which ligands are assembled 