# take_elementary_step
Given a starting structure (defined as a SMILES string) the program generates all possible structures that can be
generated via an elementary step. An elementary step, as defined by [Zimmerman](http://dx.doi.org/10.1002/jcc.23271) is 
an elementary reaction step as a chemical rearrangement "with no more than two connections breaking and two connections forming simultaneously, while maintaining the upper and lower limits for coordination number of each atom." Two structures that are related in this way are very likely to be connected by a single transition state.

The implementation is based on the atom connectivity approach by [Woo Youn Kim and coworkers](http://dx.doi.org/10.1039/C7SC03628K).
take_steep.py generates all such states, using the atom connectivity approach by [Woo Youn Kim and coworkers](http://dx.doi.org/10.1039/C7SC03628K), for a given molecule and then selects the "best ones" based on a crude energy function based on bond dissociation energies and writes out the corresponding GFN-xTB input files (modified xyz files).

rank_steps.py reads in the GFN-xTB output files and selects the "best ones" based on the computed energy.

usage
1. Define name and SMILES of starting structure and some other parameters in  take_steps.py 
2. python take_steps.py > input_smiles
3. perform the GFN-xTB (or other) calculations
4. Define directory and some other parameters in rank_steps.py
5. python rank_steps.py > output_smiles
