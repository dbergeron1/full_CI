# full_CI

Full configuration-interaction for the 2D extended Hubbard model with local, natural and Hartree-Fock orbitals.

This code computes the ground state using local, natural, and Hartree-Fock spin-orbitals, for the two-dimensional extended Hubbard model, including local, nearest neighbor, and next-nearest neighbor interaction, and next-nearest neighbor hopping.

For each type of orbitals, it plots the wave function coefficients as a function of Slater determinant (SD) index, which is ordered in increasing number of excitations with respect to the reference SD. One can then observe how fast the convergence is with respect to excitation number for the three different choices of spin-orbitals. In particular, one can observe that with natural orbitals, even at very strong coupling, the coefficient of the reference SD still has a substantially larger magnitude that all other coefficients, although smaller relatively to other SD's than at weak coupling.

The parameter file is "full_CI.dat". You can change the parameters but cannot remove any line or parameter. The order and number of parameters on each line is fixed. 

The lattice is drawn using '*' and and spaces. 
