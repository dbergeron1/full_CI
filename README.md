# full_CI
Full configuration-interaction for the 2D extended Hubbard model with natural and Hartree-Fock orbitals.

This code computes the ground state using local, natural, and Hartree-Fock spin-orbitals, for the two-dimensional extended Hubbard model, including local, nearest neighbor, and next-nearest neighbor interaction, and next-nearest neighbor hopping.

It also plots the wave-function coefficients as a function of Slater determinant (SD) index, which is ordered in increasing number of excitations with respect to the reference SD. One can therefore observe how fast the convergence is with respect to excitation number for the three different choices of spin-orbitals.

The parameter file is "full_CI.dat". You can change the parameters but do not remove any line or parameter. The order and number of parameters on each line is fixed. 

The lattice is drawn using '*' and and spaces. 
