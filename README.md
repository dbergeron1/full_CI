# full_CI

Full configuration-interaction for the 2D extended Hubbard model with local, natural and Hartree-Fock orbitals.

This code computes the ground state using local, natural, and Hartree-Fock spin-orbitals, for the two-dimensional extended Hubbard model, including local, nearest neighbor, and next-nearest neighbor interaction, and next-nearest neighbor hopping.

For each type of orbitals, it plots the wave function coefficients as a function of Slater determinant (SD) index, which is ordered in increasing number of excitations with respect to the reference SD. One can then observe how fast the convergence is with respect to excitation number for the three different choices of spin-orbitals. In particular, one can observe that with natural orbitals, even at very strong coupling, the coefficient of the reference SD still has a substantially larger magnitude that all other coefficients, although relatively smaller than at weak coupling.

The parameter file is "full_CI.dat". You can change the parameters but cannot remove any line or parameter. The order in which the lines appear and the number of parameters on each line are fixed. 

To draw the lattice, use '*' and and spaces. 

"pos_u" and "pos_d" are the indices of the spin-up and spin-down particles in the reference determinant for the local basis computation. The indices increase in the left to right and top to bottom directions of the lattice. For the Hubbard model at half-filling, an antiferromagnetic configuration typically yields a fast convergence of the Lanczos algorithm.

"Niter_Lanczos" is the maximum number of iterations in the Lanczos procedure.

"Niter_max_HF" is the maximum number of iterations in the computation of the Hartree-Fock orbitals, and "tol_HF" is the tolerance on the magnitude of the off diagonal elements of the Fock matrix. "r_init_HF" and "Dr_HF" are parameters helping the Hartree-Fock orbitals computation to converge. They must respect the relation r_init_HF+n*Dr_HF=1, where n is an integer.

The figures are displayed if "display_figures" is set to "yes"

"random_init_state" can be set to "yes" to use a random initial state in the Lanczos procedure of the local basis computation, instead of the reference determinant as the initial state, for instance if the Lanczos routine converges too slowly with the latter choice.
