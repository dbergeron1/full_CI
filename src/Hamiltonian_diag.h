//
//  Hamiltonian_diag.h
//  full_CI
//
//  Created by Dom2 on 20-04-09.
//

#ifndef Hamiltonian_diag_h
#define Hamiltonian_diag_h

#include "SO_basis_Hamiltonian.h"
#include <algorithm>
#include <set>
#include "graph_2D.h"

enum operator_type {part, hole};
enum orbital_type_T {local,natural,HartreeFock};

template<class ind_T, class val_T, class SD_T, class cfs_T>
class Hamiltonian_diag: public SO_basis_Hamiltonian<ind_T, val_T, SD_T>
{
	typedef Hamiltonian<ind_T, val_T> H_T;
	typedef SO_basis_Hamiltonian<ind_T, val_T, SD_T> SO_H_T;
public:
	Hamiltonian_diag(const array<unsigned, 2> &N_par, val_T tol_p, const coord_vec_T &lattice_par, const vector<val_T> &U_par, const vector<val_T> &t_par, const vector<val_T> &onsite_E_par);
	
	ind_T &L=H_T::L;
	array<ind_T,4> &s=H_T::s;
	unsigned &Nu=SO_H_T::Nu, &Nd=SO_H_T::Nd;
	vector<unsigned> &ph_indices=SO_H_T::ph_indices;
	vector<unsigned> &SO_ph_map=SO_H_T::SO_ph_map;
	ind_T Nph_u_max, Nph_d_max, ind_max;
	vector<vector<ind_T>> Pascal_mat;
	SD_T &ref_SD=SO_H_T::ref_SD;
	//E_ref=<ref_SD|H|ref_SD>
	val_T E_ref;
	vector<SD_T> CI_ph_basis; //CI basis expressed in the particle-hole representation relative to the reference SD. length: NSD
	vector<SD_T> CI_abs_basis; //CI basis expressed in the absolute representation. length: NSD
	map<SD_T, ind_T> CI_ph_basis_map;
	vector<short int> abs_to_ph_sign; //sign between the two representation. length: NSD
	ind_T NSD; //total number of SD's in the N-particle basis
	sparse_matrix_vec<ind_T, val_T> H_v; //Hamiltonian matrix
	vector<map<ind_T, val_T>> H_m, H_m_2; //Hamiltonian matrix used by create_H_matrix() for fast insertion of elements. Destroyed by convert_H_m_to_H_v().
	//each element in N_SD_l is a vector corresponding to a number of excitations, which elements are the cumulative numbers of SD's from 0 excitations to the pair (N_ph,N_down), inclusively. This vector is used to obtain the index of a SD.
	vector<vector<ind_T>> N_SD_l;
	
	
	vector<vector<cfs_T>> DM_up, DM_down; //single-particle density matrices
	vector<vector<val_T>> natural_basis_up, natural_basis_down;
	vector<double> DM_eig_up, DM_eig_down;
	vector<cfs_T> aL, bL; //Lanczos matrix elements
	vector<cfs_T> vL_init, vLp, vL, vLm; //Lanczos vectors
//	vector<cfs_T> vL0, vL1, vL2;
	cfs_T mod2v; //squared modulus of intermediate Lanczos vector
	ind_T iL; //Lanczos iteration index
	vector<cfs_T> GS_vec_Lanczos, GS_vec_SO; //ground-state vector in the Lanczos and spin-orbital bases, respectively
	string orbital_name;
	orbital_type_T orbital_type;
	
//	vector<cfs_T> aL_tst, bL_tst;
//	sparse_matrix_vec<ind_T, val_T> H_tst;
//	vector<cfs_T> vLp_tst, vL_tst, vLm_tst;
//	unsigned Ntst;
	
	void compute_ref_diag_E();
	
	bool create_CI_ph_basis();
	void create_list_ph(ind_T Nph_max, unsigned Ns, vector<vector<SD_T>> &list_ph);
	
	bool generate_states(vector<SD_T> &list, unsigned Nb, unsigned Nx);
	SD_T next_state(SD_T s, unsigned Nb);
	
	void print_SD(SD_T SD_p);
	ind_T get_SD_index(SD_T sd);
	ind_T get_string_index(SD_T str, unsigned Nb);
	void create_pascal_matrix();
	unsigned long binom_coef(unsigned N, unsigned n);
	
	void convert_from_ph_SD(const SD_T &SD_ph, SD_T &SD_abs, short int &sgn);
	void convert_to_ph_SD(const SD_T &SD_abs, SD_T &SD_ph, short int &sgn);
	void create_H_matrix();
	void create_H_matrix_abs();
	void create_H_matrix_ph();
	val_T elem_V(int i, int j, int k, int l);
	val_T elem_t_Fock(int i, int j);
	void insert_H_m_element_sym(ind_T i, SD_T SD_tmp, val_T val);
	void add_H_m_element(ind_T i, SD_T SD_tmp, val_T val);
	void add_H_m_element_diag(ind_T i, val_T val);
	void add_H_m_element_sym(ind_T i, SD_T SD_tmp, val_T val);
	void add_K_term_H_m(int i, int k, int l);
	void K_terms(int i, const vector<int> &v1, const vector<int> &v2);
	void add_V_term_H_m(int i, int k, int l, int m, int n);
	void V_terms(int i, const vector<int> &v1, const vector<int> &v2, const vector<int> &v3, const vector<int> &v4);
	void check_H_Hermicity();
	void convert_H_m_to_H_v();
	void compare_H_matrix();
	
	int add_bits_between(SD_T SD_tmp, int b1, int b2);
	
	bool compute_GS_vec;
	void Lanczos(unsigned NL);
	void Lanczos_init(unsigned NL, bool random_init_st=false);
	bool Lanczos_iteration();
	void eigs_Lanczos(vector<double> &eig_vals, bool compute_eig_vectors, unsigned NL);
	void compute_ground_state(orbital_type_T OT, unsigned NL_max, bool display_figs=false, bool random_init_st=false);
	void compute_natural_spin_orbitals();
	
	void test_Lanczos(unsigned Nst, unsigned NL_max, val_T a_P, val_T a_eig);
//	void Lanczos_init_tst(unsigned NL, cfs_T a_vL);
//	bool Lanczos_iteration_tst();
};

template<class ind_T, class val_T, class SD_T, class cfs_T>
Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::Hamiltonian_diag(const array<unsigned, 2> &N_par, val_T tol_p, const coord_vec_T &lattice_par, const vector<val_T> &U_par, const vector<val_T> &t_par, const vector<val_T> &onsite_E_par):SO_basis_Hamiltonian<ind_T, val_T, SD_T>(N_par, tol_p, lattice_par, U_par, t_par, onsite_E_par)
{
	ind_max=numeric_limits<ind_T>::max();
	create_pascal_matrix();
}

template<class ind_T, class val_T, class SD_T, class cfs_T>
void Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::compare_H_matrix()
{
	ind_T i,j;

	create_H_matrix_ph();
	H_m_2=H_m;
	
	auto it=H_m[0].begin();
	
	create_H_matrix();
//	create_H_matrix_abs();
	
	cout<<"H_m vs H_m_2:\n";
	for (i=0; i<NSD; i++)
	{
		for (auto &elem: H_m[i])
		{
			it=H_m_2[i].find(elem.first);
			if (it!=H_m_2[i].end())
			{
				//cout<<setw(10)<<i<<setw(10)<<elem.first<<setw(20)<<elem.second<<setw(20)<<it->second<<elem.second-it->second<<endl;
				print_SD(CI_ph_basis[i]);
				cout<<"\t\t";
				print_SD(CI_ph_basis[elem.first]);
				cout<<"\t\t"<<setw(20)<<elem.second<<setw(20)<<it->second<<elem.second-it->second<<endl;
			}
			else
			{
			//	cout<<setw(10)<<i<<setw(10)<<elem.first<<elem.second<<endl;
				print_SD(CI_ph_basis[i]);
				cout<<"\t\t";
				print_SD(CI_ph_basis[elem.first]);
				cout<<"\t\t"<<setw(20)<<elem.second<<endl;
			}
		}

	}
	
	
	cout<<"H_m_2 vs H_m:\n";
	for (i=0; i<NSD; i++)
	{
		for (auto &elem: H_m_2[i])
		{
			it=H_m[i].find(elem.first);
			if (it!=H_m[i].end())
			{
			//	cout<<setw(10)<<i<<setw(10)<<elem.first<<setw(20)<<elem.second<<setw(20)<<it->second<<elem.second-it->second<<endl;
				print_SD(CI_ph_basis[i]);
				cout<<"\t\t";
				print_SD(CI_ph_basis[elem.first]);
				cout<<"\t\t"<<setw(20)<<elem.second<<setw(20)<<it->second<<elem.second-it->second<<endl;
			}
			else
			{
			//	cout<<setw(10)<<i<<setw(10)<<elem.first<<elem.second<<endl;
				print_SD(CI_ph_basis[i]);
				cout<<"\t\t";
				print_SD(CI_ph_basis[elem.first]);
				cout<<"\t\t"<<setw(20)<<elem.second<<endl;
			}
		}
		
	}
	
//	cout<<setw(10)<<1;
//	print_SD(CI_ph_basis[1]);
//	cout<<"\t\t"<<setw(10)<<10;
//	print_SD(CI_ph_basis[10]);
//	cout<<endl;
	
}

template<class ind_T, class val_T, class SD_T, class cfs_T>
void Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::compute_natural_spin_orbitals()
{
	if (GS_vec_SO.size()!=NSD)
	{
		cerr<<"compute_natural_spin_orbitals() error: no ground state vector available\n";
		exit(EXIT_FAILURE);
	}
	
	int i,j,q;
	ind_T l, r;
	SD_T SD_l_abs, SD_l_ph;
	short int sgn_l, sgn_r;
	
	double sgn;
	cfs_T DM_tmp;
	
	DM_up.clear();
	DM_up.resize(L);
	DM_down.clear();
	DM_down.resize(L);
	for (i=0; i<L; i++)
	{
		DM_up[i].resize(L,0);
		DM_down[i].resize(L,0);
	}
	
	//spin-up off-diagonal elements
	for (i=0; i<L-1; i++)
	{
		for (j=i+1; j<L; j++)
		{
			DM_tmp=0;
			
			for (r=0; r<NSD; r++)
			{
				if (!tstbit(CI_abs_basis[r],i) && tstbit(CI_abs_basis[r],j))
				{
					sgn_r=abs_to_ph_sign[r];
					SD_l_abs=CI_abs_basis[r];
					q=add_bits_between(SD_l_abs,i,j);
					sgn=1;
					if (q%2) sgn=-1;
					clrbit(SD_l_abs,j);
					setbit(SD_l_abs,i);
					convert_to_ph_SD(SD_l_abs,SD_l_ph,sgn_l);
					l=get_SD_index(SD_l_ph);
					sgn=sgn*sgn_l*sgn_r;
					DM_tmp+=sgn*conj(GS_vec_SO[l])*GS_vec_SO[r];
				}
			}
			
			DM_up[i][j]=DM_tmp;
			DM_up[j][i]=conj(DM_tmp);
		}
	}
	
	//spin-up diagonal elements
	for (i=0; i<L; i++)
	{
		DM_tmp=0;
		
		for (r=0; r<NSD; r++)
		{
			if (tstbit(CI_abs_basis[r],i))
			{
				DM_tmp+=conj(GS_vec_SO[r])*GS_vec_SO[r];
			}
		}
		
		DM_up[i][i]=DM_tmp;
	}
	
	//spin-down off-diagonal elements
	for (i=L; i<2*L-1; i++)
	{
		for (j=i+1; j<2*L; j++)
		{
			DM_tmp=0;
			
			for (r=0; r<NSD; r++)
			{
				if (!tstbit(CI_abs_basis[r],i) && tstbit(CI_abs_basis[r],j))
				{
					sgn_r=abs_to_ph_sign[r];
					SD_l_abs=CI_abs_basis[r];
					q=add_bits_between(SD_l_abs,i,j);
					sgn=1;
					if (q%2) sgn=-1;
					clrbit(SD_l_abs,j);
					setbit(SD_l_abs,i);
					convert_to_ph_SD(SD_l_abs,SD_l_ph,sgn_l);
					l=get_SD_index(SD_l_ph);
					sgn=sgn*sgn_l*sgn_r;
					DM_tmp+=sgn*conj(GS_vec_SO[l])*GS_vec_SO[r];
				}
			}
			
			DM_down[i-L][j-L]=DM_tmp;
			DM_down[j-L][i-L]=conj(DM_tmp);
		}
	}
	
	//spin-down diagonal elements
	for (i=L; i<2*L; i++)
	{
		DM_tmp=0;
		
		for (r=0; r<NSD; r++)
		{
			if (tstbit(CI_abs_basis[r],i))
			{
				DM_tmp+=conj(GS_vec_SO[r])*GS_vec_SO[r];
			}
		}
		
		DM_down[i-L][i-L]=DM_tmp;
	}
	
//	cout<<"single-particle density matrix computed\n";
	
	cout<<"DM_up:\n";
	print_matrix(DM_up);
	
	cout<<"DM_down:\n";
	print_matrix(DM_down);
	
/*
 
 	operator_type oi, oj;
 
	//compute the single-particle density matrix
	//spin-up off-diagonal elements
	for (i=0; i<L-1; i++)
	{
		if (SO_ph_map[i]<Nu) oi=hole;
		else oi=part;
		
		for (j=i+1; j<L; j++)
		{
			DM_tmp=0;
			
			if (SO_ph_map[j]<Nu) oj=hole;
			else oj=part;
			
			if (oi==hole && oj==hole) //
			{
				for (r=0; r<NSD; r++)
				{
					if (tstbit(CI_ph_basis[r],SO_ph_map[i]) && !tstbit(CI_ph_basis[r],SO_ph_map[j]))
					{
						SD_tmp=CI_ph_basis[r];
						q=1;
						q+=add_bits_between(SD_tmp, SO_ph_map[i], SO_ph_map[j]);
						sgn=1;
						if (q%2) sgn=-1;
						clrbit(SD_tmp,SO_ph_map[i]);
						setbit(SD_tmp,SO_ph_map[j]);
						l=get_SD_index(SD_tmp);
						DM_tmp+=sgn*conj(GS_vec_SO[l])*GS_vec_SO[r];
					}
				}
			}
			else if (oi==part && oj==part)
			{
				for (r=0; r<NSD; r++)
				{
					if (tstbit(CI_ph_basis[r],SO_ph_map[j]) && !tstbit(CI_ph_basis[r],SO_ph_map[i]))
					{
						SD_tmp=CI_ph_basis[r];
						q=add_bits_between(SD_tmp, SO_ph_map[i], SO_ph_map[j]);
						sgn=1;
						if (q%2) sgn=-1;
						clrbit(SD_tmp,SO_ph_map[j]);
						setbit(SD_tmp,SO_ph_map[i]);
						l=get_SD_index(SD_tmp);
						DM_tmp+=sgn*conj(GS_vec_SO[l])*GS_vec_SO[r];
					}
				}
			}
			else if (oi==part && oj==hole)
			{
				for (r=0; r<NSD; r++)
				{
					if (!tstbit(CI_ph_basis[r],SO_ph_map[i]) && !tstbit(CI_ph_basis[r],SO_ph_map[j]))
					{
						SD_tmp=CI_ph_basis[r];
						q=add_bits_between(SD_tmp, SO_ph_map[i], SO_ph_map[j]);
						sgn=1;
						if (q%2) sgn=-1;
						setbit(SD_tmp,SO_ph_map[i]);
						setbit(SD_tmp,SO_ph_map[j]);
						l=get_SD_index(SD_tmp);
						DM_tmp+=sgn*conj(GS_vec_SO[l])*GS_vec_SO[r];
					}
				}
			}
			else
			{
				for (r=0; r<NSD; r++)
				{
					if (tstbit(CI_ph_basis[r],SO_ph_map[i]) && tstbit(CI_ph_basis[r],SO_ph_map[j]))
					{
						SD_tmp=CI_ph_basis[r];
						q=add_bits_between(SD_tmp, SO_ph_map[j], SO_ph_map[i]);
						sgn=1;
						if (q%2) sgn=-1;
						clrbit(SD_tmp,SO_ph_map[i]);
						clrbit(SD_tmp,SO_ph_map[j]);
						l=get_SD_index(SD_tmp);
						DM_tmp+=sgn*conj(GS_vec_SO[l])*GS_vec_SO[r];
					}
				}
			}
			
			DM_up[i][j]=DM_tmp;
			DM_up[j][i]=conj(DM_tmp);
		}
	}
	
	//spin-up diagonal elements
	for (i=0; i<L; i++)
	{
		if (SO_ph_map[i]<Nu) oi=hole;
		else oi=part;
		
		DM_tmp=0;
		
		for (r=0; r<NSD; r++)
		{
			if ((oi==part && tstbit(CI_ph_basis[r],SO_ph_map[i])) || (oi==hole && !tstbit(CI_ph_basis[r],SO_ph_map[i])))
			{
				DM_tmp+=conj(GS_vec_SO[r])*GS_vec_SO[r];
			}
		}
		
		DM_up[i][i]=DM_tmp;
	}
	
	
	//spin-down off-diagonal elements
	for (i=L; i<2*L-1; i++)
	{
		if ((SO_ph_map[i]-L)<Nd) oi=hole;
		else oi=part;
		
		for (j=i+1; j<2*L; j++)
		{
			DM_tmp=0;
			
			if ((SO_ph_map[j]-L)<Nd) oj=hole;
			else oj=part;
			
			if (oi==hole && oj==hole)
			{
				for (r=0; r<NSD; r++)
				{
					if (tstbit(CI_ph_basis[r],SO_ph_map[i]) && !tstbit(CI_ph_basis[r],SO_ph_map[j]))
					{
						SD_tmp=CI_ph_basis[r];
						q=1;
						q+=add_bits_between(SD_tmp, SO_ph_map[i], SO_ph_map[j]);
						sgn=1;
						if (q%2) sgn=-1;
						clrbit(SD_tmp,SO_ph_map[i]);
						setbit(SD_tmp,SO_ph_map[j]);
						l=get_SD_index(SD_tmp);
						DM_tmp+=sgn*conj(GS_vec_SO[l])*GS_vec_SO[r];
					}
				}
			}
			else if (oi==part && oj==part)
			{
				for (r=0; r<NSD; r++)
				{
					if (tstbit(CI_ph_basis[r],SO_ph_map[j]) && !tstbit(CI_ph_basis[r],SO_ph_map[i]))
					{
						SD_tmp=CI_ph_basis[r];
						q=add_bits_between(SD_tmp, SO_ph_map[i], SO_ph_map[j]);
						sgn=1;
						if (q%2) sgn=-1;
						clrbit(SD_tmp,SO_ph_map[j]);
						setbit(SD_tmp,SO_ph_map[i]);
						l=get_SD_index(SD_tmp);
						DM_tmp+=sgn*conj(GS_vec_SO[l])*GS_vec_SO[r];
					}
				}
			}
			else if (oi==part && oj==hole)
			{
				for (r=0; r<NSD; r++)
				{
					if (!tstbit(CI_ph_basis[r],SO_ph_map[i]) && !tstbit(CI_ph_basis[r],SO_ph_map[j]))
					{
						SD_tmp=CI_ph_basis[r];
						q=add_bits_between(SD_tmp, SO_ph_map[i], SO_ph_map[j]);
						sgn=1;
						if (q%2) sgn=-1;
						setbit(SD_tmp,SO_ph_map[i]);
						setbit(SD_tmp,SO_ph_map[j]);
						l=get_SD_index(SD_tmp);
						DM_tmp+=sgn*conj(GS_vec_SO[l])*GS_vec_SO[r];
					}
				}
			}
			else
			{
				for (r=0; r<NSD; r++)
				{
					if (tstbit(CI_ph_basis[r],SO_ph_map[i]) && tstbit(CI_ph_basis[r],SO_ph_map[j]))
					{
						SD_tmp=CI_ph_basis[r];
						q=add_bits_between(SD_tmp, SO_ph_map[j], SO_ph_map[i]);
						sgn=1;
						if (q%2) sgn=-1;
						clrbit(SD_tmp,SO_ph_map[i]);
						clrbit(SD_tmp,SO_ph_map[j]);
						l=get_SD_index(SD_tmp);
						DM_tmp+=sgn*conj(GS_vec_SO[l])*GS_vec_SO[r];
					}
				}
			}
			
			DM_down[i-L][j-L]=DM_tmp;
			DM_down[j-L][i-L]=conj(DM_tmp);
		}
	}
	
	//spin-down diagonal elements
	for (i=L; i<2*L; i++)
	{
		if ((SO_ph_map[i]-L)<Nd) oi=hole;
		else oi=part;
		
		DM_tmp=0;
		
		for (r=0; r<NSD; r++)
		{
			if ((oi==part && tstbit(CI_ph_basis[r],SO_ph_map[i])) || (oi==hole && !tstbit(CI_ph_basis[r],SO_ph_map[i])))
			{
				DM_tmp+=conj(GS_vec_SO[r])*GS_vec_SO[r];
			}
		}
		
		DM_down[i-L][i-L]=DM_tmp;
	}
*/
	
	//diagonalize the matrices
	char JOBZ[1];
	JOBZ[0]='V';
	char UPLO[1];
	UPLO[0]='U';
	int NA[1];
	NA[0]=L;
	int LDA[1];
	LDA[0]=L;
	double *W=new double[L];
	int LWORK[1];
	int INFO[1];
	
	graph_2D g1, g2;
	
	char xl_u[]="up natural spin-orbital index";
	char yl_u[]="up natural spin-orbital filling";
	char attr_u[]="'o',color='r',markerfacecolor='r'";
	char xl_d[]="down natural spin-orbital index";
	char yl_d[]="down natural spin-orbital filling";
	char attr_d[]="'o',color='b',markerfacecolor='b'";
	
	if (H_T::H_real)
	{
		natural_basis_up.clear();
		natural_basis_down.clear();
		
		double *A=new double[L*L];
		LWORK[0]=3*L-1;
		double *WORK=new double[LWORK[0]];
		
		//diagonalize DM_up
		for (i=0; i<L; i++)
			for (j=0; j<L; j++)
				A[i+j*L]=dval(real(DM_up[i][j]));
		
		dsyev_(JOBZ, UPLO, NA, A, LDA, W, WORK, LWORK, INFO);
		
		vector<double> SO_ind(L);
		for (i=0; i<L; i++) SO_ind[i]=(double)i;
		
		double ylims[]={0,1.0};
		
		if (!INFO[0])
		{
			vector<vector<val_T>> basis_up;
			
			DM_eig_up.resize(L);
			basis_up.resize(L);
			for (i=0; i<L; i++)
			{
				DM_eig_up[L-1-i]=W[i];
				basis_up[i].resize(L);
				for (j=0; j<L; j++)
					basis_up[i][L-1-j]=A[i+j*L];
			}
			
			g1.add_data(SO_ind.data(),DM_eig_up.data(),L);
			g1.add_attribute(attr_u);
			g1.set_axes_labels(xl_u,yl_u);
			g1.set_axes_lims(NULL,ylims);
			g1.curve_plot();
			
		//	matrix_prod(SO_H_T::SO_basis_u, basis_up, natural_basis_up);
			
			vector<vector<val_T>> Q;
			matrix_prod(SO_H_T::SO_basis_u, basis_up, Q);
			Gram_Schmidt_orthog(Q, natural_basis_up);
		}
		
		//diagonalize DM_down
		for (i=0; i<L; i++)
			for (j=0; j<L; j++)
				A[i+j*L]=dval(real(DM_down[i][j]));
		
		dsyev_(JOBZ, UPLO, NA, A, LDA, W, WORK, LWORK, INFO);
		
		if (!INFO[0])
		{
			vector<vector<val_T>> basis_down;
			
			DM_eig_down.resize(L);
			basis_down.resize(L);
			for (i=0; i<L; i++)
			{
				DM_eig_down[L-1-i]=W[i];
				basis_down[i].resize(L);
				for (j=0; j<L; j++)
					basis_down[i][L-1-j]=A[i+j*L];
			}
			
			g2.add_data(SO_ind.data(),DM_eig_down.data(),L);
			g2.add_attribute(attr_d);
			g2.set_axes_labels(xl_d,yl_d);
			g2.set_axes_lims(NULL,ylims);
			g2.curve_plot();
			
		//	matrix_prod(SO_H_T::SO_basis_d, basis_down, natural_basis_down);
			
			vector<vector<val_T>> Q;
			matrix_prod(SO_H_T::SO_basis_d, basis_down, Q);
			Gram_Schmidt_orthog(Q, natural_basis_down);
		}
		
		graph_2D::show_figures();
		
		if (A) delete [] A;
		if (WORK) delete [] WORK;
		
		if (W) delete [] W;
		
		unsigned ref_SD_tmp=0;
		
		for (i=0; i<Nu; i++) setbit(ref_SD_tmp, i);
		for (i=0; i<Nd; i++) setbit(ref_SD_tmp, i+L);
		
		cout<<"natural_basis_up:\n";
		print_matrix(natural_basis_up);
		
		cout<<"natural_basis_down:\n";
		print_matrix(natural_basis_down);
		
		SO_H_T::set_SO_basis(natural_basis_up,natural_basis_down,ref_SD_tmp);
	/*
		vector<vector<val_T>> NBut, NBdt, Id1;
		transpose_matrix(natural_basis_up, NBut);
		matrix_prod(natural_basis_up,NBut,Id1);
		cout<<"natural_basis_up*transpose(natural_basis_up):\n";
		print_matrix(Id1);
		
		transpose_matrix(natural_basis_down, NBdt);
		matrix_prod(natural_basis_down,NBdt,Id1);
		cout<<"natural_basis_down*transpose(natural_basis_down):\n";
		print_matrix(Id1);
	*/
	//	SO_H_T::transform_V(NBut,NBdt);
	//	SO_H_T::compare_V();
	
	}
}

template<class ind_T, class val_T, class SD_T, class cfs_T>
void Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::compute_ground_state(orbital_type_T OT, unsigned NL_max, bool display_figs, bool random_init_st)
{
	orbital_type=OT;
	if (orbital_type==local)
	{
		orbital_name="local";
	}
	else if (orbital_type==natural)
	{
		orbital_name="natural";
	}
	else
	{
		orbital_name="Hartree-Fock";
	}
		
	bool iterate_ref_SD=false;
	double tol_E_GS=1e-4;
	int Niter_ref_SD=1;
	
//	if (random_init_st) iterate_ref_SD=false;
	
//	create_H_matrix();
	create_H_matrix_ph();
	
	vector<double> eig_vals, eig_vals_all;
	
	unsigned NL, NL_all;
	unsigned NL_min=3;
	vector<double> NL_vec, NL_vec_all;
	
	graph_2D g1, g2, g3;
	
	char xl[]="$N_L$";
	char yl[]="E";
	char attr1[]="'o',color='r',markerfacecolor='r'";
	string orbitals=orbital_name+" orbitals";
	string ttl_E, ttl_cfs, ttl_cfs2;
	ttl_E="energies for "+orbitals;
	ttl_cfs="wave function coefficients for "+orbitals;
	ttl_cfs2="squared coefficients for "+orbitals;
	
	cout<<"computing ground state with "<<orbitals<<endl;
	
	if (!iterate_ref_SD) Niter_ref_SD=1;
	double D_E_GS, E_GS=0;
	double E_GS_prec=E_GS+2*tol_E_GS;
	
	ind_T i_max;
	cfs_T coeff2_max;
	vector<vector<val_T>> SO_basis_u_tmp, SO_basis_d_tmp;
	
	int i, iter_ref_SD=0;
	while (iter_ref_SD<Niter_ref_SD && fabs(E_GS-E_GS_prec)>tol_E_GS)
	{
		E_GS_prec=E_GS;
		
		NL_all=0;
		NL_vec.clear();
		NL_vec_all.clear();
		eig_vals.clear();
		eig_vals_all.clear();
		
		compute_GS_vec=false;
		Lanczos_init(NL_max, random_init_st);
		iL=1;
		
		//	cout<<setw(20)<<"NL"<<"E_min:"<<endl;
		while (iL<NL_max && Lanczos_iteration())
		{
			NL=iL+1;
			//	Lanczos(NL);
			
			if (NL>=NL_min)
			{
				eigs_Lanczos(eig_vals, false, NL);
				
				for (i=0; i<NL; i++) eig_vals[i]=eig_vals[i]+dval(E_ref);
				eig_vals_all.insert(eig_vals_all.end(),eig_vals.begin(),eig_vals.end());
				
				//		cout<<setw(20)<<NL<<eig_vals[0]<<endl;
				
				NL_vec.resize(NL);
				for (i=0; i<=iL; i++) NL_vec[i]=NL;
				NL_vec_all.insert(NL_vec_all.end(),NL_vec.begin(),NL_vec.end());
				
				NL_all+=NL;
			}
			
			iL++;
		}
		//	cout<<setw(20)<<NL<<eig_vals[0]<<endl;
		
		E_GS=eig_vals[0];
		
		cout<<"ground state energy:  "<<E_GS<<endl;
		
		if (display_figs)
		{
			g1.add_data(NL_vec_all.data(),eig_vals_all.data(),NL_all);
			g1.add_attribute(attr1);
			g1.set_axes_labels(xl,yl);
			g1.add_title(ttl_E.c_str());
			g1.curve_plot();
		//	graph_2D::show_figures();
		}
		
		compute_GS_vec=true;
		eigs_Lanczos(eig_vals, true, NL_max);
		Lanczos_init(NL_max, random_init_st);
		iL=1;
		while (iL<NL_max && Lanczos_iteration()) iL++;
		
		D_E_GS=eig_vals[0];
		E_GS=D_E_GS+dval(E_ref);
		
	//	cout<<"E_GS:  "<<E_GS<<endl;
	//	if (fabs(E_GS-E_GS_prec)<tol_E_GS) cout<<"gound state energy converged\n";
		
		/*
		 cfs_T vL_pr;
		 scalar_prod_vec<ind_T,cfs_T>(vL_pr,vL,vL0);
		 cout<<"vL*vL0: "<<vL_pr<<endl;
		 scalar_prod_vec<ind_T,cfs_T>(vL_pr,vL,vL1);
		 cout<<"vL*vL1: "<<vL_pr<<endl;
		 scalar_prod_vec<ind_T,cfs_T>(vL_pr,vL,vL2);
		 cout<<"vL*vL2: "<<vL_pr<<endl;
		 
		 scalar_prod_vec<ind_T,cfs_T>(vL_pr,GS_vec_SO,GS_vec_SO);
		 cout<<"GS_vec_SO*GS_vec_SO: "<<vL_pr<<endl;
		 
		 vL_pr=sqrt(vL_pr);
		 vec_div_equal<ind_T, cfs_T>(GS_vec_SO, vL_pr);
		 */
		//	cout<<"GS_vec_SO:\n";
		
		vector<cfs_T> H_GS_vec_SO;
		prod_mat_vec(H_GS_vec_SO, H_v, GS_vec_SO);
		
		i_max=0;
		coeff2_max=GS_vec_SO[0]*GS_vec_SO[0];
		
		vector<double> SO_ind(NSD), GS_vec_SO_d(NSD), GS_vec_SO_2(NSD), H_GS_vec_SO_d(NSD);
		for (ind_T j=0; j<NSD; j++)
		{
			SO_ind[j]=(double)j;
			GS_vec_SO_d[j]=dval(GS_vec_SO[j]);
			H_GS_vec_SO_d[j]=dval(H_GS_vec_SO[j])/D_E_GS;
			GS_vec_SO_2[j]=GS_vec_SO_d[j]*GS_vec_SO_d[j];
			//		cout<<GS_vec_SO[j]<<endl;
			if (GS_vec_SO_2[j]>coeff2_max)
			{
				i_max=j;
				coeff2_max=GS_vec_SO_2[j];
			}
		}
		
		if (i_max==0)
		{
			iter_ref_SD=Niter_ref_SD;
		//	cout<<"reference SD has the largest weight\n";
		//	cout<<"gound state energy converged\n";
		}
		
		char xl_GS[]="spin-orbital index";
		char yl_GS[]="GS coeff";
		char yl_GS_2[]="$(GS coeff)^2$";
		char attr_GS[]="'o',color='b',markerfacecolor='none'";
		char attr_H_GS[]="'s',color='r',markerfacecolor='none'";
		char attr_GS_2[]="'^',color='c',markerfacecolor='c'";
		
		if (display_figs)
		{
			g2.add_data(SO_ind.data(),GS_vec_SO_d.data(),NSD);
			g2.add_attribute(attr_GS);
			g2.add_data(SO_ind.data(),H_GS_vec_SO_d.data(),NSD);
			g2.add_attribute(attr_H_GS);
			g2.set_axes_labels(xl_GS,yl_GS);
			g2.add_title(ttl_cfs.c_str());
			g2.curve_plot();
			
			g3.add_data(SO_ind.data(),GS_vec_SO_2.data(),NSD);
			g3.add_attribute(attr_GS_2);
			g3.set_axes_labels(xl_GS,yl_GS_2);
			g3.add_title(ttl_cfs2.c_str());
			g3.curve_plot();
		}
		
		if (display_figs) graph_2D::show_figures();
		
		if (iterate_ref_SD && i_max!=0)
		{
			SO_basis_u_tmp=SO_H_T::SO_basis_u;
			SO_basis_d_tmp=SO_H_T::SO_basis_d;
			SO_H_T::set_SO_basis(SO_basis_u_tmp,SO_basis_d_tmp,CI_abs_basis[i_max]);
			create_CI_ph_basis();
		}
		
		iter_ref_SD++;
	}
/*
	NL_max=50;
	eig_vals.clear();
	compute_GS_vec=false;
	Lanczos_init(NL_max);
	iL=1;
	
	cout<<setw(20)<<"NL"<<"E_min:"<<endl;
	while (iL<NL_max && Lanczos_iteration())
	{
		NL=iL+1;
		//	Lanczos(NL);
		
		if (NL>=NL_min)
		{
			eigs_Lanczos(eig_vals, false, NL);
			
			for (i=0; i<NL; i++) eig_vals[i]=eig_vals[i]+E_ref;
			
			cout<<setw(20)<<NL<<eig_vals[0]<<endl;
			
		//	NL_vec.resize(NL);
		//	for (i=0; i<=iL; i++) NL_vec[i]=NL;
		
		//	g1.add_data(NL_vec.data(),eig_vals.data(),NL);
		//	g1.add_attribute(attr1);
		//	g1.set_axes_labels(xl,yl);
		//	g1.curve_plot();
		}
		
		iL++;
	}
	
	compute_GS_vec=true;
	eigs_Lanczos(eig_vals, true, NL_max);
	Lanczos_init(NL_max);
	iL=1;
	while (iL<NL_max && Lanczos_iteration()) iL++;
	
	//	cout<<"GS_vec_SO:\n";
	for (ind_T j=0; j<NSD; j++)
	{
		SO_ind[j]=(double)j;
		GS_vec_SO_2[j]=dval(GS_vec_SO[j]*GS_vec_SO[j]);
		//		cout<<GS_vec_SO[j]<<endl;
	}

	char attr_GS_1[]="'s',color='r'";
	char attr_GS_2_1[]="'v',color='m',markerfacecolor='m'";
	
	g2.add_data(SO_ind.data(),GS_vec_SO.data(),NSD);
	g2.add_attribute(attr_GS_1);
	g2.set_axes_labels(xl_GS,yl_GS);
	g2.curve_plot();
	
	g3.add_data(SO_ind.data(),GS_vec_SO_2.data(),NSD);
	g3.add_attribute(attr_GS_2_1);
	g3.set_axes_labels(xl_GS,yl_GS_2);
	g3.curve_plot();
*/
	
}

template<class ind_T, class val_T, class SD_T, class cfs_T>
void Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::eigs_Lanczos(vector<double> &eig_vals, bool compute_eig_vectors, unsigned NL)
{
//	double told=1.0e-8, told1=10*EPSILON;
//	double tolmp0=1.0e-5, tolmp1=1.0e-100;
//	cfs_T tol_z, z0, zx;
	
//	int NL=aL.size();
	
	double *a=new double[NL];
	double *b=new double[NL-1];
	
	cfs_T vtmp;
	
//	cfs_T *A=new cfs_T[NL];
//	cfs_T *B2=new cfs_T[NL-1];
//	cfs_T roots[2];
	
	a[0]=dval(aL[0]);
//	A[0]=aLp[0];
	
	int i, j;
	for (j=1; j<NL; j++)
	{
		a[j]=dval(aL[j]);
		b[j-1]=dval(bL[j]);
		//vtmp=sqrt(b2L[j]);
		//b[j-1]=dval(vtmp);

	}
	
	int INFO[1];
	int NLp[1];
	NLp[0]=NL;
	
	if (!compute_eig_vectors)
		dsterf_(NLp, a, b, INFO);
	else
	{
		double *V=new double[NL*NL];
		double *work=new double[2*(NL-1)];
		char compz[1];
		compz[0]='I';
		int ldz[1];
		ldz[0]=NL;
		
		dsteqr_(compz, NLp, a, b, V, ldz, work, INFO);
		
	//	cout<<"eigs_Lanczos() INFO: "<<INFO[0]<<endl;
	//	cout<<"eigs_Lanczos() a[0]+E_ref: "<<a[0]+E_ref<<endl;
		
	//	cout<<"GS_vec_Lanczos:\n";
		if (!INFO[0])
		{
			GS_vec_Lanczos.clear();
			GS_vec_Lanczos.resize(NL);
			for (i=0; i<NL; i++)
			{
				GS_vec_Lanczos[i]=V[i];
			//	cout<<V[i]<<endl;
			}
		}
		
		if (work) delete [] work;
		if (V) delete [] V;
	}
	
	eig_vals.resize(NL);
	for (i=0; i<NL; i++)
	{
		eig_vals[i]=a[i];
	}
	
	if (INFO[0])
	{
		cout<<"INFO: "<<INFO[0]<<endl;
		cout<<"NL: "<<NL<<endl;
	}
	
	if (a) delete [] a;
	if (b) delete [] b;
 
/*
	cfs_T Dz, z_av, z1, z2, z3;
	
	//	cout<<setprecision(15)<<setiosflags(ios::left);

	bool root_recomputed=false;
	
	cfs_T dr=a[NHL-1]-a[0];
	cfs_T rtmp;
	
	
	cout<<setprecision(15);
	if (!INFO[0])
	{
		for (j=0; j<NHL; j++)
		{
			if (!root_recomputed)
				eigs_H[j]=a[j];
			else
			{
				eigs_H[j]=roots[1];
				root_recomputed=false;
			}
			if (j<NHL-1)
			{
				z1=eigs_H[j];
				z2=a[j+1];
				//	z_av=0.5*(abs(z1)+abs(z2));
				//cout<<setw(20)<<z1<<setw(20)<<z2<<setw(20)<<z_av<<z2-z1<<endl;
				if ((z2-z1)<told*dr)
				{
					cout<<"quasi-degeneracy:\n";
					cout<<setw(20)<<"z_j"<<setw(20)<<"z_{j+1}"<<"z_{j+1}-z_j\n";
					cout<<setw(20)<<z1<<setw(20)<<z2<<z2-z1<<endl;
					cout<<"attempting to refine precision\n";
					if (j>0)
					{
						Dz=eigs_H[j]-eigs_H[j-1];
						if (j<NHL-2)
						{
							z3=a[j+2];
							if ((z3-z2)<Dz) Dz=z3-z2;
						}
					}
					if ((z2-z1)>told1*dr)
					{
						tol_z=tolmp0*(z2-z1);
					}
					else
						tol_z=tolmp1*Dz;
					z0=0.5*(z1+z2);
					extremum_det_H(z0, A, B2, NHL, tol_z, zx);
					roots_det_H(zx, A, B2, NHL, tol_z, roots);
					cout<<"refined roots:\n";
					if (roots[0]>roots[1])
					{
						rtmp=roots[0];
						roots[0]=roots[1];
						roots[1]=rtmp;
					}
					cout<<setw(20)<<roots[0]<<setw(20)<<roots[1]<<roots[1]-roots[0]<<endl;
					eigs_H[j]=roots[0];
					root_recomputed=true;
				}
			}
		}
		*/
		/*
		 cout<<"eigenvalues of HL:\n";
		 for (j=0; j<NHL; j++)
		 {
		 cout<<eigs_H[j]+Ndo_init*U<<endl;
		 }
		 *//*
	}
	cout<<setprecision(7);
	*/
	
//	if (A) delete [] A;
//	if (B2) delete [] B2;


}

template<class ind_T, class val_T, class SD_T, class cfs_T>
void Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::test_Lanczos(unsigned Ntst_p, unsigned NL_max, val_T a_P, val_T a_eig)
{
	NSD=Ntst_p;
	
	cout<<"NSD: "<<NSD<<endl;
	
	mt19937 rnd_gen;
	uniform_int_distribution<unsigned> distr;
	uniform_real_distribution<double> distr_real(0,1.0);
	
	rnd_gen.seed(time(NULL));
	
	val_T eig_tmp;
	set<val_T> eig_v;
	auto it=eig_v.begin();
	
	vector<vector<val_T>> P,Pt,Q, H_tmp;
	Q.resize(NSD);
	int i,j;
	for (i=0; i<NSD; i++)
	{
		eig_tmp=a_eig*(2*distr_real(rnd_gen)-1.0);
		eig_v.insert(eig_tmp);
		Q[i].resize(NSD);
		for (j=0; j<NSD; j++)
		{
			Q[i][j]=a_P*(2*distr_real(rnd_gen)-1.0);
		}
	}
	
	it=eig_v.begin();
	cout<<"eigs[0]: "<<*it<<endl;

	Gram_Schmidt_orthog(Q,P);

	transpose_matrix(P,Pt);
	i=0;
	for (it=eig_v.begin(); it!=eig_v.end(); it++)
	{
//		cout<<i<<endl;
		eig_tmp=*it;
		vec_times_equal<ind_T, val_T>(Pt[i], eig_tmp);
	//	cout<<eig_tmp<<endl;
		i++;
	}

	matrix_prod(P,Pt,H_tmp);
	
	H_v.resize(NSD);
	
	for (i=0; i<NSD; i++)
	{
		H_v[i].resize(NSD);
		for (j=0; j<NSD; j++)
		{
			H_v[i][j].index=j;
			H_v[i][j].value=H_tmp[i][j];
		}
	}
	
	cout<<"test_Lanczos(): Hamiltonian matrix computed\n";
	
	
	int NL, NL_min=3;
//	Lanczos_init(NL_max,true);
	Lanczos_init(NL_max,false);
	iL=1;
	while (iL<NL_max && Lanczos_iteration())
	{
		NL=iL+1;
		
		if (NL>NL_min)
		{
			double *a=new double[NL];
			double *b=new double[NL-1];
			
			cfs_T vtmp;
			
			a[0]=dval(aL[0]);
			
			int i, j;
			for (j=1; j<NL; j++)
			{
				a[j]=dval(aL[j]);
				b[j-1]=dval(bL[j]);
				//vtmp=sqrt(b2L[j]);
				//b[j-1]=dval(vtmp);
				
			}
			
			int INFO[1];
			int NLp[1];
			NLp[0]=NL;
			
			dsterf_(NLp, a, b, INFO);
			
			cout<<setw(10)<<NL<<a[0]<<endl;
			
			if (INFO[0])
			{
				cout<<"INFO: "<<INFO[0]<<endl;
				cout<<"NL: "<<NL<<endl;
			}
			
			if (a) delete [] a;
			if (b) delete [] b;
			
		}
		iL++;
	}
	
	/*
	
	H_tst.resize(Ntst);

	for (i=0; i<Ntst; i++)
	{
		H_tst[i].resize(Ntst);
		for (j=0; j<Ntst; j++)
		{
			H_tst[i][j].index=j;
			H_tst[i][j].value=H_tmp[i][j];
		}
	}
	
	cout<<"test_Lanczos(): Hamiltonian matrix computed\n";
	
	
	int NL, NL_min=3;
	Lanczos_init_tst(NL_max,a_P);
	iL=1;
	while (iL<NL_max && Lanczos_iteration_tst())
	{
		NL=iL+1;
		
		if (NL>NL_min)
		{
			double *a=new double[NL];
			double *b=new double[NL-1];
			
			cfs_T vtmp;
			
			a[0]=dval(aL_tst[0]);
			
			int i, j;
			for (j=1; j<NL; j++)
			{
				a[j]=dval(aL_tst[j]);
				b[j-1]=dval(bL_tst[j]);
				//vtmp=sqrt(b2L[j]);
				//b[j-1]=dval(vtmp);
				
			}
			
			int INFO[1];
			int NLp[1];
			NLp[0]=NL;
			
			
			dsterf_(NLp, a, b, INFO);
			
			cout<<setw(10)<<NL<<a[0]<<endl;
			
			if (INFO[0])
			{
				cout<<"INFO: "<<INFO[0]<<endl;
				cout<<"NL: "<<NL<<endl;
			}
			
			if (a) delete [] a;
			if (b) delete [] b;
			
		}
		iL++;
	}
	
	*/
}

/*
template<class ind_T, class val_T, class SD_T, class cfs_T>
void Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::Lanczos_init_tst(unsigned NL, cfs_T a_vL)
{
	cfs_T vHv;
	
	aL_tst.clear();
	aL_tst.resize(NL,0);
	bL_tst.clear();
	bL_tst.resize(NL,0);
	
	vL_tst.clear();
	vL_tst.resize(Ntst,0);
	vLm_tst=vL_tst;
	
	//i=0
*/
	/*
	 mt19937 rnd_gen;
	 uniform_int_distribution<unsigned> distr;
	 uniform_real_distribution<double> distr_real(0,1.0);
	for (int i=0; i<Ntst; i++)
	{
		vL_tst[i]=a_vL*(2*distr_real(rnd_gen)-1.0);
	}
	scalar_prod_vec<ind_T, cfs_T>(mod2v, vL_tst, vL_tst);
	vec_times_equal<ind_T, cfs_T>(vL_tst,1.0/sqrt(mod2v));
	 */
/*
	vL_tst[0]=1;
	
	//keep first Lanczos to check orthogonality
//	vL0=vL;
	
	//a[iL]=<vL|H|vL>
	prod_mat_vec(vLp_tst, H_tst, vL_tst);
	scalar_prod_vec<ind_T, cfs_T>(vHv, vL_tst, vLp_tst);
	aL_tst[0]=real(vHv);
	//iL=0: vLp=H*vL-a[0]*vL
	vec_plus_equal<ind_T, cfs_T>(vLp_tst, vL_tst, -aL_tst[0]);
	

}

template<class ind_T, class val_T, class SD_T, class cfs_T>
bool Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::Lanczos_iteration_tst()
{
	cfs_T vHv;
	
	vLm_tst=vL_tst;
	vL_tst=vLp_tst;
	
	//mod2v=(vL^+)*vL
	scalar_prod_vec<ind_T, cfs_T>(mod2v, vL_tst, vL_tst);
	
	if (mod2v)
	{
		//bL[iL]=sqrt((vL^+)*vL)
		bL_tst[iL]=sqrt(real(mod2v));
		//vL=vL/sqrt((vL^+)*vL)
		vec_times_equal<ind_T, cfs_T>(vL_tst,1.0/bL_tst[iL]);
		
	//	if (iL==1) vL1=vL;
	//	else if (iL==2) vL2=vL;
 
		//vLp=H*vL
		prod_mat_vec(vLp_tst, H_tst, vL_tst);
		//aL[iL]=(vL^+)*H*vL
		scalar_prod_vec<ind_T, cfs_T>(vHv, vL_tst, vLp_tst);
		aL_tst[iL]=real(vHv);
		//vLp=H*vL-aL[iL]*vL-bL[iL]*vLm
		vec_plus_equal<ind_T, cfs_T>(vLp_tst, vL_tst, -aL_tst[iL]);
		vec_plus_equal<ind_T, cfs_T>(vLp_tst, vLm_tst, -bL_tst[iL]);
		
	//	cout<<setw(10)<<iL<<setw(20)<<aL_tst[iL]<<bL_tst[iL]<<endl;
		return true;
	}
	else
		return false;
	
}
*/

template<class ind_T, class val_T, class SD_T, class cfs_T>
void Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::Lanczos_init(unsigned NL, bool random_init_st)
{
	cfs_T vHv;
	
	aL.clear();
	aL.resize(NL,0);
	bL.clear();
	bL.resize(NL,0);
	
	vL.clear();
	vL.resize(NSD,0);
	vLm=vL;
	
	//i=0
	if (compute_GS_vec)
	{
		vL=vL_init;
	}
	else if (random_init_st)
	{
	//	cout<<"Lanczos_init(): using random initial state.\n";
		mt19937 rnd_gen;
		uniform_int_distribution<unsigned> distr;
		uniform_real_distribution<double> distr_real(0,1.0);
		
		for (int i=0; i<NSD; i++)
		{
			vL[i]=2*distr_real(rnd_gen)-1.0;
		}
		scalar_prod_vec<ind_T, cfs_T>(mod2v, vL, vL);
		vec_times_equal<ind_T, cfs_T>(vL,1.0/sqrt(mod2v));
		vL_init=vL;
	}
	else
	{
	//	cout<<"Lanczos_init(): using reference SD as initial state.\n";
		vL[0]=1;
		vL_init=vL;
	}
	
	//keep first Lanczos to check orthogonality
//	vL0=vL;
	
	//vLp=H*vL-a[0]*vL
	prod_mat_vec(vLp, H_v, vL);
	scalar_prod_vec<ind_T, cfs_T>(vHv, vL, vLp);
	aL[0]=real(vHv);
	//iL=0: vLp=H*vL-a[0]*vL
	vec_plus_equal<ind_T, cfs_T>(vLp, vL, -aL[0]);
	if (compute_GS_vec)
	{
		GS_vec_SO.clear();
		GS_vec_SO.resize(NSD,0);
		vec_plus_equal<ind_T, cfs_T>(GS_vec_SO, vL, GS_vec_Lanczos[0]);
		//	cout<<"GS_vec_SO:\n";
		//	for (ind_T j=0; j<NSD; j++) cout<<setw(10)<<j<<GS_vec_SO[j]<<endl;
	}
}

template<class ind_T, class val_T, class SD_T, class cfs_T>
bool Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::Lanczos_iteration()
{
	cfs_T vHv;
	
	vLm=vL;
	vL=vLp;
	
	//mod2v=(vL^+)*vL
	scalar_prod_vec<ind_T, cfs_T>(mod2v, vL, vL);
	
	if (mod2v)
	{
		//bL[iL]=sqrt((vL^+)*vL)
		bL[iL]=sqrt(real(mod2v));
		//vL=vL/sqrt((vL^+)*vL)
		vec_times_equal<ind_T, cfs_T>(vL,1.0/bL[iL]);
		
	//	if (iL==1) vL1=vL;
	//	else if (iL==2) vL2=vL;
		
		if (compute_GS_vec)
		{
		//	cout<<"iL: "<<iL<<endl;
		//	cout<<"GS_vec_Lanczos[iL]: "<<GS_vec_Lanczos[iL]<<endl;
		//	cout<<"mod2v: "<<mod2v<<endl;
			//	for (ind_T j=0; j<NSD; j++) cout<<setw(10)<<j<<GS_vec_SO[j]<<endl;
			vec_plus_equal<ind_T, cfs_T>(GS_vec_SO, vL, GS_vec_Lanczos[iL]);
		}
		
		//vLp=H*vL
		prod_mat_vec(vLp, H_v, vL);
		//aL[iL]=(vL^+)*H*vL
		scalar_prod_vec<ind_T, cfs_T>(vHv, vL, vLp);
		aL[iL]=real(vHv);
		//vLp=H*vL-aL[iL]*vL-bL[iL]*vLm
		vec_plus_equal<ind_T, cfs_T>(vLp, vL, -aL[iL]);
		vec_plus_equal<ind_T, cfs_T>(vLp, vLm, -bL[iL]);
		
	//	cout<<setw(10)<<iL<<setw(20)<<aL[iL]<<bL[iL]<<endl;
		return true;
	}
	else
		return false;
	
}

template<class ind_T, class val_T, class SD_T, class cfs_T>
void Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::Lanczos(unsigned NL)
{
	Lanczos_init(NL);
	
	iL=1;
	while (iL<NL && Lanczos_iteration()) iL++;
}

//version without the sign of the absolute representation
template<class ind_T, class val_T, class SD_T, class cfs_T>
void Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::create_H_matrix_abs()
{
	cout<<"Defining Hamiltonian matrix...\n";
	
	compute_ref_diag_E();
	
	H_m.clear();
	H_m.resize(NSD);
	
	SD_T SD_r_abs, SD_l_abs, SD_l_ph;
	ind_T i, ind_V;
	double sgn;
	short int sgn_r, sgn_l;
	int k,l,m,n;
	int q;
	val_T val;
	auto V_it=SO_H_T::V_u_SO_basis.begin();
	
	for (i=0; i<NSD-1; i++)
	{
		SD_r_abs=CI_abs_basis[i];
	//	sgn_r=abs_to_ph_sign[i];
		
		add_H_m_element_diag(i, -E_ref);
		
		//spin-up hopping
		for (l=0; l<L; l++)
		{
			if (tstbit(SD_r_abs,l))
			{
				//diagonal element
				val=SO_H_T::t_SO_basis_u[l][l];
				add_H_m_element_diag(i, val);
				//off-diagonal elements
				for (k=0; k<L; k++)
				{
					if (!tstbit(SD_r_abs,k))
					{
						SD_l_abs=SD_r_abs;
						q=add_bits_between(SD_r_abs,k,l);
						sgn=1;
						if (q%2) sgn=-1;
						clrbit(SD_l_abs,l);
						setbit(SD_l_abs,k);
						convert_to_ph_SD(SD_l_abs,SD_l_ph,sgn_l);
					//	sgn=sgn*sgn_l*sgn_r;
						val=sgn*SO_H_T::t_SO_basis_u[k][l];
						add_H_m_element_sym(i, SD_l_ph, val);
					}
				}
			}
		}
		//spin-down hopping
		for (l=L; l<2*L; l++)
		{
			if (tstbit(SD_r_abs,l))
			{
				//diagonal element
				val=SO_H_T::t_SO_basis_d[l-L][l-L];
				add_H_m_element_diag(i, val);
				//off-diagonal elements
				for (k=L; k<2*L; k++)
				{
					if (!tstbit(SD_r_abs,k))
					{
						SD_l_abs=SD_r_abs;
						q=add_bits_between(SD_r_abs,k,l);
						sgn=1;
						if (q%2) sgn=-1;
						clrbit(SD_l_abs,l);
						setbit(SD_l_abs,k);
						convert_to_ph_SD(SD_l_abs,SD_l_ph,sgn_l);
					//	sgn=sgn*sgn_l*sgn_r;
						val=sgn*SO_H_T::t_SO_basis_d[k-L][l-L];
						add_H_m_element_sym(i, SD_l_ph, val);
					}
				}
			}
		}
		
		//up-up interaction if finite
		if (SO_H_T::V_u_SO_basis.size())
		{
			for (n=0; n<L-1; n++)
			{
				if (tstbit(SD_r_abs,n))
				{
					for (m=n+1; m<L; m++)
					{
						if (tstbit(SD_r_abs,m))
						{
							for (l=0; l<L-1; l++)
							{
								if (!tstbit(SD_r_abs,l) || l==m || l==n)
								{
									for (k=l+1; k<L; k++)
									{
										if (!tstbit(SD_r_abs,k) || k==n || k==m)
										{
											ind_V=k+n*s[1]+l*s[2]+m*s[3];
											V_it=SO_H_T::V_u_SO_basis.find(ind_V);
											if (V_it!=SO_H_T::V_u_SO_basis.end())
											{
												SD_l_abs=SD_r_abs;
												if (n!=l)
												{
													q=0;
													if (l!=m)
													{
														q+=add_bits_between(SD_l_abs,l,m);
														clrbit(SD_l_abs,m);
														setbit(SD_l_abs,l);
													}
													if (k!=n)
													{
														q+=add_bits_between(SD_l_abs,k,n);
														clrbit(SD_l_abs,n);
														setbit(SD_l_abs,k);
													}
													if (SD_l_abs!=SD_r_abs)
													{
														sgn=1;
														if (q%2) sgn=-1;
														convert_to_ph_SD(SD_l_abs,SD_l_ph,sgn_l);
														//sgn=sgn*sgn_l*sgn_r;
														val=sgn*V_it->second;
														add_H_m_element_sym(i, SD_l_ph, val);
													}
													else
													{
														val=V_it->second;
														add_H_m_element_diag(i, val);
													}
												}
												else
												{
													q=1;
													if (k!=m)
													{
														q+=add_bits_between(SD_l_abs,k,m);
														clrbit(SD_l_abs,m);
														setbit(SD_l_abs,k);
													}
													sgn=1;
													if (q%2) sgn=-1;
													if (SD_l_abs!=SD_r_abs)
													{
														convert_to_ph_SD(SD_l_abs,SD_l_ph,sgn_l);
														//sgn=sgn*sgn_l*sgn_r;
														val=sgn*V_it->second;
														add_H_m_element_sym(i, SD_l_ph, val);
													}
													else
													{
														val=sgn*V_it->second;
														add_H_m_element_diag(i, val);
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		//down-down interaction if finite
		if (SO_H_T::V_d_SO_basis.size())
		{
			for (n=L; n<2*L-1; n++)
			{
				if (tstbit(SD_r_abs,n))
				{
					for (m=n+1; m<2*L; m++) //include only distinct pairs (m,n)
					{
						if (tstbit(SD_r_abs,m))
						{
							for (l=L; l<2*L-1; l++)
							{
								if (!tstbit(SD_r_abs,l) || l==m || l==n)
								{
									for (k=l+1; k<2*L; k++)
									{
										if (!tstbit(SD_r_abs,k) || k==n || k==m)
										{
											ind_V=k+n*s[1]+l*s[2]+m*s[3];
											V_it=SO_H_T::V_d_SO_basis.find(ind_V);
											if (V_it!=SO_H_T::V_d_SO_basis.end())
											{
												SD_l_abs=SD_r_abs;
												if (n!=l)
												{
													q=0;
													if (l!=m)
													{
														q+=add_bits_between(SD_l_abs,l,m);
														clrbit(SD_l_abs,m);
														setbit(SD_l_abs,l);
													}
													if (k!=n)
													{
														q+=add_bits_between(SD_l_abs,k,n);
														clrbit(SD_l_abs,n);
														setbit(SD_l_abs,k);
													}
													if (SD_l_abs!=SD_r_abs)
													{
														sgn=1;
														if (q%2) sgn=-1;
														convert_to_ph_SD(SD_l_abs,SD_l_ph,sgn_l);
														//sgn=sgn*sgn_l*sgn_r;
														val=sgn*V_it->second;
														add_H_m_element_sym(i, SD_l_ph, val);
													}
													else
													{
														val=V_it->second;
														add_H_m_element_diag(i, val);
													}
												}
												else
												{
													q=1;
													if (k!=m)
													{
														q+=add_bits_between(SD_l_abs,k,m);
														clrbit(SD_l_abs,m);
														setbit(SD_l_abs,k);
													}
													sgn=1;
													if (q%2) sgn=-1;
													if (SD_l_abs!=SD_r_abs)
													{
														convert_to_ph_SD(SD_l_abs,SD_l_ph,sgn_l);
														//sgn=sgn*sgn_l*sgn_r;
														val=sgn*V_it->second;
														add_H_m_element_sym(i, SD_l_ph, val);
													}
													else
													{
														val=sgn*V_it->second;
														add_H_m_element_diag(i, val);
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		//up-down interaction
		for (n=0; n<L; n++)
		{
			if (tstbit(SD_r_abs,n))
			{
				for (m=L; m<2*L; m++)
				{
					if (tstbit(SD_r_abs,m))
					{
						for (l=L; l<2*L; l++)
						{
							if (!tstbit(SD_r_abs,l) || l==m)
							{
								for (k=0; k<L; k++)
								{
									if (!tstbit(SD_r_abs,k) || k==n)
									{
										ind_V=k+n*s[1]+l*s[2]+m*s[3];
										V_it=SO_H_T::V_ud_SO_basis.find(ind_V);
										if (V_it!=SO_H_T::V_ud_SO_basis.end())
										{
											SD_l_abs=SD_r_abs;
											
											q=0;
											if (l!=m)
											{
												q+=add_bits_between(SD_l_abs,l,m);
												clrbit(SD_l_abs,m);
												setbit(SD_l_abs,l);
											}
											if (k!=n)
											{
												q+=add_bits_between(SD_l_abs,k,n);
												clrbit(SD_l_abs,n);
												setbit(SD_l_abs,k);
											}
											if (SD_l_abs!=SD_r_abs)
											{
												sgn=1;
												if (q%2) sgn=-1;
												convert_to_ph_SD(SD_l_abs,SD_l_ph,sgn_l);
												//sgn=sgn*sgn_l*sgn_r;
												val=sgn*V_it->second;
												add_H_m_element_sym(i, SD_l_ph, val);
											}
											else
											{
												val=V_it->second;
												add_H_m_element_diag(i, val);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
	}
	
	
	i=NSD-1;
	SD_r_abs=CI_abs_basis[i];
	
	add_H_m_element_diag(i, -E_ref);
	
	//spin-up hopping
	for (l=0; l<L; l++)
	{
		if (tstbit(SD_r_abs,l))
		{
			//diagonal element
			val=SO_H_T::t_SO_basis_u[l][l];
			add_H_m_element_diag(i, val);
		}
	}
	//spin-down hopping
	for (l=L; l<2*L; l++)
	{
		if (tstbit(SD_r_abs,l))
		{
			//diagonal element
			val=SO_H_T::t_SO_basis_d[l-L][l-L];
			add_H_m_element_diag(i, val);
		}
	}
	
	//up-up interaction if finite
	if (SO_H_T::V_u_SO_basis.size())
	{
		for (n=0; n<L-1; n++)
		{
			if (tstbit(SD_r_abs,n))
			{
				for (m=n+1; m<L; m++)
				{
					if (tstbit(SD_r_abs,m))
					{
						ind_V=n+n*s[1]+m*s[2]+m*s[3];
						V_it=SO_H_T::V_u_SO_basis.find(ind_V);
						if (V_it!=SO_H_T::V_u_SO_basis.end())
						{
							val=V_it->second;
							add_H_m_element_diag(i, val);
						}
					}
				}
			}
		}
	}
	
	
	//down-down interaction if finite
	if (SO_H_T::V_d_SO_basis.size())
	{
		for (n=L; n<2*L-1; n++)
		{
			if (tstbit(SD_r_abs,n))
			{
				for (m=n+1; m<2*L; m++)
				{
					if (tstbit(SD_r_abs,m))
					{
						ind_V=n+n*s[1]+m*s[2]+m*s[3];
						V_it=SO_H_T::V_d_SO_basis.find(ind_V);
						if (V_it!=SO_H_T::V_d_SO_basis.end())
						{
							val=V_it->second;
							add_H_m_element_diag(i, val);
						}
					}
				}
			}
		}
	}
	
	//up-down interaction
	for (n=0; n<L; n++)
	{
		if (tstbit(SD_r_abs,n))
		{
			for (m=L; m<2*L; m++)
			{
				if (tstbit(SD_r_abs,m))
				{
					ind_V=n+n*s[1]+m*s[2]+m*s[3];
					V_it=SO_H_T::V_ud_SO_basis.find(ind_V);
					if (V_it!=SO_H_T::V_ud_SO_basis.end())
					{
						val=V_it->second;
						add_H_m_element_diag(i, val);
					}
				}
			}
		}
	}
	
	cout<<"Hamiltonian matrix defined.\n";
	
	//	check_H_Hermicity();
	convert_H_m_to_H_v();
	
	cout<<"H[0][0]:\n";
	cout<<"index: "<<H_v[0][0].index<<endl;
	cout<<"value: "<<H_v[0][0].value<<endl;
	
}


//version with the sign of the particle-hole representation
template<class ind_T, class val_T, class SD_T, class cfs_T>
void Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::create_H_matrix()
{
	cout<<"Defining Hamiltonian matrix...\n";
	
	compute_ref_diag_E();
	
	H_m.clear();
	H_m.resize(NSD);
	
	SD_T SD_r_abs, SD_l_abs, SD_l_ph;
	ind_T i, ind_V;
	double sgn;
	short int sgn_r, sgn_l;
	int k,l,m,n;
	int q;
	val_T val;
	auto V_it=SO_H_T::V_u_SO_basis.begin();
	
	for (i=0; i<NSD-1; i++)
	{
		SD_r_abs=CI_abs_basis[i];
		sgn_r=abs_to_ph_sign[i];
		
		add_H_m_element_diag(i, -E_ref);
		
		//spin-up hopping
		for (l=0; l<L; l++)
		{
			if (tstbit(SD_r_abs,l))
			{
				//diagonal element
				val=SO_H_T::t_SO_basis_u[l][l];
				add_H_m_element_diag(i, val);
				//off-diagonal elements
				for (k=0; k<L; k++)
				{
					if (!tstbit(SD_r_abs,k))
					{
						SD_l_abs=SD_r_abs;
						q=add_bits_between(SD_r_abs,k,l);
						sgn=1;
						if (q%2) sgn=-1;
						clrbit(SD_l_abs,l);
						setbit(SD_l_abs,k);
						convert_to_ph_SD(SD_l_abs,SD_l_ph,sgn_l);
						sgn=sgn*sgn_l*sgn_r;
						val=sgn*SO_H_T::t_SO_basis_u[k][l];
						add_H_m_element_sym(i, SD_l_ph, val);
					}
				}
			}
		}
		//spin-down hopping
		for (l=L; l<2*L; l++)
		{
			if (tstbit(SD_r_abs,l))
			{
				//diagonal element
				val=SO_H_T::t_SO_basis_d[l-L][l-L];
				add_H_m_element_diag(i, val);
				//off-diagonal elements
				for (k=L; k<2*L; k++)
				{
					if (!tstbit(SD_r_abs,k))
					{
						SD_l_abs=SD_r_abs;
						q=add_bits_between(SD_r_abs,k,l);
						sgn=1;
						if (q%2) sgn=-1;
						clrbit(SD_l_abs,l);
						setbit(SD_l_abs,k);
						convert_to_ph_SD(SD_l_abs,SD_l_ph,sgn_l);
						sgn=sgn*sgn_l*sgn_r;
						val=sgn*SO_H_T::t_SO_basis_d[k-L][l-L];
						add_H_m_element_sym(i, SD_l_ph, val);
					}
				}
			}
		}
		
		//up-up interaction if finite
		if (SO_H_T::V_u_SO_basis.size())
		{
			for (n=0; n<L-1; n++)
			{
				if (tstbit(SD_r_abs,n))
				{
					for (m=n+1; m<L; m++)
					{
						if (tstbit(SD_r_abs,m))
						{
							for (l=0; l<L-1; l++)
							{
								if (!tstbit(SD_r_abs,l) || l==m || l==n)
								{
									for (k=l+1; k<L; k++)
									{
										if (!tstbit(SD_r_abs,k) || k==n || k==m)
										{
											ind_V=k+n*s[1]+l*s[2]+m*s[3];
											V_it=SO_H_T::V_u_SO_basis.find(ind_V);
											if (V_it!=SO_H_T::V_u_SO_basis.end())
											{
												SD_l_abs=SD_r_abs;
												if (n!=l)
												{
													q=0;
													if (l!=m)
													{
														q+=add_bits_between(SD_l_abs,l,m);
														clrbit(SD_l_abs,m);
														setbit(SD_l_abs,l);
													}
													if (k!=n)
													{
														q+=add_bits_between(SD_l_abs,k,n);
														clrbit(SD_l_abs,n);
														setbit(SD_l_abs,k);
													}
													if (SD_l_abs!=SD_r_abs)
													{
														sgn=1;
														if (q%2) sgn=-1;
														convert_to_ph_SD(SD_l_abs,SD_l_ph,sgn_l);
														sgn=sgn*sgn_l*sgn_r;
														val=sgn*V_it->second;
														add_H_m_element_sym(i, SD_l_ph, val);
													}
													else
													{
														val=V_it->second;
														add_H_m_element_diag(i, val);
													}
												}
												else
												{
													q=1;
													if (k!=m)
													{
														q+=add_bits_between(SD_l_abs,k,m);
														clrbit(SD_l_abs,m);
														setbit(SD_l_abs,k);
													}
													sgn=1;
													if (q%2) sgn=-1;
													if (SD_l_abs!=SD_r_abs)
													{
														convert_to_ph_SD(SD_l_abs,SD_l_ph,sgn_l);
														sgn=sgn*sgn_l*sgn_r;
														val=sgn*V_it->second;
														add_H_m_element_sym(i, SD_l_ph, val);
													}
													else
													{
														val=sgn*V_it->second;
														add_H_m_element_diag(i, val);
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		//down-down interaction if finite
		if (SO_H_T::V_d_SO_basis.size())
		{
			for (n=L; n<2*L-1; n++)
			{
				if (tstbit(SD_r_abs,n))
				{
					for (m=n+1; m<2*L; m++) //include only distinct pairs (m,n)
					{
						if (tstbit(SD_r_abs,m))
						{
							for (l=L; l<2*L-1; l++)
							{
								if (!tstbit(SD_r_abs,l) || l==m || l==n)
								{
									for (k=l+1; k<2*L; k++)
									{
										if (!tstbit(SD_r_abs,k) || k==n || k==m)
										{
											ind_V=k+n*s[1]+l*s[2]+m*s[3];
											V_it=SO_H_T::V_d_SO_basis.find(ind_V);
											if (V_it!=SO_H_T::V_d_SO_basis.end())
											{
												SD_l_abs=SD_r_abs;
												if (n!=l)
												{
													q=0;
													if (l!=m)
													{
														q+=add_bits_between(SD_l_abs,l,m);
														clrbit(SD_l_abs,m);
														setbit(SD_l_abs,l);
													}
													if (k!=n)
													{
														q+=add_bits_between(SD_l_abs,k,n);
														clrbit(SD_l_abs,n);
														setbit(SD_l_abs,k);
													}
													if (SD_l_abs!=SD_r_abs)
													{
														sgn=1;
														if (q%2) sgn=-1;
														convert_to_ph_SD(SD_l_abs,SD_l_ph,sgn_l);
														sgn=sgn*sgn_l*sgn_r;
														val=sgn*V_it->second;
														add_H_m_element_sym(i, SD_l_ph, val);
													}
													else
													{
														val=V_it->second;
														add_H_m_element_diag(i, val);
													}
												}
												else
												{
													q=1;
													if (k!=m)
													{
														q+=add_bits_between(SD_l_abs,k,m);
														clrbit(SD_l_abs,m);
														setbit(SD_l_abs,k);
													}
													sgn=1;
													if (q%2) sgn=-1;
													if (SD_l_abs!=SD_r_abs)
													{
														convert_to_ph_SD(SD_l_abs,SD_l_ph,sgn_l);
														sgn=sgn*sgn_l*sgn_r;
														val=sgn*V_it->second;
														add_H_m_element_sym(i, SD_l_ph, val);
													}
													else
													{
														val=sgn*V_it->second;
														add_H_m_element_diag(i, val);
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		//up-down interaction
		for (n=0; n<L; n++)
		{
			if (tstbit(SD_r_abs,n))
			{
				for (m=L; m<2*L; m++)
				{
					if (tstbit(SD_r_abs,m))
					{
						for (l=L; l<2*L; l++)
						{
							if (!tstbit(SD_r_abs,l) || l==m)
							{
								for (k=0; k<L; k++)
								{
									if (!tstbit(SD_r_abs,k) || k==n)
									{
										ind_V=k+n*s[1]+l*s[2]+m*s[3];
										V_it=SO_H_T::V_ud_SO_basis.find(ind_V);
										if (V_it!=SO_H_T::V_ud_SO_basis.end())
										{
											SD_l_abs=SD_r_abs;
											
											q=0;
											if (l!=m)
											{
												q+=add_bits_between(SD_l_abs,l,m);
												clrbit(SD_l_abs,m);
												setbit(SD_l_abs,l);
											}
											if (k!=n)
											{
												q+=add_bits_between(SD_l_abs,k,n);
												clrbit(SD_l_abs,n);
												setbit(SD_l_abs,k);
											}
											if (SD_l_abs!=SD_r_abs)
											{
												sgn=1;
												if (q%2) sgn=-1;
												convert_to_ph_SD(SD_l_abs,SD_l_ph,sgn_l);
												sgn=sgn*sgn_l*sgn_r;
												val=sgn*V_it->second;
												add_H_m_element_sym(i, SD_l_ph, val);
											}
											else
											{
												val=V_it->second;
												add_H_m_element_diag(i, val);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
	}

	
	i=NSD-1;
	SD_r_abs=CI_abs_basis[i];
	
	add_H_m_element_diag(i, -E_ref);
	
	//spin-up hopping
	for (l=0; l<L; l++)
	{
		if (tstbit(SD_r_abs,l))
		{
			//diagonal element
			val=SO_H_T::t_SO_basis_u[l][l];
			add_H_m_element_diag(i, val);
		}
	}
	//spin-down hopping
	for (l=L; l<2*L; l++)
	{
		if (tstbit(SD_r_abs,l))
		{
			//diagonal element
			val=SO_H_T::t_SO_basis_d[l-L][l-L];
			add_H_m_element_diag(i, val);
		}
	}
	
	//up-up interaction if finite
	if (SO_H_T::V_u_SO_basis.size())
	{
		for (n=0; n<L-1; n++)
		{
			if (tstbit(SD_r_abs,n))
			{
				for (m=n+1; m<L; m++)
				{
					if (tstbit(SD_r_abs,m))
					{
						ind_V=n+n*s[1]+m*s[2]+m*s[3];
						V_it=SO_H_T::V_u_SO_basis.find(ind_V);
						if (V_it!=SO_H_T::V_u_SO_basis.end())
						{
							val=V_it->second;
							add_H_m_element_diag(i, val);
						}
					}
				}
			}
		}
	}

	
	//down-down interaction if finite
	if (SO_H_T::V_d_SO_basis.size())
	{
		for (n=L; n<2*L-1; n++)
		{
			if (tstbit(SD_r_abs,n))
			{
				for (m=n+1; m<2*L; m++)
				{
					if (tstbit(SD_r_abs,m))
					{
						ind_V=n+n*s[1]+m*s[2]+m*s[3];
						V_it=SO_H_T::V_d_SO_basis.find(ind_V);
						if (V_it!=SO_H_T::V_d_SO_basis.end())
						{
							val=V_it->second;
							add_H_m_element_diag(i, val);
						}
					}
				}
			}
		}
	}
	
	//up-down interaction
	for (n=0; n<L; n++)
	{
		if (tstbit(SD_r_abs,n))
		{
			for (m=L; m<2*L; m++)
			{
				if (tstbit(SD_r_abs,m))
				{
					ind_V=n+n*s[1]+m*s[2]+m*s[3];
					V_it=SO_H_T::V_ud_SO_basis.find(ind_V);
					if (V_it!=SO_H_T::V_ud_SO_basis.end())
					{
						val=V_it->second;
						add_H_m_element_diag(i, val);
					}
				}
			}
		}
	}
	
	cout<<"Hamiltonian matrix defined.\n";
	
//	check_H_Hermicity();
	convert_H_m_to_H_v();
	
	cout<<"H[0][0]:\n";
	cout<<"index: "<<H_v[0][0].index<<endl;
	cout<<"value: "<<H_v[0][0].value<<endl;
	
}

template<class ind_T, class val_T, class SD_T, class cfs_T>
void Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::print_SD(SD_T SD_p)
{
	int j;
	for (j=2*L-1; j>=L; j--)
	{
		cout<<tstbit(SD_p,j);
	}
	cout<<' ';
	for (j=L-1; j>=0; j--)
	{
		cout<<tstbit(SD_p,j);
	}
}

template<class ind_T, class val_T, class SD_T, class cfs_T>
void Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::convert_to_ph_SD(const SD_T &SD_abs, SD_T &SD_ph, short int &sgn)
{
	SD_ph=0;
	
	int i;
	
	SD_T XO=SD_abs^ref_SD;
	
	SD_T pos_h=XO&ref_SD;
	SD_T pos_p=XO&SD_abs;
	
	for (i=0; i<2*L; i++)
	{
		if (tstbit(pos_h,i)) setbit(SD_ph,SO_ph_map[i]);
		else if (tstbit(pos_p,i)) setbit(SD_ph,SO_ph_map[i]);
	}
	
	ind_T j=get_SD_index(SD_ph);
	sgn=abs_to_ph_sign[j];
	
}

template<class ind_T, class val_T, class SD_T, class cfs_T>
void Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::convert_from_ph_SD(const SD_T &SD_ph, SD_T &SD_abs, short int &sgn)
{
	SD_abs=ref_SD;
	
	sgn=1;
	int m=0;
	int i, j;
	for (i=Nu-1; i>=0; i--)
	{
		if (tstbit(SD_ph,i))
		{
			m+=add_bits(SD_abs, 0,ph_indices[i]-1);
			clrbit(SD_abs,ph_indices[i]);
		}
	}
	for (i=L-1; i>=Nu; i--)
	{
		if (tstbit(SD_ph,i))
		{
			m+=add_bits(SD_abs,0,ph_indices[i]-1);
			setbit(SD_abs,ph_indices[i]);
		}
	}
	for (i=L+Nd-1; i>=L; i--)
	{
		if (tstbit(SD_ph,i))
		{
			m+=add_bits(SD_abs, L, ph_indices[i]-1);
			clrbit(SD_abs,ph_indices[i]);
		}
	}
	for (i=2*L-1; i>=L+Nd; i--)
	{
		if (tstbit(SD_ph,i))
		{
			m+=add_bits(SD_abs,L,ph_indices[i]-1);
			setbit(SD_abs,ph_indices[i]);
		}
	}
	
	if (m%2) sgn=-1;
	
}

template<class ind_T, class val_T, class SD_T, class cfs_T>
void Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::convert_H_m_to_H_v()
{
	H_v.clear();
	H_v.resize(NSD);
	
	ind_T i, j;
	auto it=H_m[0].begin();
	matrix_element<ind_T, val_T> elem;
	
	for (i=0; i<NSD; i++)
	{
		H_v[i].resize(H_m[i].size());
		j=0;
		for (it=H_m[i].begin(); it!=H_m[i].end(); it++)
		{
			elem.index=it->first;
			elem.value=it->second;
			H_v[i][j]=elem;
			j++;
		}
	}
//	H_m.clear();
	
}

template<class ind_T, class val_T, class SD_T, class cfs_T>
void Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::add_H_m_element_diag(ind_T i, val_T val)
{
	//check if elements exist at position (i,j). If so, add val to its value. Otherwise insert a new element.
	auto H_it=H_m[i].find(i);
	if (H_it!=H_m[i].end()) H_it->second+=val;
	else H_m[i].insert(pair<ind_T, val_T>(i,val));
}

template<class ind_T, class val_T, class SD_T, class cfs_T>
void Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::add_H_m_element(ind_T i, SD_T SD_tmp, val_T val)
{
	ind_T j=get_SD_index(SD_tmp);
	//check if elements exist at position (i,j). If so, add val to its value. Otherwise insert a new element.
	auto H_it=H_m[i].find(j);
	if (H_it!=H_m[i].end()) H_it->second+=val;
	else H_m[i].insert(pair<ind_T, val_T>(j,val));
}

template<class ind_T, class val_T, class SD_T, class cfs_T>
void Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::add_H_m_element_sym(ind_T i, SD_T SD_tmp, val_T val)
{
	ind_T j=get_SD_index(SD_tmp);
	
//	cout<<"add_H_m_element_sym:  "<<setw(20)<<SD_tmp<<setw(20)<<CI_ph_basis[j]<<CI_ph_basis[j]-SD_tmp<<endl;
	
	if (j>i)
	{
		//check if elements exist at positions (i,j) and (j,i) in H_m. If so, add val to their value. Otherwise insert new elements.
		auto H_it=H_m[i].find(j);
		if (H_it!=H_m[i].end())
		{
			val_T val_tmp=H_it->second;
			H_it->second+=conj(val);
			H_it=H_m[j].find(i);
			if (H_it!=H_m[j].end())
			{
				if (H_it->second!=conj(val_tmp))
				{
					cerr<<"add_H_m_element_sym() error: H is not Hermitian\n";
					cout<<setw(10)<<i;
					print_binary_string(CI_ph_basis[i],2*L);
					cout<<"\t"<<setw(10)<<j;
					print_binary_string(CI_ph_basis[j],2*L);
					cout<<endl;
					//	cout<<setw(20)<<val_tmp<<H_it->second<<endl;
				}
				H_it->second+=val;
			}
			else
			{
				cerr<<"add_H_m_element_sym() error: H is not Hermitian\n";
				//	cout<<setw(10)<<i;
				//	print_binary_string(CI_ph_basis[i],2*L);
				//	cout<<"\t"<<setw(10)<<j;
				//	print_binary_string(CI_ph_basis[j],2*L);
				//	cout<<endl;
			}
		}
		else //insert the element at (i,j) and (j,i)
		{
			H_m[i].insert(pair<ind_T, val_T>(j,conj(val)));
			H_m[j].insert(pair<ind_T, val_T>(i,val));
			/*
			 if (add_bits(CI_ph_basis[i], 0, 2*L-1)==add_bits(CI_ph_basis[j], 0, 2*L-1))
			 {
			 cout<<setw(10)<<i;
			 print_binary_string(CI_ph_basis[i],2*L);
			 cout<<"\t"<<setw(10)<<j;
			 print_binary_string(CI_ph_basis[j],2*L);
			 cout<<endl;
			 }
			 */
		}
	}
}

template<class ind_T, class val_T, class SD_T, class cfs_T>
void Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::insert_H_m_element_sym(ind_T i, SD_T SD_tmp, val_T val)
{
	ind_T j=get_SD_index(SD_tmp);
	H_m[i].insert(pair<ind_T, val_T>(j,conj(val)));
	H_m[j].insert(pair<ind_T, val_T>(i,val));
}

template<class ind_T, class val_T, class SD_T, class cfs_T>
val_T Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::elem_V(int i, int j, int k, int l)
{
	ind_T ind=ph_indices[i]+ph_indices[j]*s[1]+ph_indices[k]*s[2]+ph_indices[l]*s[3];
	
	val_T V_tmp=0;
//	if (i<L && k<L)
	if (i<L && j<L && k<L && l<L)
	{
		auto V_it=SO_H_T::V_u_SO_basis.find(ind);
		if (V_it!=SO_H_T::V_u_SO_basis.end()) V_tmp=V_it->second;
	}
//	else if (i>=L && k>=L)
	else if (i>=L && j>=L && k>=L && l>=L)
	{
		auto V_it=SO_H_T::V_d_SO_basis.find(ind);
		if (V_it!=SO_H_T::V_d_SO_basis.end()) V_tmp=V_it->second;
	}
	else
	{
		auto V_it=SO_H_T::V_ud_SO_basis.find(ind);
		if (V_it!=SO_H_T::V_ud_SO_basis.end()) V_tmp=V_it->second;
	}
	
	return V_tmp;
}

template<class ind_T, class val_T, class SD_T, class cfs_T>
val_T Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::elem_t_Fock(int i, int j)
{
	if (i<L && j<L)
	{
		return SO_H_T::t_Fock_u[ph_indices[i]][ph_indices[j]];
	}
	else if (i>=L && j>=L)
	{
		return SO_H_T::t_Fock_d[ph_indices[i]-L][ph_indices[j]-L];
	}
	else
	{
		cerr<<"elem_t_Fock() warning: Fock matrix is diagonal with respect to spin\n";
		return 0;
	}
		
}

template<class ind_T, class val_T, class SD_T, class cfs_T>
void Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::create_H_matrix_ph()
{
	cout<<"Defining Hamiltonian matrix...\n";
	
	compute_ref_diag_E();
	
	H_m.clear();
	H_m.resize(NSD);
	
	SD_T SD_tmp;
	ind_T i;
	int q, k, l, m, n, nu, nd, ind;
	val_T t_tmp, V_tmp, sgn;
	auto V_it=SO_H_T::V_u_SO_basis.begin();
	
	//for the reference SD
	i=0;
	
	//creation of a single particle-hole pair by the one-body term
	
	//creation of a up-down pair
	for (l=0; l<Nu; l++)
	{
		for (k=Nu; k<L; k++)
		{
			SD_tmp=CI_ph_basis[i];
			
			t_tmp=elem_t_Fock(k,l);
			//	t_tmp=SO_H_T::t_Fock_u[ph_indices[k]][ph_indices[l]];
			
			if (t_tmp)
			{
				setbit(SD_tmp,k);
				setbit(SD_tmp,l);
				add_H_m_element_sym(i, SD_tmp, t_tmp);
			}
		}
	}
	
	//creation of a down-up pair
	for (l=L; l<L+Nd; l++)
	{
		for (k=L+Nd; k<2*L; k++)
		{
			SD_tmp=CI_ph_basis[i];
			
			t_tmp=elem_t_Fock(k,l);
			//	t_tmp=SO_H_T::t_Fock_d[ph_indices[k]-L][ph_indices[l]-L];
			
			if (t_tmp)
			{
				setbit(SD_tmp,k);
				setbit(SD_tmp,l);
				add_H_m_element_sym(i, SD_tmp, t_tmp);
			}
		}
	}
	
	//Creation of two particle-hole pairs
	//Each distinct combination appears only once, so that there is no prefactor.
	
	//creation of up-down and one down-up pairs
	for (n=0; n<Nu; n++)
	{
		for (m=L; m<L+Nd; m++)
		{
			for (l=L+Nd; l<2*L; l++)
			{
				for (k=Nu; k<L; k++)
				{
					SD_tmp=CI_ph_basis[i];
					
					V_tmp=elem_V(k,n,l,m);
					
					if (V_tmp)
					{
						setbit(SD_tmp,l);
						setbit(SD_tmp,m);
						setbit(SD_tmp,k);
						setbit(SD_tmp,n);
						add_H_m_element_sym(i, SD_tmp, V_tmp);
					}
				}
			}
		}
	}
	
	//creation of two up-down pairs
	for (n=0; n<Nu-1; n++)
	{
		for (m=n+1; m<Nu; m++)
		{
			for (l=Nu; l<L-1; l++)
			{
				for (k=l+1; k<L; k++)
				{
					SD_tmp=CI_ph_basis[i];
					
					V_tmp=elem_V(k,n,l,m);
					
					if (V_tmp)
					{
						setbit(SD_tmp,l);
						setbit(SD_tmp,m);
						q=add_bits_between(SD_tmp,k,n);
						setbit(SD_tmp,k);
						setbit(SD_tmp,n);
						sgn=1;
						if (q%2) sgn=-1;
						V_tmp=sgn*V_tmp;
						add_H_m_element_sym(i, SD_tmp, V_tmp);
					}
				}
			}
		}
	}
	
	//creation of two down-up pairs
	for (n=L; n<L+Nd-1; n++)
	{
		for (m=n+1; m<L+Nd; m++)
		{
			for (l=L+Nd; l<2*L-1; l++)
			{
				for (k=l+1; k<2*L; k++)
				{
					SD_tmp=CI_ph_basis[i];
					
					V_tmp=elem_V(k,n,l,m);
					
					if (V_tmp)
					{
						setbit(SD_tmp,l);
						setbit(SD_tmp,m);
						q=add_bits_between(SD_tmp,k,n);
						setbit(SD_tmp,k);
						setbit(SD_tmp,n);
						sgn=1;
						if (q%2) sgn=-1;
						V_tmp=sgn*V_tmp;
						add_H_m_element_sym(i, SD_tmp, V_tmp);
					}
				}
			}
		}
	}
	
	//for all other SD's
	for (i=1; i<NSD; i++)
	{
		SD_tmp=CI_ph_basis[i];
		
		nu=add_bits(SD_tmp, 0, L-1)/2;  //number of excitations in the spin-up part of the SD
		nd=add_bits(SD_tmp, L, 2*L-1)/2; //number of excitations in the spin-down part
		
		//creation of an up-down pair or hopping of spin-down hole
		for (l=0; l<Nu; l++)
		{
			if (!tstbit(CI_ph_basis[i],l)) //creation of a pair
			{
				for (k=Nu; k<L; k++)
				{
					if (!tstbit(CI_ph_basis[i],k))
					{
						SD_tmp=CI_ph_basis[i];
						
						t_tmp=elem_t_Fock(k,l);
						//	t_tmp=SO_H_T::t_Fock_u[ph_indices[k]][ph_indices[l]];
						
						if (t_tmp)
						{
							q=add_bits_between(SD_tmp,k,l);
							setbit(SD_tmp,k);
							setbit(SD_tmp,l);
							sgn=1;
							if (q%2) sgn=-1;
							t_tmp=sgn*t_tmp;
							add_H_m_element_sym(i, SD_tmp, t_tmp);
						}
					}
				}
			}
			else //hopping of a hole if the orbital is occupied
			{
				for (k=0; k<Nu; k++)
				{
					if (!tstbit(CI_ph_basis[i],k)) //hopping
					{
						SD_tmp=CI_ph_basis[i];
						
						t_tmp=-conj(elem_t_Fock(k,l));
						//	t_tmp=-conj(SO_H_T::t_Fock_u[ph_indices[k]][ph_indices[l]]);
						
						if (t_tmp)
						{
							q=add_bits_between(SD_tmp,k,l);
							setbit(SD_tmp,k);
							clrbit(SD_tmp,l);
							sgn=1;
							if (q%2) sgn=-1;
							t_tmp=sgn*t_tmp;
							add_H_m_element_sym(i, SD_tmp, t_tmp);
						}
					}
				}
				//onsite energy
				t_tmp=-elem_t_Fock(l,l);
				//	t_tmp=-SO_H_T::t_Fock_u[ph_indices[l]][ph_indices[l]];
				
				if (t_tmp)
				{
					add_H_m_element_diag(i, t_tmp);
				}
			}
		}
		
		//creation of an down-up pair or hopping of spin-up hole
		for (l=L; l<L+Nd; l++)
		{
			if (!tstbit(CI_ph_basis[i],l)) //creation of a pair
			{
				for (k=L+Nd; k<2*L; k++)
				{
					if (!tstbit(CI_ph_basis[i],k))
					{
						SD_tmp=CI_ph_basis[i];
						
						t_tmp=elem_t_Fock(k,l);
						//	t_tmp=SO_H_T::t_Fock_d[ph_indices[k]-L][ph_indices[l]-L];
						
						if (t_tmp)
						{
							q=add_bits_between(SD_tmp,k,l);
							setbit(SD_tmp,k);
							setbit(SD_tmp,l);
							sgn=1;
							if (q%2) sgn=-1;
							t_tmp=sgn*t_tmp;
							add_H_m_element_sym(i, SD_tmp, t_tmp);
						}
					}
				}
			}
			else //hopping of a hole if the orbital is occupied
			{
				for (k=L; k<L+Nd; k++)
				{
					if (!tstbit(CI_ph_basis[i],k)) //hopping
					{
						SD_tmp=CI_ph_basis[i];
						
						t_tmp=-conj(elem_t_Fock(k,l));
						//	t_tmp=-conj(SO_H_T::t_Fock_d[ph_indices[k]-L][ph_indices[l]-L]);
						
						if (t_tmp)
						{
							q=add_bits_between(SD_tmp,k,l);
							setbit(SD_tmp,k);
							clrbit(SD_tmp,l);
							sgn=1;
							if (q%2) sgn=-1;
							t_tmp=sgn*t_tmp;
							add_H_m_element_sym(i, SD_tmp, t_tmp);
						}
					}
				}
				//onsite energy
				t_tmp=-elem_t_Fock(l,l);
				//	t_tmp=-SO_H_T::t_Fock_d[ph_indices[l]-L][ph_indices[l]-L];
				
				if (t_tmp)
				{
					add_H_m_element_diag(i, t_tmp);
				}
			}
		}
		
		
		if (nu>0)
		{
			//hopping of a spin-up particle
			for (l=Nu; l<L; l++)
			{
				if (tstbit(CI_ph_basis[i],l))
				{
					for (k=Nu; k<L; k++)
					{
						if (!tstbit(CI_ph_basis[i],k))
						{
							SD_tmp=CI_ph_basis[i];
							
							t_tmp=elem_t_Fock(k,l);
							//	t_tmp=SO_H_T::t_Fock_u[ph_indices[k]][ph_indices[l]];
							
							if (t_tmp)
							{
								q=add_bits_between(SD_tmp,k,l);
								setbit(SD_tmp,k);
								clrbit(SD_tmp,l);
								sgn=1;
								if (q%2) sgn=-1;
								t_tmp=sgn*t_tmp;
								add_H_m_element_sym(i, SD_tmp, t_tmp);
							}
						}
					}
					//onsite energy
					t_tmp=elem_t_Fock(l,l);
					//	t_tmp=SO_H_T::t_Fock_u[ph_indices[l]][ph_indices[l]];
					
					if (t_tmp)
					{
						add_H_m_element_diag(i, t_tmp);
					}
				}
			}
		}
		
		if (nd>0)
		{
			//hopping of a spin-down particle
			for (l=L+Nd; l<2*L; l++)
			{
				if (tstbit(CI_ph_basis[i],l))
				{
					for (k=L+Nd; k<2*L; k++)
					{
						if (!tstbit(CI_ph_basis[i],k))
						{
							SD_tmp=CI_ph_basis[i];
							
							t_tmp=elem_t_Fock(k,l);
							//	t_tmp=SO_H_T::t_Fock_d[ph_indices[k]-L][ph_indices[l]-L];
							
							if (t_tmp)
							{
								q=add_bits_between(SD_tmp,k,l);
								setbit(SD_tmp,k);
								clrbit(SD_tmp,l);
								sgn=1;
								if (q%2) sgn=-1;
								t_tmp=sgn*t_tmp;
								add_H_m_element_sym(i, SD_tmp, t_tmp);
							}
						}
					}
					//onsite energy
					t_tmp=elem_t_Fock(l,l);
					//	t_tmp=SO_H_T::t_Fock_d[ph_indices[l]-L][ph_indices[l]-L];
					
					if (t_tmp)
					{
						add_H_m_element_diag(i, t_tmp);
					}
				}
			}
		}
		
		
		//creation of two particle-hole pairs by the two-body term V_knlm(p^+_k)(p^+_l)(h^+_m)(h^+_n)=V_knlm(p^+_k)(h^+_n)(p^+_l)(h^+_m)
		//Each distinct combination appears only once, so that there is no prefactor.
		
		if (nu<Nph_u_max && nd<Nph_d_max)
		{
			//creation of one up-down and one down-up pairs
			for (n=0; n<Nu; n++)
			{
				if (!tstbit(CI_ph_basis[i],n)) //unoccupied spin-down hole orbital
				{
					for (m=L; m<L+Nd; m++)
					{
						if (!tstbit(CI_ph_basis[i],m)) //unoccupied spin-up hole orbital
						{
							for (l=L+Nd; l<2*L; l++)
							{
								if (!tstbit(CI_ph_basis[i],l)) //unoccupied spin-down particle orbital
								{
									for (k=Nu; k<L; k++)
									{
										if (!tstbit(CI_ph_basis[i],k)) //unoccupied spin-up particle orbital
										{
											SD_tmp=CI_ph_basis[i];
											
											V_tmp=elem_V(k,n,l,m);
											
											if (V_tmp)
											{
												q=add_bits_between(SD_tmp,l,m);
												setbit(SD_tmp,l);
												setbit(SD_tmp,m);
												q+=add_bits_between(SD_tmp,k,n);
												setbit(SD_tmp,k);
												setbit(SD_tmp,n);
												sgn=1;
												if (q%2) sgn=-1;
												V_tmp=sgn*V_tmp;
												add_H_m_element_sym(i, SD_tmp, V_tmp);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		if (nu<Nph_u_max-1)
		{
			//creation of two up-down pairs
			for (n=0; n<Nu-1; n++) //spin-down hole orbital
			{
				if (!tstbit(CI_ph_basis[i],n)) //unoccupied
				{
					for (m=n+1; m<Nu; m++) //spin-down hole orbital
					{
						if (!tstbit(CI_ph_basis[i],m)) //unoccupied
						{
							for (l=Nu; l<L-1; l++) //spin-up particle orbital
							{
								if (!tstbit(CI_ph_basis[i],l)) //unoccupied
								{
									for (k=l+1; k<L; k++) //spin-up particle orbital
									{
										if (!tstbit(CI_ph_basis[i],k)) //unoccupied
										{
											SD_tmp=CI_ph_basis[i];
											
											V_tmp=elem_V(k,n,l,m);
											
											if (V_tmp)
											{
												q=add_bits_between(SD_tmp,l,m);
												setbit(SD_tmp,l);
												setbit(SD_tmp,m);
												q+=add_bits_between(SD_tmp,k,n);
												setbit(SD_tmp,k);
												setbit(SD_tmp,n);
												sgn=1;
												if (q%2) sgn=-1;
												V_tmp=sgn*V_tmp;
												add_H_m_element_sym(i, SD_tmp, V_tmp);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		if (nd<Nph_d_max-1)
		{
			//creation of two down-up pairs
			for (n=L; n<L+Nd-1; n++)
			{
				if (!tstbit(CI_ph_basis[i],n))
				{
					for (m=n+1; m<L+Nd; m++)
					{
						if (!tstbit(CI_ph_basis[i],m))
						{
							for (l=L+Nd; l<2*L-1; l++)
							{
								if (!tstbit(CI_ph_basis[i],l))
								{
									for (k=l+1; k<2*L; k++)
									{
										if (!tstbit(CI_ph_basis[i],k))
										{
											SD_tmp=CI_ph_basis[i];
											
											V_tmp=elem_V(k,n,l,m);
											
											if (V_tmp)
											{
												q=add_bits_between(SD_tmp,l,m);
												setbit(SD_tmp,l);
												setbit(SD_tmp,m);
												q+=add_bits_between(SD_tmp,k,n);
												setbit(SD_tmp,k);
												setbit(SD_tmp,n);
												sgn=1;
												if (q%2) sgn=-1;
												V_tmp=sgn*V_tmp;
												add_H_m_element_sym(i, SD_tmp, V_tmp);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		//creation of one particle-hole pair and hopping of a particle
		
		//spin-up particle and up-down pair
		if (nu>0 && nu<Nph_u_max)
		{
			for (n=Nu; n<L; n++)
			{
				if (tstbit(CI_ph_basis[i],n))
				{
					for (m=0; m<Nu; m++)
					{
						if (!tstbit(CI_ph_basis[i],m))
						{
							for (l=Nu; l<L-1; l++)
							{
								for (k=l+1; k<L; k++)
								{
									if ( (!tstbit(CI_ph_basis[i],k) && !tstbit(CI_ph_basis[i],l)) || (!tstbit(CI_ph_basis[i],k) && l==n) || (!tstbit(CI_ph_basis[i],l) && k==n) )
									{
										SD_tmp=CI_ph_basis[i];
										
										V_tmp=elem_V(k,n,l,m);
										
										if (V_tmp)
										{
											if (!tstbit(CI_ph_basis[i],l))
											{
												q=add_bits_between(SD_tmp,l,m);
												setbit(SD_tmp,l);
												setbit(SD_tmp,m);
												if (k!=n)
												{
													q+=add_bits_between(SD_tmp,k,n);
													setbit(SD_tmp,k);
													clrbit(SD_tmp,n);
												}
											}
											else
											{
												q=1;
												q+=add_bits_between(SD_tmp,k,m);
												setbit(SD_tmp,k);
												setbit(SD_tmp,m);
											}
											sgn=1;
											if (q%2) sgn=-1;
											V_tmp=sgn*V_tmp;
											add_H_m_element_sym(i, SD_tmp, V_tmp);
										}
									}
								}
								
							}
						}
					}
				}
			}
		}
		
		//spin-down particle and down-up pair
		if (nd>0 && nd<Nph_d_max)
		{
			for (n=L+Nd; n<2*L; n++)
			{
				if (tstbit(CI_ph_basis[i],n))
				{
					for (m=L; m<L+Nd; m++)
					{
						if (!tstbit(CI_ph_basis[i],m))
						{
							for (l=L+Nd; l<2*L-1; l++)
							{
								for (k=l+1; k<2*L; k++)
								{
									if ( (!tstbit(CI_ph_basis[i],k) && !tstbit(CI_ph_basis[i],l)) || (!tstbit(CI_ph_basis[i],k) && l==n) || (!tstbit(CI_ph_basis[i],l) && k==n) )
									{
										SD_tmp=CI_ph_basis[i];
										
										V_tmp=elem_V(k,n,l,m);
										
										if (V_tmp)
										{
											if (!tstbit(CI_ph_basis[i],l))
											{
												q=add_bits_between(SD_tmp,l,m);
												setbit(SD_tmp,l);
												setbit(SD_tmp,m);
												if (k!=n)
												{
													q+=add_bits_between(SD_tmp,k,n);
													setbit(SD_tmp,k);
													clrbit(SD_tmp,n);
												}
											}
											else
											{
												q=1;
												q+=add_bits_between(SD_tmp,k,m);
												setbit(SD_tmp,k);
												setbit(SD_tmp,m);
											}
											sgn=1;
											if (q%2) sgn=-1;
											V_tmp=sgn*V_tmp;
											add_H_m_element_sym(i, SD_tmp, V_tmp);
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		//spin-up particle and down-up pair
		if (nu>0 && nd<Nph_d_max)
		{
			for (n=Nu; n<L; n++)
			{
				if (tstbit(CI_ph_basis[i],n))
				{
					for (m=L; m<L+Nd; m++)
					{
						if (!tstbit(CI_ph_basis[i],m))
						{
							for (l=L+Nd; l<2*L; l++)
							{
								if (!tstbit(CI_ph_basis[i],l))
								{
									for (k=Nu; k<L; k++)
									{
										if ( !tstbit(CI_ph_basis[i],k) || k==n )
										{
											SD_tmp=CI_ph_basis[i];
											
											V_tmp=elem_V(k,n,l,m);
											
											if (V_tmp)
											{
												q=add_bits_between(SD_tmp,l,m);
												setbit(SD_tmp,l);
												setbit(SD_tmp,m);
												if (k!=n)
												{
													q+=add_bits_between(SD_tmp,k,n);
													setbit(SD_tmp,k);
													clrbit(SD_tmp,n);
												}
												sgn=1;
												if (q%2) sgn=-1;
												V_tmp=sgn*V_tmp;
												add_H_m_element_sym(i, SD_tmp, V_tmp);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		//spin-down particle and up-down pair
		if (nd>0 && nu<Nph_u_max)
		{
			for (n=L+Nd; n<2*L; n++)
			{
				if (tstbit(CI_ph_basis[i],n))
				{
					for (m=0; m<Nu; m++)
					{
						if (!tstbit(CI_ph_basis[i],m))
						{
							for (l=Nu; l<L; l++)
							{
								if (!tstbit(CI_ph_basis[i],l))
								{
									for (k=L+Nd; k<2*L; k++)
									{
										if ( !tstbit(CI_ph_basis[i],k) || k==n )
										{
											SD_tmp=CI_ph_basis[i];
											
											V_tmp=elem_V(k,n,l,m);
											
											if (V_tmp)
											{
												q=add_bits_between(SD_tmp,l,m);
												setbit(SD_tmp,l);
												setbit(SD_tmp,m);
												if (k!=n)
												{
													q+=add_bits_between(SD_tmp,k,n);
													setbit(SD_tmp,k);
													clrbit(SD_tmp,n);
												}
												sgn=1;
												if (q%2) sgn=-1;
												V_tmp=sgn*V_tmp;
												add_H_m_element_sym(i, SD_tmp, V_tmp);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		//creation of one particle-hole pair and hopping of a hole
		
		//spin-down hole and up-down pair
		if (nu>0 && nu<Nph_u_max)
		{
			for (n=0; n<Nu; n++)
			{
				if (tstbit(CI_ph_basis[i],n))
				{
					for (l=Nu; l<L; l++)
					{
						if (!tstbit(CI_ph_basis[i],l))
						{
							for (m=0; m<Nu-1; m++)
							{
								for (k=m+1; k<Nu; k++)
								{
									if ( (!tstbit(CI_ph_basis[i],k) && !tstbit(CI_ph_basis[i],m)) || (!tstbit(CI_ph_basis[i],k) && m==n) || (!tstbit(CI_ph_basis[i],m) && k==n) )
									{
										SD_tmp=CI_ph_basis[i];
										
										V_tmp=-conj(elem_V(k,n,m,l));
										
										if (V_tmp)
										{
											if (!tstbit(CI_ph_basis[i],m))
											{
												q=0;
												if (k!=n)
												{
													q+=add_bits_between(SD_tmp,k,n);
													setbit(SD_tmp,k);
													clrbit(SD_tmp,n);
												}
												q+=add_bits_between(SD_tmp,l,m);
												setbit(SD_tmp,l);
												setbit(SD_tmp,m);
											}
											else
											{
												q=1;
												q+=add_bits_between(SD_tmp,l,k);
												setbit(SD_tmp,k);
												setbit(SD_tmp,l);
											}
											sgn=1;
											if (q%2) sgn=-1;
											V_tmp=sgn*V_tmp;
											add_H_m_element_sym(i, SD_tmp, V_tmp);
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		//spin-up hole and down-up pair
		if (nd>0 && nd<Nph_d_max)
		{
			for (n=L; n<L+Nd; n++)
			{
				if (tstbit(CI_ph_basis[i],n))
				{
					for (l=L+Nd; l<2*L; l++)
					{
						if (!tstbit(CI_ph_basis[i],l))
						{
							for (m=L; m<L+Nd-1; m++)
							{
								for (k=m+1; k<L+Nd; k++)
								{
									if ( (!tstbit(CI_ph_basis[i],k) && !tstbit(CI_ph_basis[i],m)) || (!tstbit(CI_ph_basis[i],k) && m==n) || (!tstbit(CI_ph_basis[i],m) && k==n) )
									{
										SD_tmp=CI_ph_basis[i];
										
										V_tmp=-conj(elem_V(k,n,m,l));
										
										if (V_tmp)
										{
											if (!tstbit(CI_ph_basis[i],m))
											{
												q=0;
												if (k!=n)
												{
													q+=add_bits_between(SD_tmp,k,n);
													setbit(SD_tmp,k);
													clrbit(SD_tmp,n);
												}
												q+=add_bits_between(SD_tmp,l,m);
												setbit(SD_tmp,l);
												setbit(SD_tmp,m);
											}
											else
											{
												q=1;
												q+=add_bits_between(SD_tmp,l,k);
												setbit(SD_tmp,k);
												setbit(SD_tmp,l);
											}
											sgn=1;
											if (q%2) sgn=-1;
											V_tmp=sgn*V_tmp;
											add_H_m_element_sym(i, SD_tmp, V_tmp);
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		//spin-down hole and down-up pair
		if (nu>0 && nd<Nph_d_max)
		{
			for (n=0; n<Nu; n++)
			{
				if (tstbit(CI_ph_basis[i],n))
				{
					for (l=L+Nd; l<2*L; l++)
					{
						if (!tstbit(CI_ph_basis[i],l))
						{
							for (m=L; m<L+Nd; m++)
							{
								if (!tstbit(CI_ph_basis[i],m))
								{
									for (k=0; k<Nu; k++)
									{
										if ( !tstbit(CI_ph_basis[i],k) || k==n )
										{
											SD_tmp=CI_ph_basis[i];
											
											V_tmp=-conj(elem_V(k,n,m,l));
											
											if (V_tmp)
											{
												q=0;
												if (k!=n)
												{
													q+=add_bits_between(SD_tmp,k,n);
													setbit(SD_tmp,k);
													clrbit(SD_tmp,n);
												}
												q+=add_bits_between(SD_tmp,l,m);
												setbit(SD_tmp,l);
												setbit(SD_tmp,m);
												sgn=1;
												if (q%2) sgn=-1;
												V_tmp=sgn*V_tmp;
												add_H_m_element_sym(i, SD_tmp, V_tmp);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		//spin-up hole and up-down pair
		if (nd>0 && nu<Nph_u_max)
		{
			for (n=L; n<L+Nd; n++)
			{
				if (tstbit(CI_ph_basis[i],n))
				{
					for (l=Nu; l<L; l++)
					{
						if (!tstbit(CI_ph_basis[i],l))
						{
							for (m=0; m<Nu; m++)
							{
								if (!tstbit(CI_ph_basis[i],m))
								{
									for (k=L; k<L+Nd; k++)
									{
										if ( !tstbit(CI_ph_basis[i],k) || k==n )
										{
											SD_tmp=CI_ph_basis[i];
											
											V_tmp=-conj(elem_V(k,n,m,l));
											
											if (V_tmp)
											{
												q=0;
												if (k!=n)
												{
													q+=add_bits_between(SD_tmp,k,n);
													setbit(SD_tmp,k);
													clrbit(SD_tmp,n);
												}
												q+=add_bits_between(SD_tmp,l,m);
												setbit(SD_tmp,l);
												setbit(SD_tmp,m);
												sgn=1;
												if (q%2) sgn=-1;
												V_tmp=sgn*V_tmp;
												add_H_m_element_sym(i, SD_tmp, V_tmp);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		
		//hopping of two particles
		
		//two spin-up particles
		if (nu>1)
		{
			for (n=Nu; n<L-1; n++)
			{
				if (tstbit(CI_ph_basis[i],n))
				{
					for (m=n+1; m<L; m++)
					{
						if (tstbit(CI_ph_basis[i],m))
						{
							for (l=Nu; l<L-1; l++)
							{
								if ( !tstbit(CI_ph_basis[i],l) || l==m || l==n)
								{
									for (k=l+1; k<L; k++)
									{
										if ( !tstbit(CI_ph_basis[i],k) || k==m || k==n)
										{
											SD_tmp=CI_ph_basis[i];
											
											V_tmp=elem_V(k,n,l,m);
											
											if (V_tmp)
											{
												if (l!=n)
												{
													q=0;
													if (l!=m)
													{
														q+=add_bits_between(SD_tmp,l,m);
														setbit(SD_tmp,l);
														clrbit(SD_tmp,m);
													}
													if (k!=n)
													{
														q+=add_bits_between(SD_tmp,k,n);
														setbit(SD_tmp,k);
														clrbit(SD_tmp,n);
													}
												}
												else
												{
													q=1;
													if (k!=m)
													{
														q+=add_bits_between(SD_tmp,k,m);
														setbit(SD_tmp,k);
														clrbit(SD_tmp,m);
													}
												}
												sgn=1;
												if (q%2) sgn=-1;
												V_tmp=sgn*V_tmp;
												if (SD_tmp!=CI_ph_basis[i]) add_H_m_element_sym(i, SD_tmp, V_tmp);
												else add_H_m_element_diag(i, V_tmp);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		//two spin-down particles
		if (nd>1)
		{
			for (n=L+Nd; n<2*L-1; n++)
			{
				if (tstbit(CI_ph_basis[i],n))
				{
					for (m=n+1; m<2*L; m++)
					{
						if (tstbit(CI_ph_basis[i],m))
						{
							for (l=L+Nd; l<2*L-1; l++)
							{
								if ( !tstbit(CI_ph_basis[i],l) || l==m || l==n)
								{
									for (k=l+1; k<2*L; k++)
									{
										if ( !tstbit(CI_ph_basis[i],k) || k==m || k==n)
										{
											SD_tmp=CI_ph_basis[i];
											
											V_tmp=elem_V(k,n,l,m);
											
											if (V_tmp)
											{
												if (l!=n)
												{
													q=0;
													if (l!=m)
													{
														q+=add_bits_between(SD_tmp,l,m);
														setbit(SD_tmp,l);
														clrbit(SD_tmp,m);
													}
													if (k!=n)
													{
														q+=add_bits_between(SD_tmp,k,n);
														setbit(SD_tmp,k);
														clrbit(SD_tmp,n);
													}
												}
												else
												{
													q=1;
													if (k!=m)
													{
														q+=add_bits_between(SD_tmp,k,m);
														setbit(SD_tmp,k);
														clrbit(SD_tmp,m);
													}
												}
												sgn=1;
												if (q%2) sgn=-1;
												V_tmp=sgn*V_tmp;
												if (SD_tmp!=CI_ph_basis[i]) add_H_m_element_sym(i, SD_tmp, V_tmp);
												else add_H_m_element_diag(i, V_tmp);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		//one spin-up and one spin-down particles
		if (nu>0 && nd>0)
		{
			for (n=L+Nd; n<2*L; n++)
			{
				if (tstbit(CI_ph_basis[i],n))
				{
					for (m=Nu; m<L; m++)
					{
						if (tstbit(CI_ph_basis[i],m))
						{
							for (l=Nu; l<L; l++)
							{
								if ( !tstbit(CI_ph_basis[i],l) || l==m)
								{
									for (k=L+Nd; k<2*L; k++)
									{
										if ( !tstbit(CI_ph_basis[i],k) || k==n)
										{
											SD_tmp=CI_ph_basis[i];
											
											V_tmp=elem_V(k,n,l,m);
											
											if (V_tmp)
											{
												q=0;
												if (l!=m)
												{
													q+=add_bits_between(SD_tmp,l,m);
													setbit(SD_tmp,l);
													clrbit(SD_tmp,m);
												}
												if (k!=n)
												{
													q+=add_bits_between(SD_tmp,k,n);
													setbit(SD_tmp,k);
													clrbit(SD_tmp,n);
												}
												sgn=1;
												if (q%2) sgn=-1;
												V_tmp=sgn*V_tmp;
												if (SD_tmp!=CI_ph_basis[i]) add_H_m_element_sym(i, SD_tmp, V_tmp);
												else add_H_m_element_diag(i, V_tmp);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		
		//hopping of two holes
		
		//two spin-down holes
		if (nu>1)
		{
			for (n=0; n<Nu-1; n++)
			{
				if (tstbit(CI_ph_basis[i],n))
				{
					for (m=n+1; m<Nu; m++)
					{
						if (tstbit(CI_ph_basis[i],m))
						{
							for (l=0; l<Nu-1; l++)
							{
								if ( !tstbit(CI_ph_basis[i],l) || l==m || l==n)
								{
									for (k=l+1; k<Nu; k++)
									{
										if ( !tstbit(CI_ph_basis[i],k) || k==m || k==n)
										{
											SD_tmp=CI_ph_basis[i];
											
											V_tmp=conj(elem_V(k,n,l,m));
											
											if (V_tmp)
											{
												if (l!=n)
												{
													q=0;
													if (l!=m)
													{
														q+=add_bits_between(SD_tmp,l,m);
														setbit(SD_tmp,l);
														clrbit(SD_tmp,m);
													}
													if (k!=n)
													{
														q+=add_bits_between(SD_tmp,k,n);
														setbit(SD_tmp,k);
														clrbit(SD_tmp,n);
													}
												}
												else
												{
													q=1;
													if (k!=m)
													{
														q+=add_bits_between(SD_tmp,k,m);
														setbit(SD_tmp,k);
														clrbit(SD_tmp,m);
													}
												}
												sgn=1;
												if (q%2) sgn=-1;
												V_tmp=sgn*V_tmp;
												if (SD_tmp!=CI_ph_basis[i]) add_H_m_element_sym(i, SD_tmp, V_tmp);
												else add_H_m_element_diag(i, V_tmp);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		//two spin-up holes
		if (nd>1)
		{
			for (n=L; n<L+Nd-1; n++)
			{
				if (tstbit(CI_ph_basis[i],n))
				{
					for (m=n+1; m<L+Nd; m++)
					{
						if (tstbit(CI_ph_basis[i],m))
						{
							for (l=L; l<L+Nd-1; l++)
							{
								if ( !tstbit(CI_ph_basis[i],l) || l==m || l==n)
								{
									for (k=l+1; k<L+Nd; k++)
									{
										if ( !tstbit(CI_ph_basis[i],k) || k==m || k==n)
										{
											SD_tmp=CI_ph_basis[i];
											
											V_tmp=conj(elem_V(k,n,l,m));
											
											if (V_tmp)
											{
												if (l!=n)
												{
													q=0;
													if (l!=m)
													{
														q+=add_bits_between(SD_tmp,l,m);
														setbit(SD_tmp,l);
														clrbit(SD_tmp,m);
													}
													if (k!=n)
													{
														q+=add_bits_between(SD_tmp,k,n);
														setbit(SD_tmp,k);
														clrbit(SD_tmp,n);
													}
												}
												else
												{
													q=1;
													if (k!=m)
													{
														q+=add_bits_between(SD_tmp,k,m);
														setbit(SD_tmp,k);
														clrbit(SD_tmp,m);
													}
												}
												sgn=1;
												if (q%2) sgn=-1;
												V_tmp=sgn*V_tmp;
												if (SD_tmp!=CI_ph_basis[i]) add_H_m_element_sym(i, SD_tmp, V_tmp);
												else add_H_m_element_diag(i, V_tmp);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		//one spin-down and one spin-up holes
		if (nu>0 && nd>0)
		{
			for (n=0; n<Nu; n++)
			{
				if (tstbit(CI_ph_basis[i],n))
				{
					for (m=L; m<L+Nd; m++)
					{
						if (tstbit(CI_ph_basis[i],m))
						{
							for (l=L; l<L+Nd; l++)
							{
								if ( !tstbit(CI_ph_basis[i],l) || l==m )
								{
									for (k=0; k<Nu; k++)
									{
										if ( !tstbit(CI_ph_basis[i],k) || k==n)
										{
											SD_tmp=CI_ph_basis[i];
											
											V_tmp=conj(elem_V(k,n,l,m));
											
											if (V_tmp)
											{
												q=0;
												if (l!=m)
												{
													q+=add_bits_between(SD_tmp,l,m);
													setbit(SD_tmp,l);
													clrbit(SD_tmp,m);
												}
												if (k!=n)
												{
													q+=add_bits_between(SD_tmp,k,n);
													setbit(SD_tmp,k);
													clrbit(SD_tmp,n);
												}
												sgn=1;
												if (q%2) sgn=-1;
												V_tmp=sgn*V_tmp;
												if (SD_tmp!=CI_ph_basis[i]) add_H_m_element_sym(i, SD_tmp, V_tmp);
												else add_H_m_element_diag(i, V_tmp);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		
		//hopping of one particle and one hole
		
		//one spin-up particle and one spin-down hole
		if (nu>0)
		{
			for (n=Nu; n<L; n++)
			{
				if (tstbit(CI_ph_basis[i],n))
				{
					for (m=0; m<Nu; m++)
					{
						if (tstbit(CI_ph_basis[i],m))
						{
							for (l=0; l<Nu; l++)
							{
								if ( !tstbit(CI_ph_basis[i],l) || l==m )
								{
									for (k=Nu; k<L; k++)
									{
										if ( !tstbit(CI_ph_basis[i],k) || k==n)
										{
											SD_tmp=CI_ph_basis[i];
											
											V_tmp=-elem_V(k,n,m,l);
											
											if (V_tmp)
											{
												q=0;
												if (l!=m)
												{
													q+=add_bits_between(SD_tmp,l,m);
													setbit(SD_tmp,l);
													clrbit(SD_tmp,m);
												}
												if (k!=n)
												{
													q+=add_bits_between(SD_tmp,k,n);
													setbit(SD_tmp,k);
													clrbit(SD_tmp,n);
												}
												sgn=1;
												if (q%2) sgn=-1;
												V_tmp=sgn*V_tmp;
												if (SD_tmp!=CI_ph_basis[i]) add_H_m_element_sym(i, SD_tmp, V_tmp);
												else add_H_m_element_diag(i, V_tmp);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		//one spin-down particle and one spin-up hole
		if (nd>0)
		{
			for (n=L+Nd; n<2*L; n++)
			{
				if (tstbit(CI_ph_basis[i],n))
				{
					for (m=L; m<L+Nd; m++)
					{
						if (tstbit(CI_ph_basis[i],m))
						{
							for (l=L; l<L+Nd; l++)
							{
								if ( !tstbit(CI_ph_basis[i],l) || l==m )
								{
									for (k=L+Nd; k<2*L; k++)
									{
										if ( !tstbit(CI_ph_basis[i],k) || k==n)
										{
											SD_tmp=CI_ph_basis[i];
											
											V_tmp=-elem_V(k,n,m,l);
											
											if (V_tmp)
											{
												q=0;
												if (l!=m)
												{
													q+=add_bits_between(SD_tmp,l,m);
													setbit(SD_tmp,l);
													clrbit(SD_tmp,m);
												}
												if (k!=n)
												{
													q+=add_bits_between(SD_tmp,k,n);
													setbit(SD_tmp,k);
													clrbit(SD_tmp,n);
												}
												sgn=1;
												if (q%2) sgn=-1;
												V_tmp=sgn*V_tmp;
												if (SD_tmp!=CI_ph_basis[i]) add_H_m_element_sym(i, SD_tmp, V_tmp);
												else add_H_m_element_diag(i, V_tmp);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		//one spin-up particle and one spin-up hole
		if (nu>0 && nd>0)
		{
			for (n=Nu; n<L; n++)
			{
				if (tstbit(CI_ph_basis[i],n))
				{
					for (m=L; m<L+Nd; m++)
					{
						if (tstbit(CI_ph_basis[i],m))
						{
							for (l=L; l<L+Nd; l++)
							{
								if ( !tstbit(CI_ph_basis[i],l) || l==m )
								{
									for (k=Nu; k<L; k++)
									{
										if ( !tstbit(CI_ph_basis[i],k) || k==n)
										{
											SD_tmp=CI_ph_basis[i];
											
											V_tmp=-elem_V(k,n,m,l);
											
											if (V_tmp)
											{
												q=0;
												if (l!=m)
												{
													q+=add_bits_between(SD_tmp,l,m);
													setbit(SD_tmp,l);
													clrbit(SD_tmp,m);
												}
												if (k!=n)
												{
													q+=add_bits_between(SD_tmp,k,n);
													setbit(SD_tmp,k);
													clrbit(SD_tmp,n);
												}
												sgn=1;
												if (q%2) sgn=-1;
												V_tmp=sgn*V_tmp;
												if (SD_tmp!=CI_ph_basis[i]) add_H_m_element_sym(i, SD_tmp, V_tmp);
												else add_H_m_element_diag(i, V_tmp);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		//one spin-down particle and one spin-down hole
		if (nu>0 && nd>0)
		{
			for (n=L+Nd; n<2*L; n++)
			{
				if (tstbit(CI_ph_basis[i],n))
				{
					for (m=0; m<Nu; m++)
					{
						if (tstbit(CI_ph_basis[i],m))
						{
							for (l=0; l<Nu; l++)
							{
								if ( !tstbit(CI_ph_basis[i],l) || l==m )
								{
									for (k=L+Nd; k<2*L; k++)
									{
										if ( !tstbit(CI_ph_basis[i],k) || k==n)
										{
											SD_tmp=CI_ph_basis[i];
											
											V_tmp=-elem_V(k,n,m,l);
											
											if (V_tmp)
											{
												q=0;
												if (l!=m)
												{
													q+=add_bits_between(SD_tmp,l,m);
													setbit(SD_tmp,l);
													clrbit(SD_tmp,m);
												}
												if (k!=n)
												{
													q+=add_bits_between(SD_tmp,k,n);
													setbit(SD_tmp,k);
													clrbit(SD_tmp,n);
												}
												sgn=1;
												if (q%2) sgn=-1;
												V_tmp=sgn*V_tmp;
												if (SD_tmp!=CI_ph_basis[i]) add_H_m_element_sym(i, SD_tmp, V_tmp);
												else add_H_m_element_diag(i, V_tmp);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		//destruction of a up-down pair and creation of a down-up pair
		if (nu>0 && nd<Nph_d_max)
		{
			for (n=Nu; n<L; n++)
			{
				if (tstbit(CI_ph_basis[i],n))
				{
					for (m=0; m<Nu; m++)
					{
						if (tstbit(CI_ph_basis[i],m))
						{
							for (l=L; l<L+Nd; l++)
							{
								if ( !tstbit(CI_ph_basis[i],l) )
								{
									for (k=L+Nd; k<2*L; k++)
									{
										if ( !tstbit(CI_ph_basis[i],k) )
										{
											SD_tmp=CI_ph_basis[i];
											
											V_tmp=-elem_V(k,n,m,l);
											
											if (V_tmp)
											{
												q=add_bits_between(SD_tmp,n,m);
												clrbit(SD_tmp,m);
												clrbit(SD_tmp,n);
												q+=add_bits_between(SD_tmp,k,l);
												setbit(SD_tmp,k);
												setbit(SD_tmp,l);
												sgn=1;
												if (q%2) sgn=-1;
												V_tmp=sgn*V_tmp;
												add_H_m_element_sym(i, SD_tmp, V_tmp);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		//destruction of a down-up pair and creation of a up-down pair
		if (nd>0 && nu<Nph_u_max)
		{
			for (n=L+Nd; n<2*L; n++)
			{
				if (tstbit(CI_ph_basis[i],n))
				{
					for (m=L; m<L+Nd; m++)
					{
						if (tstbit(CI_ph_basis[i],m))
						{
							for (l=0; l<Nu; l++)
							{
								if ( !tstbit(CI_ph_basis[i],l) )
								{
									for (k=Nu; k<L; k++)
									{
										if ( !tstbit(CI_ph_basis[i],k) )
										{
											SD_tmp=CI_ph_basis[i];
											
											V_tmp=-elem_V(k,n,m,l);
											
											if (V_tmp)
											{
												q=add_bits_between(SD_tmp,n,m);
												clrbit(SD_tmp,m);
												clrbit(SD_tmp,n);
												q+=add_bits_between(SD_tmp,k,l);
												setbit(SD_tmp,k);
												setbit(SD_tmp,l);
												sgn=1;
												if (q%2) sgn=-1;
												V_tmp=sgn*V_tmp;
												add_H_m_element_sym(i, SD_tmp, V_tmp);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
	}//end of loop over basis SD's
	
	cout<<"Hamiltonian matrix defined\n";
	
	convert_H_m_to_H_v();
}

/*
template<class ind_T, class val_T, class SD_T, class cfs_T>
void Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::create_H_matrix_ph()
{
	cout<<"Defining Hamiltonian matrix...\n";
	
	compute_ref_diag_E();
	
	H_m.clear();
	H_m.resize(NSD);
	
	SD_T SD_tmp;
	ind_T i;
	int q, k, l, m, n, nu, nd, ind;
	val_T t_tmp, V_tmp, sgn;
	auto V_it=SO_H_T::V_u_SO_basis.begin();
	
	//for the reference SD
	i=0;
	
	//creation of a single particle-hole pair by the one-body term
	
	//creation of a up-down pair
	for (l=0; l<Nu; l++)
	{
		for (k=Nu; k<L; k++)
		{
			SD_tmp=CI_ph_basis[i];
			
			t_tmp=elem_t_Fock(k,l);
		//	t_tmp=SO_H_T::t_Fock_u[ph_indices[k]][ph_indices[l]];
			
			if (t_tmp)
			{
				setbit(SD_tmp,k);
				setbit(SD_tmp,l);
				add_H_m_element_sym(i, SD_tmp, t_tmp);
			}
		}
	}
	
	//creation of a down-up pair
	for (l=L; l<L+Nd; l++)
	{
		for (k=L+Nd; k<2*L; k++)
		{
			SD_tmp=CI_ph_basis[i];
			
			t_tmp=elem_t_Fock(k,l);
		//	t_tmp=SO_H_T::t_Fock_d[ph_indices[k]-L][ph_indices[l]-L];
			
			if (t_tmp)
			{
				setbit(SD_tmp,k);
				setbit(SD_tmp,l);
				add_H_m_element_sym(i, SD_tmp, t_tmp);
			}
		}
	}
	
	//Creation of two particle-hole pairs
	//Each distinct combination appears only once, so that there is no prefactor.
	
	//creation of up-down and one down-up pairs
	for (n=0; n<Nu; n++)
	{
		for (m=L; m<L+Nd; m++)
		{
			for (l=L+Nd; l<2*L; l++)
			{
				for (k=Nu; k<L; k++)
				{
					SD_tmp=CI_ph_basis[i];
					
					V_tmp=elem_V(k,n,l,m);
 
					if (V_tmp)
					{
						setbit(SD_tmp,l);
						setbit(SD_tmp,m);
						setbit(SD_tmp,k);
						setbit(SD_tmp,n);
						add_H_m_element_sym(i, SD_tmp, V_tmp);
					}
				}
			}
		}
	}
	
	//creation of two up-down pairs
	for (n=0; n<Nu-1; n++)
	{
		for (m=n+1; m<Nu; m++)
		{
			for (l=Nu; l<L-1; l++)
			{
				for (k=l+1; k<L; k++)
				{
					SD_tmp=CI_ph_basis[i];
					
					V_tmp=elem_V(k,n,l,m);
				
					if (V_tmp)
					{
						setbit(SD_tmp,l);
						setbit(SD_tmp,m);
						q=add_bits_between(SD_tmp,k,n);
						setbit(SD_tmp,k);
						setbit(SD_tmp,n);
						sgn=1;
						if (q%2) sgn=-1;
						V_tmp=sgn*V_tmp;
						add_H_m_element_sym(i, SD_tmp, V_tmp);
					}
				}
			}
		}
	}
	
	//creation of two down-up pairs
	for (n=L; n<L+Nd-1; n++)
	{
		for (m=n+1; m<L+Nd; m++)
		{
			for (l=L+Nd; l<2*L-1; l++)
			{
				for (k=l+1; k<2*L; k++)
				{
					SD_tmp=CI_ph_basis[i];
					
					V_tmp=elem_V(k,n,l,m);
				
					if (V_tmp)
					{
						setbit(SD_tmp,l);
						setbit(SD_tmp,m);
						q=add_bits_between(SD_tmp,k,n);
						setbit(SD_tmp,k);
						setbit(SD_tmp,n);
						sgn=1;
						if (q%2) sgn=-1;
						V_tmp=sgn*V_tmp;
						add_H_m_element_sym(i, SD_tmp, V_tmp);
					}
				}
			}
		}
	}
	
	//for all other SD's
	for (i=1; i<NSD; i++)
	{
		SD_tmp=CI_ph_basis[i];
		
		nu=add_bits(SD_tmp, 0, L-1)/2;  //number of excitations in the spin-up part of the SD
		nd=add_bits(SD_tmp, L, 2*L-1)/2; //number of excitations in the spin-down part
		
		//creation of an up-down pair or hopping of spin-down hole
		for (l=0; l<Nu; l++)
		{
			if (!tstbit(CI_ph_basis[i],l)) //creation of a pair
			{
				for (k=Nu; k<L; k++)
				{
					if (!tstbit(CI_ph_basis[i],k))
					{
						SD_tmp=CI_ph_basis[i];
						
						t_tmp=elem_t_Fock(k,l);
					//	t_tmp=SO_H_T::t_Fock_u[ph_indices[k]][ph_indices[l]];
						
						if (t_tmp)
						{
							q=add_bits_between(SD_tmp,k,l);
							setbit(SD_tmp,k);
							setbit(SD_tmp,l);
							sgn=1;
							if (q%2) sgn=-1;
							t_tmp=sgn*t_tmp;
							add_H_m_element_sym(i, SD_tmp, t_tmp);
						}
					}
				}
			}
			else //hopping of a hole if the orbital is occupied
			{
				for (k=0; k<Nu; k++)
				{
					if (!tstbit(CI_ph_basis[i],k)) //hopping
					{
						SD_tmp=CI_ph_basis[i];
						
						t_tmp=-conj(elem_t_Fock(k,l));
					//	t_tmp=-conj(SO_H_T::t_Fock_u[ph_indices[k]][ph_indices[l]]);
						
						if (t_tmp)
						{
							q=add_bits_between(SD_tmp,k,l);
							setbit(SD_tmp,k);
							clrbit(SD_tmp,l);
							sgn=1;
							if (q%2) sgn=-1;
							t_tmp=sgn*t_tmp;
							add_H_m_element_sym(i, SD_tmp, t_tmp);
						}
					}
				}
				//onsite energy
				t_tmp=-elem_t_Fock(l,l);
			//	t_tmp=-SO_H_T::t_Fock_u[ph_indices[l]][ph_indices[l]];
				
				if (t_tmp)
				{
					add_H_m_element_diag(i, t_tmp);
				}
			}
		}
		
		//creation of an down-up pair or hopping of spin-up hole
		for (l=L; l<L+Nd; l++)
		{
			if (!tstbit(CI_ph_basis[i],l)) //creation of a pair
			{
				for (k=L+Nd; k<2*L; k++)
				{
					if (!tstbit(CI_ph_basis[i],k))
					{
						SD_tmp=CI_ph_basis[i];
						
						t_tmp=elem_t_Fock(k,l);
					//	t_tmp=SO_H_T::t_Fock_d[ph_indices[k]-L][ph_indices[l]-L];
						
						if (t_tmp)
						{
							q=add_bits_between(SD_tmp,k,l);
							setbit(SD_tmp,k);
							setbit(SD_tmp,l);
							sgn=1;
							if (q%2) sgn=-1;
							t_tmp=sgn*t_tmp;
							add_H_m_element_sym(i, SD_tmp, t_tmp);
						}
					}
				}
			}
			else //hopping of a hole if the orbital is occupied
			{
				for (k=L; k<L+Nd; k++)
				{
					if (!tstbit(CI_ph_basis[i],k)) //hopping
					{
						SD_tmp=CI_ph_basis[i];
						
						t_tmp=-conj(elem_t_Fock(k,l));
					//	t_tmp=-conj(SO_H_T::t_Fock_d[ph_indices[k]-L][ph_indices[l]-L]);
						
						if (t_tmp)
						{
							q=add_bits_between(SD_tmp,k,l);
							setbit(SD_tmp,k);
							clrbit(SD_tmp,l);
							sgn=1;
							if (q%2) sgn=-1;
							t_tmp=sgn*t_tmp;
							add_H_m_element_sym(i, SD_tmp, t_tmp);
						}
					}
				}
				//onsite energy
				t_tmp=-elem_t_Fock(l,l);
			//	t_tmp=-SO_H_T::t_Fock_d[ph_indices[l]-L][ph_indices[l]-L];
				
				if (t_tmp)
				{
					add_H_m_element_diag(i, t_tmp);
				}
			}
		}
		
		
		if (nu>0)
		{
			//hopping of a spin-up particle
			for (l=Nu; l<L; l++)
			{
				if (tstbit(CI_ph_basis[i],l))
				{
					for (k=Nu; k<L; k++)
					{
						if (!tstbit(CI_ph_basis[i],k))
						{
							SD_tmp=CI_ph_basis[i];
							
							t_tmp=elem_t_Fock(k,l);
						//	t_tmp=SO_H_T::t_Fock_u[ph_indices[k]][ph_indices[l]];
							
							if (t_tmp)
							{
								q=add_bits_between(SD_tmp,k,l);
								setbit(SD_tmp,k);
								clrbit(SD_tmp,l);
								sgn=1;
								if (q%2) sgn=-1;
								t_tmp=sgn*t_tmp;
								add_H_m_element_sym(i, SD_tmp, t_tmp);
							}
						}
					}
					//onsite energy
					t_tmp=elem_t_Fock(l,l);
				//	t_tmp=SO_H_T::t_Fock_u[ph_indices[l]][ph_indices[l]];
					
					if (t_tmp)
					{
						add_H_m_element_diag(i, t_tmp);
					}
				}
			}
		}
		
		if (nd>0)
		{
			//hopping of a spin-down particle
			for (l=L+Nd; l<2*L; l++)
			{
				if (tstbit(CI_ph_basis[i],l))
				{
					for (k=L+Nd; k<2*L; k++)
					{
						if (!tstbit(CI_ph_basis[i],k))
						{
							SD_tmp=CI_ph_basis[i];
							
							t_tmp=elem_t_Fock(k,l);
						//	t_tmp=SO_H_T::t_Fock_d[ph_indices[k]-L][ph_indices[l]-L];
							
							if (t_tmp)
							{
								q=add_bits_between(SD_tmp,k,l);
								setbit(SD_tmp,k);
								clrbit(SD_tmp,l);
								sgn=1;
								if (q%2) sgn=-1;
								t_tmp=sgn*t_tmp;
								add_H_m_element_sym(i, SD_tmp, t_tmp);
							}
						}
					}
					//onsite energy
					t_tmp=elem_t_Fock(l,l);
				//	t_tmp=SO_H_T::t_Fock_d[ph_indices[l]-L][ph_indices[l]-L];
					
					if (t_tmp)
					{
						add_H_m_element_diag(i, t_tmp);
					}
				}
			}
		}
		
		
		//creation of two particle-hole pairs by the two-body term V_knlm(p^+_k)(p^+_l)(h^+_m)(h^+_n)=V_knlm(p^+_k)(h^+_n)(p^+_l)(h^+_m)
		//Each distinct combination appears only once, so that there is no prefactor.
		
		if (nu<Nph_u_max && nd<Nph_d_max)
		{
			//creation of one up-down and one down-up pairs
			for (n=0; n<Nu; n++)
			{
				if (!tstbit(CI_ph_basis[i],n)) //unoccupied spin-down hole orbital
				{
					for (m=L; m<L+Nd; m++)
					{
						if (!tstbit(CI_ph_basis[i],m)) //unoccupied spin-up hole orbital
						{
							for (l=L+Nd; l<2*L; l++)
							{
								if (!tstbit(CI_ph_basis[i],l)) //unoccupied spin-down particle orbital
								{
									for (k=Nu; k<L; k++)
									{
										if (!tstbit(CI_ph_basis[i],k)) //unoccupied spin-up particle orbital
										{
											SD_tmp=CI_ph_basis[i];
											
											V_tmp=elem_V(k,n,l,m);
										
											if (V_tmp)
											{
												q=add_bits_between(SD_tmp,l,m);
												setbit(SD_tmp,l);
												setbit(SD_tmp,m);
												q+=add_bits_between(SD_tmp,k,n);
												setbit(SD_tmp,k);
												setbit(SD_tmp,n);
												sgn=1;
												if (q%2) sgn=-1;
												V_tmp=sgn*V_tmp;
												add_H_m_element_sym(i, SD_tmp, V_tmp);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		if (nu<Nph_u_max-1)
		{
			//creation of two up-down pairs
			for (n=0; n<Nu-1; n++) //spin-down hole orbital
			{
				if (!tstbit(CI_ph_basis[i],n)) //unoccupied
				{
					for (m=n+1; m<Nu; m++) //spin-down hole orbital
					{
						if (!tstbit(CI_ph_basis[i],m)) //unoccupied
						{
							for (l=Nu; l<L-1; l++) //spin-up particle orbital
							{
								if (!tstbit(CI_ph_basis[i],l)) //unoccupied
								{
									for (k=l+1; k<L; k++) //spin-up particle orbital
									{
										if (!tstbit(CI_ph_basis[i],k)) //unoccupied
										{
											SD_tmp=CI_ph_basis[i];
											
											V_tmp=elem_V(k,n,l,m);
										
											if (V_tmp)
											{
												q=add_bits_between(SD_tmp,l,m);
												setbit(SD_tmp,l);
												setbit(SD_tmp,m);
												q+=add_bits_between(SD_tmp,k,n);
												setbit(SD_tmp,k);
												setbit(SD_tmp,n);
												sgn=1;
												if (q%2) sgn=-1;
												V_tmp=sgn*V_tmp;
												add_H_m_element_sym(i, SD_tmp, V_tmp);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		if (nd<Nph_d_max-1)
		{
			//creation of two down-up pairs
			for (n=L; n<L+Nd-1; n++)
			{
				if (!tstbit(CI_ph_basis[i],n))
				{
					for (m=n+1; m<L+Nd; m++)
					{
						if (!tstbit(CI_ph_basis[i],m))
						{
							for (l=L+Nd; l<2*L-1; l++)
							{
								if (!tstbit(CI_ph_basis[i],l))
								{
									for (k=l+1; k<2*L; k++)
									{
										if (!tstbit(CI_ph_basis[i],k))
										{
											SD_tmp=CI_ph_basis[i];
											
											V_tmp=elem_V(k,n,l,m);
										
											if (V_tmp)
											{
												q=add_bits_between(SD_tmp,l,m);
												setbit(SD_tmp,l);
												setbit(SD_tmp,m);
												q+=add_bits_between(SD_tmp,k,n);
												setbit(SD_tmp,k);
												setbit(SD_tmp,n);
												sgn=1;
												if (q%2) sgn=-1;
												V_tmp=sgn*V_tmp;
												add_H_m_element_sym(i, SD_tmp, V_tmp);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		//creation of one particle-hole pair and hopping of a particle
		
		//spin-up particle and up-down pair
		if (nu>0 && nu<Nph_u_max)
		{
			for (n=Nu; n<L; n++)
			{
				if (tstbit(CI_ph_basis[i],n))
				{
					for (m=0; m<Nu; m++)
					{
						if (!tstbit(CI_ph_basis[i],m))
						{
							for (l=Nu; l<L-1; l++)
							{
								for (k=l+1; k<L; k++)
								{
									if ( (!tstbit(CI_ph_basis[i],k) && !tstbit(CI_ph_basis[i],l)) || (!tstbit(CI_ph_basis[i],k) && l==n) || (!tstbit(CI_ph_basis[i],l) && k==n) )
									{
										SD_tmp=CI_ph_basis[i];
										
										V_tmp=elem_V(k,n,l,m);
									
										if (V_tmp)
										{
											if (!tstbit(CI_ph_basis[i],l))
											{
												q=add_bits_between(SD_tmp,l,m);
												setbit(SD_tmp,l);
												setbit(SD_tmp,m);
												if (k!=n)
												{
													q+=add_bits_between(SD_tmp,k,n);
													setbit(SD_tmp,k);
													clrbit(SD_tmp,n);
												}
											}
											else
											{
												q=1;
												q+=add_bits_between(SD_tmp,k,m);
												setbit(SD_tmp,k);
												setbit(SD_tmp,m);
											}
											sgn=1;
											if (q%2) sgn=-1;
											V_tmp=sgn*V_tmp;
											add_H_m_element_sym(i, SD_tmp, V_tmp);
										}
									}
								}
								
							}
						}
					}
				}
			}
		}
		
		//spin-down particle and down-up pair
		if (nd>0 && nd<Nph_d_max)
		{
			for (n=L+Nd; n<2*L; n++)
			{
				if (tstbit(CI_ph_basis[i],n))
				{
					for (m=L; m<L+Nd; m++)
					{
						if (!tstbit(CI_ph_basis[i],m))
						{
							for (l=L+Nd; l<2*L-1; l++)
							{
								for (k=l+1; k<2*L; k++)
								{
									if ( (!tstbit(CI_ph_basis[i],k) && !tstbit(CI_ph_basis[i],l)) || (!tstbit(CI_ph_basis[i],k) && l==n) || (!tstbit(CI_ph_basis[i],l) && k==n) )
									{
										SD_tmp=CI_ph_basis[i];
										
										V_tmp=elem_V(k,n,l,m);
									
										if (V_tmp)
										{
											if (!tstbit(CI_ph_basis[i],l))
											{
												q=add_bits_between(SD_tmp,l,m);
												setbit(SD_tmp,l);
												setbit(SD_tmp,m);
												if (k!=n)
												{
													q+=add_bits_between(SD_tmp,k,n);
													setbit(SD_tmp,k);
													clrbit(SD_tmp,n);
												}
											}
											else
											{
												q=1;
												q+=add_bits_between(SD_tmp,k,m);
												setbit(SD_tmp,k);
												setbit(SD_tmp,m);
											}
											sgn=1;
											if (q%2) sgn=-1;
											V_tmp=sgn*V_tmp;
											add_H_m_element_sym(i, SD_tmp, V_tmp);
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		//spin-up particle and down-up pair
		if (nu>0 && nd<Nph_d_max)
		{
			for (n=Nu; n<L; n++)
			{
				if (tstbit(CI_ph_basis[i],n))
				{
					for (m=L; m<L+Nd; m++)
					{
						if (!tstbit(CI_ph_basis[i],m))
						{
							for (l=L+Nd; l<2*L; l++)
							{
								if (!tstbit(CI_ph_basis[i],l))
								{
									for (k=Nu; k<L; k++)
									{
										if ( !tstbit(CI_ph_basis[i],k) || k==n )
										{
											SD_tmp=CI_ph_basis[i];
											
											V_tmp=elem_V(k,n,l,m);
										
											if (V_tmp)
											{
												q=add_bits_between(SD_tmp,l,m);
												setbit(SD_tmp,l);
												setbit(SD_tmp,m);
												if (k!=n)
												{
													q+=add_bits_between(SD_tmp,k,n);
													setbit(SD_tmp,k);
													clrbit(SD_tmp,n);
												}
												sgn=1;
												if (q%2) sgn=-1;
												V_tmp=sgn*V_tmp;
												add_H_m_element_sym(i, SD_tmp, V_tmp);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		//spin-down particle and up-down pair
		if (nd>0 && nu<Nph_u_max)
		{
			for (n=L+Nd; n<2*L; n++)
			{
				if (tstbit(CI_ph_basis[i],n))
				{
					for (m=0; m<Nu; m++)
					{
						if (!tstbit(CI_ph_basis[i],m))
						{
							for (l=Nu; l<L; l++)
							{
								if (!tstbit(CI_ph_basis[i],l))
								{
									for (k=L+Nd; k<2*L; k++)
									{
										if ( !tstbit(CI_ph_basis[i],k) || k==n )
										{
											SD_tmp=CI_ph_basis[i];
											
											V_tmp=elem_V(k,n,l,m);
										
											if (V_tmp)
											{
												q=add_bits_between(SD_tmp,l,m);
												setbit(SD_tmp,l);
												setbit(SD_tmp,m);
												if (k!=n)
												{
													q+=add_bits_between(SD_tmp,k,n);
													setbit(SD_tmp,k);
													clrbit(SD_tmp,n);
												}
												sgn=1;
												if (q%2) sgn=-1;
												V_tmp=sgn*V_tmp;
												add_H_m_element_sym(i, SD_tmp, V_tmp);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		//creation of one particle-hole pair and hopping of a hole
		
		//spin-down hole and up-down pair
		if (nu>0 && nu<Nph_u_max)
		{
			for (n=0; n<Nu; n++)
			{
				if (tstbit(CI_ph_basis[i],n))
				{
					for (l=Nu; l<L; l++)
					{
						if (!tstbit(CI_ph_basis[i],l))
						{
							for (m=0; m<Nu-1; m++)
							{
								for (k=m+1; k<Nu; k++)
								{
									if ( (!tstbit(CI_ph_basis[i],k) && !tstbit(CI_ph_basis[i],m)) || (!tstbit(CI_ph_basis[i],k) && m==n) || (!tstbit(CI_ph_basis[i],m) && k==n) )
									{
										SD_tmp=CI_ph_basis[i];
										
										V_tmp=-conj(elem_V(k,n,m,l));
									
										if (V_tmp)
										{
											if (!tstbit(CI_ph_basis[i],m))
											{
												q=0;
												if (k!=n)
												{
													q+=add_bits_between(SD_tmp,k,n);
													setbit(SD_tmp,k);
													clrbit(SD_tmp,n);
												}
												q+=add_bits_between(SD_tmp,l,m);
												setbit(SD_tmp,l);
												setbit(SD_tmp,m);
											}
											else
											{
												q=1;
												q+=add_bits_between(SD_tmp,l,k);
												setbit(SD_tmp,k);
												setbit(SD_tmp,l);
											}
											sgn=1;
											if (q%2) sgn=-1;
											V_tmp=sgn*V_tmp;
											add_H_m_element_sym(i, SD_tmp, V_tmp);
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		//spin-up hole and down-up pair
		if (nd>0 && nd<Nph_d_max)
		{
			for (n=L; n<L+Nd; n++)
			{
				if (tstbit(CI_ph_basis[i],n))
				{
					for (l=L+Nd; l<2*L; l++)
					{
						if (!tstbit(CI_ph_basis[i],l))
						{
							for (m=L; m<L+Nd-1; m++)
							{
								for (k=m+1; k<L+Nd; k++)
								{
									if ( (!tstbit(CI_ph_basis[i],k) && !tstbit(CI_ph_basis[i],m)) || (!tstbit(CI_ph_basis[i],k) && m==n) || (!tstbit(CI_ph_basis[i],m) && k==n) )
									{
										SD_tmp=CI_ph_basis[i];
										
										V_tmp=-conj(elem_V(k,n,m,l));
									
										if (V_tmp)
										{
											if (!tstbit(CI_ph_basis[i],m))
											{
												q=0;
												if (k!=n)
												{
													q+=add_bits_between(SD_tmp,k,n);
													setbit(SD_tmp,k);
													clrbit(SD_tmp,n);
												}
												q+=add_bits_between(SD_tmp,l,m);
												setbit(SD_tmp,l);
												setbit(SD_tmp,m);
											}
											else
											{
												q=1;
												q+=add_bits_between(SD_tmp,l,k);
												setbit(SD_tmp,k);
												setbit(SD_tmp,l);
											}
											sgn=1;
											if (q%2) sgn=-1;
											V_tmp=sgn*V_tmp;
											add_H_m_element_sym(i, SD_tmp, V_tmp);
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		//spin-down hole and down-up pair
		if (nu>0 && nd<Nph_d_max)
		{
			for (n=0; n<Nu; n++)
			{
				if (tstbit(CI_ph_basis[i],n))
				{
					for (l=L+Nd; l<2*L; l++)
					{
						if (!tstbit(CI_ph_basis[i],l))
						{
							for (m=L; m<L+Nd; m++)
							{
								if (!tstbit(CI_ph_basis[i],m))
								{
									for (k=0; k<Nu; k++)
									{
										if ( !tstbit(CI_ph_basis[i],k) || k==n )
										{
											SD_tmp=CI_ph_basis[i];
											
											V_tmp=-conj(elem_V(k,n,m,l));
										
											if (V_tmp)
											{
												q=0;
												if (k!=n)
												{
													q+=add_bits_between(SD_tmp,k,n);
													setbit(SD_tmp,k);
													clrbit(SD_tmp,n);
												}
												q+=add_bits_between(SD_tmp,l,m);
												setbit(SD_tmp,l);
												setbit(SD_tmp,m);
												sgn=1;
												if (q%2) sgn=-1;
												V_tmp=sgn*V_tmp;
												add_H_m_element_sym(i, SD_tmp, V_tmp);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		//spin-up hole and up-down pair
		if (nd>0 && nu<Nph_u_max)
		{
			for (n=L; n<L+Nd; n++)
			{
				if (tstbit(CI_ph_basis[i],n))
				{
					for (l=Nu; l<L; l++)
					{
						if (!tstbit(CI_ph_basis[i],l))
						{
							for (m=0; m<Nu; m++)
							{
								if (!tstbit(CI_ph_basis[i],m))
								{
									for (k=L; k<L+Nd; k++)
									{
										if ( !tstbit(CI_ph_basis[i],k) || k==n )
										{
											SD_tmp=CI_ph_basis[i];
											
											V_tmp=-conj(elem_V(k,n,m,l));
										
											if (V_tmp)
											{
												q=0;
												if (k!=n)
												{
													q+=add_bits_between(SD_tmp,k,n);
													setbit(SD_tmp,k);
													clrbit(SD_tmp,n);
												}
												q+=add_bits_between(SD_tmp,l,m);
												setbit(SD_tmp,l);
												setbit(SD_tmp,m);
												sgn=1;
												if (q%2) sgn=-1;
												V_tmp=sgn*V_tmp;
												add_H_m_element_sym(i, SD_tmp, V_tmp);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		
		//hopping of two particles
		
		//two spin-up particles
		if (nu>1)
		{
			for (n=Nu; n<L-1; n++)
			{
				if (tstbit(CI_ph_basis[i],n))
				{
					for (m=n+1; m<L; m++)
					{
						if (tstbit(CI_ph_basis[i],m))
						{
							for (l=Nu; l<L-1; l++)
							{
								if ( !tstbit(CI_ph_basis[i],l) || l==m || l==n)
								{
									for (k=l+1; k<L; k++)
									{
										if ( !tstbit(CI_ph_basis[i],k) || k==m || k==n)
										{
											SD_tmp=CI_ph_basis[i];
											
											V_tmp=elem_V(k,n,l,m);
										
											if (V_tmp)
											{
												if (l!=n)
												{
													q=0;
													if (l!=m)
													{
														q+=add_bits_between(SD_tmp,l,m);
														setbit(SD_tmp,l);
														clrbit(SD_tmp,m);
													}
													if (k!=n)
													{
														q+=add_bits_between(SD_tmp,k,n);
														setbit(SD_tmp,k);
														clrbit(SD_tmp,n);
													}
												}
												else
												{
													q=1;
													if (k!=m)
													{
														q+=add_bits_between(SD_tmp,k,m);
														setbit(SD_tmp,k);
														clrbit(SD_tmp,m);
													}
												}
												sgn=1;
												if (q%2) sgn=-1;
												V_tmp=sgn*V_tmp;
												if (SD_tmp!=CI_ph_basis[i]) add_H_m_element_sym(i, SD_tmp, V_tmp);
												else add_H_m_element_diag(i, V_tmp);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		//two spin-down particles
		if (nd>1)
		{
			for (n=L+Nd; n<2*L-1; n++)
			{
				if (tstbit(CI_ph_basis[i],n))
				{
					for (m=n+1; m<2*L; m++)
					{
						if (tstbit(CI_ph_basis[i],m))
						{
							for (l=L+Nd; l<2*L-1; l++)
							{
								if ( !tstbit(CI_ph_basis[i],l) || l==m || l==n)
								{
									for (k=l+1; k<2*L; k++)
									{
										if ( !tstbit(CI_ph_basis[i],k) || k==m || k==n)
										{
											SD_tmp=CI_ph_basis[i];
											
											V_tmp=elem_V(k,n,l,m);
										
											if (V_tmp)
											{
												if (l!=n)
												{
													q=0;
													if (l!=m)
													{
														q+=add_bits_between(SD_tmp,l,m);
														setbit(SD_tmp,l);
														clrbit(SD_tmp,m);
													}
													if (k!=n)
													{
														q+=add_bits_between(SD_tmp,k,n);
														setbit(SD_tmp,k);
														clrbit(SD_tmp,n);
													}
												}
												else
												{
													q=1;
													if (k!=m)
													{
														q+=add_bits_between(SD_tmp,k,m);
														setbit(SD_tmp,k);
														clrbit(SD_tmp,m);
													}
												}
												sgn=1;
												if (q%2) sgn=-1;
												V_tmp=sgn*V_tmp;
												if (SD_tmp!=CI_ph_basis[i]) add_H_m_element_sym(i, SD_tmp, V_tmp);
												else add_H_m_element_diag(i, V_tmp);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		//one spin-up and one spin-down particles
		if (nu>0 && nd>0)
		{
			for (n=L+Nd; n<2*L; n++)
			{
				if (tstbit(CI_ph_basis[i],n))
				{
					for (m=Nu; m<L; m++)
					{
						if (tstbit(CI_ph_basis[i],m))
						{
							for (l=Nu; l<L; l++)
							{
								if ( !tstbit(CI_ph_basis[i],l) || l==m)
								{
									for (k=L+Nd; k<2*L; k++)
									{
										if ( !tstbit(CI_ph_basis[i],k) || k==n)
										{
											SD_tmp=CI_ph_basis[i];
											
											V_tmp=elem_V(k,n,l,m);
										
											if (V_tmp)
											{
												q=0;
												if (l!=m)
												{
													q+=add_bits_between(SD_tmp,l,m);
													setbit(SD_tmp,l);
													clrbit(SD_tmp,m);
												}
												if (k!=n)
												{
													q+=add_bits_between(SD_tmp,k,n);
													setbit(SD_tmp,k);
													clrbit(SD_tmp,n);
												}
												sgn=1;
												if (q%2) sgn=-1;
												V_tmp=sgn*V_tmp;
												if (SD_tmp!=CI_ph_basis[i]) add_H_m_element_sym(i, SD_tmp, V_tmp);
												else add_H_m_element_diag(i, V_tmp);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		
		//hopping of two holes
		
		//two spin-down holes
		if (nu>1)
		{
			for (n=0; n<Nu-1; n++)
			{
				if (tstbit(CI_ph_basis[i],n))
				{
					for (m=n+1; m<Nu; m++)
					{
						if (tstbit(CI_ph_basis[i],m))
						{
							for (l=0; l<Nu-1; l++)
							{
								if ( !tstbit(CI_ph_basis[i],l) || l==m || l==n)
								{
									for (k=l+1; k<Nu; k++)
									{
										if ( !tstbit(CI_ph_basis[i],k) || k==m || k==n)
										{
											SD_tmp=CI_ph_basis[i];
											
											V_tmp=conj(elem_V(k,n,l,m));
										
											if (V_tmp)
											{
												if (l!=n)
												{
													q=0;
													if (l!=m)
													{
														q+=add_bits_between(SD_tmp,l,m);
														setbit(SD_tmp,l);
														clrbit(SD_tmp,m);
													}
													if (k!=n)
													{
														q+=add_bits_between(SD_tmp,k,n);
														setbit(SD_tmp,k);
														clrbit(SD_tmp,n);
													}
												}
												else
												{
													q=1;
													if (k!=m)
													{
														q+=add_bits_between(SD_tmp,k,m);
														setbit(SD_tmp,k);
														clrbit(SD_tmp,m);
													}
												}
												sgn=1;
												if (q%2) sgn=-1;
												V_tmp=sgn*V_tmp;
												if (SD_tmp!=CI_ph_basis[i]) add_H_m_element_sym(i, SD_tmp, V_tmp);
												else add_H_m_element_diag(i, V_tmp);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		//two spin-up holes
		if (nd>1)
		{
			for (n=L; n<L+Nd-1; n++)
			{
				if (tstbit(CI_ph_basis[i],n))
				{
					for (m=n+1; m<L+Nd; m++)
					{
						if (tstbit(CI_ph_basis[i],m))
						{
							for (l=L; l<L+Nd-1; l++)
							{
								if ( !tstbit(CI_ph_basis[i],l) || l==m || l==n)
								{
									for (k=l+1; k<L+Nd; k++)
									{
										if ( !tstbit(CI_ph_basis[i],k) || k==m || k==n)
										{
											SD_tmp=CI_ph_basis[i];
											
											V_tmp=conj(elem_V(k,n,l,m));
										
											if (V_tmp)
											{
												if (l!=n)
												{
													q=0;
													if (l!=m)
													{
														q+=add_bits_between(SD_tmp,l,m);
														setbit(SD_tmp,l);
														clrbit(SD_tmp,m);
													}
													if (k!=n)
													{
														q+=add_bits_between(SD_tmp,k,n);
														setbit(SD_tmp,k);
														clrbit(SD_tmp,n);
													}
												}
												else
												{
													q=1;
													if (k!=m)
													{
														q+=add_bits_between(SD_tmp,k,m);
														setbit(SD_tmp,k);
														clrbit(SD_tmp,m);
													}
												}
												sgn=1;
												if (q%2) sgn=-1;
												V_tmp=sgn*V_tmp;
												if (SD_tmp!=CI_ph_basis[i]) add_H_m_element_sym(i, SD_tmp, V_tmp);
												else add_H_m_element_diag(i, V_tmp);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		//one spin-down and one spin-up holes
		if (nu>0 && nd>0)
		{
			for (n=0; n<Nu; n++)
			{
				if (tstbit(CI_ph_basis[i],n))
				{
					for (m=L; m<L+Nd; m++)
					{
						if (tstbit(CI_ph_basis[i],m))
						{
							for (l=L; l<L+Nd; l++)
							{
								if ( !tstbit(CI_ph_basis[i],l) || l==m )
								{
									for (k=0; k<Nu; k++)
									{
										if ( !tstbit(CI_ph_basis[i],k) || k==n)
										{
											SD_tmp=CI_ph_basis[i];
											
											V_tmp=conj(elem_V(k,n,l,m));
										
											if (V_tmp)
											{
												q=0;
												if (l!=m)
												{
													q+=add_bits_between(SD_tmp,l,m);
													setbit(SD_tmp,l);
													clrbit(SD_tmp,m);
												}
												if (k!=n)
												{
													q+=add_bits_between(SD_tmp,k,n);
													setbit(SD_tmp,k);
													clrbit(SD_tmp,n);
												}
												sgn=1;
												if (q%2) sgn=-1;
												V_tmp=sgn*V_tmp;
												if (SD_tmp!=CI_ph_basis[i]) add_H_m_element_sym(i, SD_tmp, V_tmp);
												else add_H_m_element_diag(i, V_tmp);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		
		//hopping of one particle and one hole
		
		//one spin-up particle and one spin-down hole
		if (nu>0)
		{
			for (n=Nu; n<L; n++)
			{
				if (tstbit(CI_ph_basis[i],n))
				{
					for (m=0; m<Nu; m++)
					{
						if (tstbit(CI_ph_basis[i],m))
						{
							for (l=0; l<Nu; l++)
							{
								if ( !tstbit(CI_ph_basis[i],l) || l==m )
								{
									for (k=Nu; k<L; k++)
									{
										if ( !tstbit(CI_ph_basis[i],k) || k==n)
										{
											SD_tmp=CI_ph_basis[i];
											
											V_tmp=-elem_V(k,n,m,l);
										
											if (V_tmp)
											{
												q=0;
												if (l!=m)
												{
													q+=add_bits_between(SD_tmp,l,m);
													setbit(SD_tmp,l);
													clrbit(SD_tmp,m);
												}
												if (k!=n)
												{
													q+=add_bits_between(SD_tmp,k,n);
													setbit(SD_tmp,k);
													clrbit(SD_tmp,n);
												}
												sgn=1;
												if (q%2) sgn=-1;
												V_tmp=sgn*V_tmp;
												if (SD_tmp!=CI_ph_basis[i]) add_H_m_element_sym(i, SD_tmp, V_tmp);
												else add_H_m_element_diag(i, V_tmp);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		//one spin-down particle and one spin-up hole
		if (nd>0)
		{
			for (n=L+Nd; n<2*L; n++)
			{
				if (tstbit(CI_ph_basis[i],n))
				{
					for (m=L; m<L+Nd; m++)
					{
						if (tstbit(CI_ph_basis[i],m))
						{
							for (l=L; l<L+Nd; l++)
							{
								if ( !tstbit(CI_ph_basis[i],l) || l==m )
								{
									for (k=L+Nd; k<2*L; k++)
									{
										if ( !tstbit(CI_ph_basis[i],k) || k==n)
										{
											SD_tmp=CI_ph_basis[i];
											
											V_tmp=-elem_V(k,n,m,l);
										
											if (V_tmp)
											{
												q=0;
												if (l!=m)
												{
													q+=add_bits_between(SD_tmp,l,m);
													setbit(SD_tmp,l);
													clrbit(SD_tmp,m);
												}
												if (k!=n)
												{
													q+=add_bits_between(SD_tmp,k,n);
													setbit(SD_tmp,k);
													clrbit(SD_tmp,n);
												}
												sgn=1;
												if (q%2) sgn=-1;
												V_tmp=sgn*V_tmp;
												if (SD_tmp!=CI_ph_basis[i]) add_H_m_element_sym(i, SD_tmp, V_tmp);
												else add_H_m_element_diag(i, V_tmp);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		//one spin-up particle and one spin-up hole
		if (nu>0 && nd>0)
		{
			for (n=Nu; n<L; n++)
			{
				if (tstbit(CI_ph_basis[i],n))
				{
					for (m=L; m<L+Nd; m++)
					{
						if (tstbit(CI_ph_basis[i],m))
						{
							for (l=L; l<L+Nd; l++)
							{
								if ( !tstbit(CI_ph_basis[i],l) || l==m )
								{
									for (k=Nu; k<L; k++)
									{
										if ( !tstbit(CI_ph_basis[i],k) || k==n)
										{
											SD_tmp=CI_ph_basis[i];
											
											V_tmp=-elem_V(k,n,m,l);
										
											if (V_tmp)
											{
												q=0;
												if (l!=m)
												{
													q+=add_bits_between(SD_tmp,l,m);
													setbit(SD_tmp,l);
													clrbit(SD_tmp,m);
												}
												if (k!=n)
												{
													q+=add_bits_between(SD_tmp,k,n);
													setbit(SD_tmp,k);
													clrbit(SD_tmp,n);
												}
												sgn=1;
												if (q%2) sgn=-1;
												V_tmp=sgn*V_tmp;
												if (SD_tmp!=CI_ph_basis[i]) add_H_m_element_sym(i, SD_tmp, V_tmp);
												else add_H_m_element_diag(i, V_tmp);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		//one spin-down particle and one spin-down hole
		if (nu>0 && nd>0)
		{
			for (n=L+Nd; n<2*L; n++)
			{
				if (tstbit(CI_ph_basis[i],n))
				{
					for (m=0; m<Nu; m++)
					{
						if (tstbit(CI_ph_basis[i],m))
						{
							for (l=0; l<Nu; l++)
							{
								if ( !tstbit(CI_ph_basis[i],l) || l==m )
								{
									for (k=L+Nd; k<2*L; k++)
									{
										if ( !tstbit(CI_ph_basis[i],k) || k==n)
										{
											SD_tmp=CI_ph_basis[i];
											
											V_tmp=-elem_V(k,n,m,l);
										
											if (V_tmp)
											{
												q=0;
												if (l!=m)
												{
													q+=add_bits_between(SD_tmp,l,m);
													setbit(SD_tmp,l);
													clrbit(SD_tmp,m);
												}
												if (k!=n)
												{
													q+=add_bits_between(SD_tmp,k,n);
													setbit(SD_tmp,k);
													clrbit(SD_tmp,n);
												}
												sgn=1;
												if (q%2) sgn=-1;
												V_tmp=sgn*V_tmp;
												if (SD_tmp!=CI_ph_basis[i]) add_H_m_element_sym(i, SD_tmp, V_tmp);
												else add_H_m_element_diag(i, V_tmp);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
	}//end of loop over basis SD's
	
	cout<<"Hamiltonian matrix defined\n";
	
	convert_H_m_to_H_v();
}
*/
/*
template<class ind_T, class val_T, class SD_T, class cfs_T>
void Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::create_H_matrix_ph()
{
	cout<<"Defining Hamiltonian matrix...\n";
	
	compute_ref_diag_E();
	
	H_m.clear();
	H_m.resize(NSD);
	
	SD_T SD_tmp;
	ind_T i;
	int k, nu, nd; //, nu_h, nu_p, nd_h, nd_p;
	
	
	vector<int> occ_u_h;  //vector of positions of holes in the spin-up part of the SD
	vector<int> free_u_h; //vector of positions of free hole orbitals in the spin-up part of the SD
	vector<int> occ_u_p; //vector of positions of particles in the spin-up part of the SD
	vector<int> free_u_p; //vector of positions of free particle orbitals in the spin-up part of the SD
	
	vector<int> occ_d_h; //vector of positions of holes in the spin-down part of the SD
	vector<int> free_d_h; //vector of positions of free hole orbitals in the spin-down part of the SD
	vector<int> occ_d_p; //vector of positions of particles in the spin-down part of the SD
	vector<int> free_d_p; //vector of positions of free particle orbitals in the spin-down part of the SD
	
	//for the reference SD
	i=0;
	
	for (k=0; k<Nu; k++)
	{
		free_u_h.push_back(k);
	}
	for (k=Nu; k<L; k++)
	{
		free_u_p.push_back(k);
	}
	for (k=L; k<L+Nd; k++)
	{
		free_d_h.push_back(k);
	}
	for (k=L+Nd; k<2*L; k++)
	{
		free_d_p.push_back(k);
	}
	
	//creation of a single particle-hole pair by the one-body term
	
	//up-down pair
	K_terms(i,free_u_p,free_u_h);
	
	//down-up pair
	K_terms(i,free_d_p,free_d_h);
	
	
	//Creation of two particle-hole pairs by the two-body term V_knlm(p^+_k)(p^+_l)(h^+_m)(h^+_n)
	//Each distinct combination appears only once, so that there is no prefactor.
	
	//two up-down pairs
	V_terms(i, free_u_p, free_u_p, free_u_h, free_u_h);
	
	//two down-up pairs
	V_terms(i, free_d_p, free_d_p, free_d_h, free_d_h);

	//cone up-down and one down-up pairs
	V_terms(i, free_u_p, free_d_p, free_d_h, free_u_h);
	
	auto H_it=H_m[0].begin();
	
	//for all other SD's
	for (i=1; i<NSD; i++)
	{
		
		SD_tmp=CI_ph_basis[i];
		
		nu=add_bits(SD_tmp, 0, L-1)/2;  //number of excitations in the spin-up part of the SD
		nd=add_bits(SD_tmp, L, 2*L-1)/2; //number of excitations in the spin-down part
//		nu_h=add_bits(SD_tmp, 0, Nu-1); //number of holes in the spin-up part (spin-down holes)
//		nu_p=add_bits(SD_tmp, Nu, L-1); //number of particles in the spin-up part
//		nd_h=add_bits(SD_tmp, L, L+Nd-1); //number of holes in the spin-down part (spin-up holes)
//		nd_p=add_bits(SD_tmp, L+Nd, 2*L-1); //number of particles in the spin-down part
		
//		print_binary_string(SD_tmp,2*L);
//		cout<<endl;
//		cout<<"nd, nu: "<<setw(10)<<nd<<nu<<endl;
//		cout<<"nd_p, nd_h, nu_p, nu_h: "<<setw(10)<<nd_p<<setw(10)<<nd_h<<setw(10)<<nu_p<<setw(10)<<nu_h<<endl;
		
		occ_u_h.clear();
		free_u_h.clear();
		occ_u_p.clear();
		free_u_p.clear();
		
		occ_d_h.clear();
		free_d_h.clear();
		occ_d_p.clear();
		free_d_p.clear();
		
		for (k=0; k<Nu; k++)
		{
			if (tstbit(SD_tmp,k)) occ_u_h.push_back(k);
			else free_u_h.push_back(k);
		}
		for (k=Nu; k<L; k++)
		{
			if (tstbit(SD_tmp,k)) occ_u_p.push_back(k);
			else free_u_p.push_back(k);
		}
		for (k=L; k<L+Nd; k++)
		{
			if (tstbit(SD_tmp,k)) occ_d_h.push_back(k);
			else free_d_h.push_back(k);
		}
		for (k=L+Nd; k<2*L; k++)
		{
			if (tstbit(SD_tmp,k)) occ_d_p.push_back(k);
			else free_d_p.push_back(k);
		}
		
		//creation of a single particle-hole pair by the one-body term
		
		//up-down pair
		if (nu<Nph_u_max) K_terms(i,free_u_p,free_u_h);
		
		//down-up pair
		if (nd<Nph_d_max) K_terms(i,free_d_p,free_d_h);
		
		
		//hopping and onsite energy
		
		//spin-down hole
		if (nu>0) K_terms(i,free_u_h,occ_u_h);
		
		//spin-up particle
		if (nu>0) K_terms(i,free_u_p,occ_u_p);
		
		//spin-up hole
		if (nd>0) K_terms(i,free_d_h,occ_d_h);
		
		//spin-down particle
		if (nd>0) K_terms(i,free_d_p,occ_d_p);
		
	
		//creation of two particle-hole pairs by the two-body term V_knlm(p^+_k)(p^+_l)(h^+_m)(h^+_n)=V_knlm(p^+_k)(h^+_n)(p^+_l)(h^+_m)
		//Each distinct combination appears only once, so that there is no prefactor.
		
		//two up-down pairs
		if (nu<Nph_u_max-1) V_terms(i, free_u_p, free_u_p, free_u_h, free_u_h);
		
		//two down-up pairs
		if (nd<Nph_d_max-1) V_terms(i, free_d_p, free_d_p, free_d_h, free_d_h);
		
		//one up-down and one down-up pairs
		if (nd<Nph_d_max && nu<Nph_u_max) V_terms(i, free_u_p, free_d_p, free_d_h, free_u_h);
		
	
		//creation of one particle-hole pair and hopping of a particle
		
		//up-down pair and hopping of a spin-up particle
		if (nu>0 && nu<Nph_u_max) V_terms(i, free_u_p, free_u_p, free_u_h, occ_u_p);
		
		//up-down particle-hole pair and hopping of a spin-down particle
		if (nd>0 && nu<Nph_u_max) V_terms(i, free_d_p, free_u_p, free_u_h, occ_d_p);
		
		//one down-up particle-hole pair and hopping of a spin-down particle
		if (nd>0 && nd<Nph_d_max) V_terms(i, free_d_p, free_d_p, free_d_h, occ_d_p);
		
		//one down-up particle-hole pair and hopping of a spin-up particle
		if (nu>0 && nd<Nph_d_max) V_terms(i, free_u_p, free_d_p, free_d_h, occ_u_p);
		
	
		//creation of one particle-hole pair and hopping of a hole
		
		//up-down pair and hopping of a spin-down hole
		if (nu>0 && nu<Nph_u_max) V_terms(i, free_u_h, free_u_p, free_u_h, occ_u_h);
		
		//up-down pair and hopping of a spin-up hole
		if (nd>0 && nu<Nph_u_max) V_terms(i, free_d_h, free_u_p, free_u_h, occ_d_h);
		
		//down-up pair and hopping of a spin-up hole
		if (nd>0 && nd<Nph_d_max) V_terms(i, free_d_h, free_d_p, free_d_h, occ_d_h);
		
		//down-up pair and hopping of a spin-down hole
		if (nu>0 && nd<Nph_d_max) V_terms(i, free_u_h, free_d_p, free_d_h, occ_u_h);
		
	
		//hopping of two particles
	
		//two spin-up particles
		if (nu>1) V_terms(i, free_u_p, free_u_p, occ_u_p, occ_u_p);
	
		//two spin-down particles
		if (nd>1) V_terms(i, free_d_p, free_d_p, occ_d_p, occ_d_p);
		
		//one spin-up and one spin-down particles
		if (nu>0 && nd>0) V_terms(i, free_u_p, free_d_p, occ_d_p, occ_u_p);
		
	
		//hopping of two holes
		
		//two spin-down holes
		if (nu>1) V_terms(i, free_u_h, free_u_h, occ_u_h, occ_u_h);
		
		//two spin-up holes
		if (nd>1) V_terms(i, free_d_h, free_d_h, occ_d_h, occ_d_h);
		
		//one spin-down and one spin-up holes
		if (nu>0 && nd>0) V_terms(i, free_u_h, free_d_h, occ_d_h, occ_u_h);
		
		
		//hopping of one particle and one hole
		
		//spin-up particle and spin-down hole
		if (nu>0 && nu>0) V_terms(i, free_u_p, free_u_h, occ_u_h, occ_u_p);
		
		//spin-up particle and spin-up hole
		if (nu>0 && nd>0) V_terms(i, free_u_p, free_d_h, occ_d_h, occ_u_p);
		
		//spin-down particle and spin-down hole
		if (nd>0 && nu>0) V_terms(i, free_d_p, free_u_h, occ_u_h, occ_d_p);
		
		//spin-down particle and spin-up hole
		if (nd>0 && nd>0) V_terms(i, free_d_p, free_d_h, occ_d_h, occ_d_p);
		
	}//end of loop over basis SD's
	
	cout<<"Hamiltonian matrix defined\n";
	
	convert_H_m_to_H_v();
}
*/

template<class ind_T, class val_T, class SD_T, class cfs_T>
void Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::check_H_Hermicity()
{
	ind_T i, j;
	auto it1=H_m[0].begin();
	auto it2=H_m[0].begin();
	bool H_Hermitian=true;
	
	for (i=1; i<NSD; i++)
	{
		for (it1=H_m[i].begin(); it1!=H_m[i].end(); it1++)
		{
			j=it1->first;
			it2=H_m[j].find(i);
			if (it2->second!=conj(it1->second))
			{
				H_Hermitian=false;
				i=NSD;
				break;
			}
		}
	}
	
	if (H_Hermitian) cout<<"Hamiltonian matrix is Hermitian\n";
	else cout<<"Hamiltonian matrix is not Hermitian\n";
}

template<class ind_T, class val_T, class SD_T, class cfs_T>
void Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::K_terms(int i, const vector<int> &v1, const vector<int> &v2)
{
	bool hopping=false;
	
	if (tstbit(CI_ph_basis[i],v2[0])) hopping=true;
	
	int a,b,k,l;
	
	for (b=0; b<v2.size(); b++)
	{
		l=v2[b];
		
		if (hopping) add_K_term_H_m(i, l, l);
		
		for (a=0; a<v1.size(); a++)
		{
			k=v1[a];
			
			add_K_term_H_m(i, k, l);
		}
		
	}
}

//add kinetic energy term to H_m matrix. If k!=l, the transposed matrix element is also added so this function must be called only once for each pair of spin-orbital indices
template<class ind_T, class val_T, class SD_T, class cfs_T>
void Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::add_K_term_H_m(int i, int k, int l)
{
	val_T t_tmp, sgn;
	int q;
	
	SD_T SD_tmp=CI_ph_basis[i];
	
	operator_type ok, ol;
	
	if ( k<Nu ||  (k>=L && k<L+Nd) ) ok=hole;
	else ok=part;
	
	if ( l<Nu ||  (l>=L && l<L+Nd) ) ol=hole;
	else ol=part;
	
	if (k<L && l<L) t_tmp=SO_H_T::t_Fock_u[ph_indices[k]][ph_indices[l]];
	else if (k>=L && l>=L) t_tmp=SO_H_T::t_Fock_d[ph_indices[k]-L][ph_indices[l]-L];
	else cerr<<"add_K_term() error: spins must be equal\n";
	
	if (t_tmp)
	{
		if (ok==hole && ol==hole) t_tmp=-conj(t_tmp);
		if (l!=k)
		{
			sgn=1;
			q=add_bits_between(SD_tmp, k, l);
			setbit(SD_tmp,k);
			if (!tstbit(CI_ph_basis[i],l))  setbit(SD_tmp,l);
			else clrbit(SD_tmp,l);
			if (q%2) sgn=-1;
			t_tmp=sgn*t_tmp;
			add_H_m_element_sym(i, SD_tmp, t_tmp);
		}
		else add_H_m_element_diag(i, t_tmp);
	}
}

template<class ind_T, class val_T, class SD_T, class cfs_T>
void Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::V_terms(int i, const vector<int> &v1, const vector<int> &v2, const vector<int> &v3, const vector<int> &v4)
{
	bool one_hopping=false, two_hoppings=false;
	
	if (tstbit(CI_ph_basis[i],v4[0])) one_hopping=true;
	if (tstbit(CI_ph_basis[i],v3[0])) two_hoppings=true;
	
	int a,b,c,d,k,l,m,n;
	int c0, a0;
	
	for (d=0; d<v4.size(); d++)
	{
		n=v4[d];
		
		c0=0;
		if (&v3==&v4) c0=d+1;
		for (c=c0; c<v3.size(); c++)
		{
			m=v3[c];
		
			if (two_hoppings) //if both m and n are occupied orbitals, add the density-density terms
			{
				add_V_term_H_m(i, n, m, m, n);
			}
		
			for (b=0; b<v2.size(); b++)
			{
				l=v2[b];
			
				a0=0;
				if (&v1==&v2) a0=b+1;
				for (a=a0; a<v1.size(); a++)
				{
					k=v1[a];
					
					add_V_term_H_m(i, k, l, m, n);
				}
			}
		}
	}
	
	//correlated hopping and correlated pair creation terms
	if (one_hopping)
	{
		if (&v3!=&v4)
		{
			//correlated hopping of the form Vkklm(p+_lp_m p+_kp_k) or correlated pair creation of the form Vkklm(p+_lh+_m p+_kp_k)
 			for (d=0; d<v4.size(); d++)
			{
				n=v4[d];
				k=n;
				for (c=0; c<v3.size(); c++)
				{
					m=v3[c];
					for (b=0; b<v2.size(); b++)
					{
						l=v2[b];
						add_V_term_H_m(i, k, l, m, n);
					}
					
				}
			}
			if (two_hoppings) //if both v3 and v4 are occupied spin-orbitals (and are not the same), add the correlated hopping terms of the form Vknll(p+_kp_n p+_lp_l)
			{
				for (c=0; c<v3.size(); c++)
				{
					m=v3[c];
					l=m;
					for (d=0; d<v4.size(); d++)
					{
						n=v4[d];
						for (a=0; a<v1.size(); a++)
						{
							k=v1[a];
							
							add_V_term_H_m(i, k, l, m, n);
						}
					}
				}
			}
		}
		else  //add the correlated hopping terms of the form Vkklm(p+_lp_m p+_kp_k) where k goes through all the occupied spin-orbitals and m and l take all the other possible values
		{
			for (d=0; d<v4.size(); d++)
			{
				n=v4[d];
				k=n;
				for (c=0; c<d; c++)
				{
					m=v3[c];
					for (b=0; b<v2.size(); b++)
					{
						l=v2[b];
						add_V_term_H_m(i, k, l, m, n);
					}
				}
				for (c=d+1; c<v3.size(); c++)
				{
					m=v3[c];
					for (b=0; b<v2.size(); b++)
					{
						l=v2[b];
						add_V_term_H_m(i, k, l, m, n);
					}
				}
			}
		}
	}
}

//add potential energy term to H_m matrix. If the term creates at least one pair, the transposed matrix element is also added
template<class ind_T, class val_T, class SD_T, class cfs_T>
void Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::add_V_term_H_m(int i, int k, int l, int m, int n)
{
	SD_T SD_tmp=CI_ph_basis[i];

	operator_type ok, ol, om, on;
	
	if ( k<Nu ||  (k>=L && k<L+Nd) ) ok=hole;
	else ok=part;
	
	if ( l<Nu ||  (l>=L && l<L+Nd) ) ol=hole;
	else ol=part;
	
	if ( m<Nu ||  (m>=L && m<L+Nd)) om=hole;
	else om=part;
	
	if ( n<Nu ||  (n>=L && n<L+Nd)) on=hole;
	else on=part;
	
	int ind=ph_indices[k]+ph_indices[n]*s[1]+ph_indices[l]*s[2]+ph_indices[m]*s[3];
	if ( (ok==hole && ol==part) || (ok==part && ol==hole))
		ind=ph_indices[k]+ph_indices[n]*s[1]+ph_indices[m]*s[2]+ph_indices[l]*s[3];
	
	if ((k<L && n>=L) || (k>=L && n<L)) cerr<<"add_V_term_H_m() error: first and last spins must be equal\n";
	if ((l<L && m>=L) || (l>=L && m<L)) cerr<<"add_V_term_H_m() error: second and third spins must be equal\n";
	
	val_T V_tmp=0, sgn;
	auto V_it=SO_H_T::V_u_SO_basis.begin();
	if (k<L && l<L)
	{
		V_it=SO_H_T::V_u_SO_basis.find(ind);
		if (V_it!=SO_H_T::V_u_SO_basis.end()) V_tmp=V_it->second;
	}
	else if (k>L && l>L)
	{
		V_it=SO_H_T::V_d_SO_basis.find(ind);
		if (V_it!=SO_H_T::V_d_SO_basis.end()) V_tmp=V_it->second;
	}
	else
	{
		V_it=SO_H_T::V_ud_SO_basis.find(ind);
		if (V_it!=SO_H_T::V_ud_SO_basis.end()) V_tmp=V_it->second;
	}
	
	int q;
	if (V_tmp)
	{
		if ( (ok==hole || ol==hole) && om==hole && on==hole) V_tmp=conj(V_tmp);
		if ( (ok==hole && ol==part) || (ok==part && ol==hole)) V_tmp=-V_tmp;
		sgn=1;
		if (n!=l)
		{
			q=0;
			if (l!=m)
			{
				q+=add_bits_between(SD_tmp, l, m);
				setbit(SD_tmp,l);
				if (!tstbit(SD_tmp,m)) setbit(SD_tmp,m);
				else clrbit(SD_tmp,m);
			}
			if (n!=k)
			{
				q+=add_bits_between(SD_tmp, k, n);
				setbit(SD_tmp,k);
				if (!tstbit(SD_tmp,n)) setbit(SD_tmp,n);
				else clrbit(SD_tmp,n);
			}
		}
		else
		{
			q=1;
			if (m!=k)
			{
				q+=add_bits_between(SD_tmp, k, m);
				setbit(SD_tmp,k);
				if (!tstbit(SD_tmp,m)) setbit(SD_tmp,m);
				else clrbit(SD_tmp,m);
			}
		}
		if (q%2) sgn=-1;
		V_tmp=sgn*V_tmp;
		if (SD_tmp!=CI_ph_basis[i]) add_H_m_element_sym(i, SD_tmp, V_tmp);
		else add_H_m_element_diag(i, V_tmp);
	}
}

template<class ind_T, class val_T, class SD_T, class cfs_T>
int Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::add_bits_between(SD_T SD_tmp, int b1, int b2)
{
	int s=std::min(b1,b2);
	int g=std::max(b1,b2);
	
	int n=0;
	if (g>s+1) n=add_bits(SD_tmp, s+1, g-1);
	
	return n;
}

template<class ind_T, class val_T, class SD_T, class cfs_T>
ind_T Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::get_SD_index(SD_T sd)
{
	auto it=CI_ph_basis_map.find(sd);
	if (it!=CI_ph_basis_map.end()) return it->second;
	else
	{
		cerr<<"get_SD_index() error\n";
		exit(EXIT_FAILURE);
	}
}
/*
template<class ind_T, class val_T, class SD_T, class cfs_T>
ind_T Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::get_SD_index(SD_T sd)
{
	unsigned nd, nu, n, nd_min;
	
	nu=add_bits(sd, 0, L-1)/2;
	nd=add_bits(sd, L, 2*L-1)/2;
	n=nu+nd;
	
	ind_T N_u_h=binom_coef(Nu,nu);
	ind_T N_u_p=binom_coef(L-Nu,nu);
	ind_T N_d_h=binom_coef(Nd,nd);
	
	nd_min=std::max(0,(int)n-Nph_u_max);
	
	ind_T ref_ind=N_SD_l[n][nd-nd_min];
	
	SD_T str=sd;
	ind_T i,j,k,l;
	
	i=get_string_index(str,Nu);
	str=str>>Nu;
	j=get_string_index(str,L-Nu);
	str=str>>L-Nu;
	k=get_string_index(str,Nd);
	str=str>>Nd;
	l=get_string_index(str,L-Nd);
	
	return ref_ind+i+N_u_h*(j+N_u_p*(k+N_d_h*l));
}
*/

template<class ind_T, class val_T, class SD_T, class cfs_T>
ind_T Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::get_string_index(SD_T str, unsigned Nb)
{
	ind_T ind=0;
	
	int i=Nb-1;
	unsigned s;
	
	while (i>=0 && !tstbit(str,i)) i--;
	if (i<0) return 0;
	
	s=1;
	i--;
	while (i>=0)
	{
		if (tstbit(str,i)) s++;
		else ind+=binom_coef(Nb-i-1, s-1);
		i--;
	}
	
	return ind;
}

template<class ind_T, class val_T, class SD_T, class cfs_T>
void Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::create_list_ph(ind_T Nph_max, unsigned Ns, vector<vector<SD_T>> &list_ph)
{
	int i;
	ind_T j, k, Nst_h, Nst_p, Nst_ph;
	SD_T st_tmp;
	vector<SD_T> list_h, list_p;
	
	list_ph.resize(Nph_max+1);
	list_ph[0].resize(1);
	list_ph[0][0]=0;
	for (i=1; i<=Nph_max; i++)
	{
		//generate all the hole strings with the numbers of set bits (excitations) i
		//after the call to generate_states(), list_h is a vector of bit strings of type SD_T
		generate_states(list_h, Ns, i);
		Nst_h=binom_coef(Ns,i);
		
		//generate all the particle strings with the numbers of set bits (excitations) i
		//after the call to generate_states(), list_p is a vector of bit strings of type SD_T
		if ((L-Ns)!=Ns)
		{
			generate_states(list_p, L-Ns, i);
			Nst_p=binom_coef(L-Ns,i);
		}
		else
		{
			list_p=list_h;
			Nst_p=Nst_h;
		}
		
		//check that list_h and list_p have the correct sizes
		if (Nst_h!=list_h.size())
		{
			cerr<<"create_list_ph() Nst_h error\n";
		}
		if (Nst_p!=list_p.size())
		{
			cerr<<"create_list_ph() Nst_p error\n";
		}
		
		//number of different spin-up particle-hole bit strings for i excitations
		Nst_ph=Nst_h*Nst_p;
		
		//fill list_ph[i] with the Nst_ph combinations of hole and particle strings for i particle-hole excitations
		list_ph[i].resize(Nst_ph);
		for (j=0; j<Nst_p; j++)
		{
			st_tmp=list_p[j];
			//shift the particle string left by Ns bits
			st_tmp=st_tmp<<Ns;
			for (k=0; k<Nst_h; k++)
			{
				list_ph[i][k+j*Nst_h]=list_h[k]+st_tmp;
			}
		}
	}
	
}

template<class ind_T, class val_T, class SD_T, class cfs_T>
bool Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::create_CI_ph_basis()
{
//	define_SO_ph_map();
	
	unsigned long NSD_u=binom_coef(L, Nu);
	unsigned long NSD_d=binom_coef(L, Nd);
	
	NSD=(ind_T)NSD_u*NSD_d;
	
	cout<<"total number of SD's in the basis: "<<NSD<<endl;
	
	//determine the maximum numbers of excitations for each spin and total
	Nph_u_max=Nu;
	Nph_d_max=Nd;
	if ((L-Nu)<Nu) Nph_u_max=L-Nu;
	if ((L-Nd)<Nd) Nph_d_max=L-Nd;
	int Nph_max=Nph_d_max+Nph_u_max;
	
	//prepare the basis vector for insertion of the basis SD's and insert the reference SD (0 excitation SD)
	CI_ph_basis_map.clear();
	CI_ph_basis.clear();
	CI_ph_basis.resize(NSD);
	CI_abs_basis.clear();
	CI_abs_basis.resize(NSD);
	abs_to_ph_sign.clear();
	abs_to_ph_sign.resize(NSD);
	
	CI_ph_basis[0]=0;
	CI_abs_basis[0]=ref_SD;
	abs_to_ph_sign[0]=1;
	ind_T Nst_sum=1;
	
	CI_ph_basis_map.insert(pair<SD_T,ind_T>(0,0));
	auto it=CI_ph_basis_map.begin();
	
	//the bit strings required to generate the SD's. list_h and list_p are vectors of hole and particle strings
	//list_ph_u and list_ph_d are vectors in which each element is a vector of strings with a given number of excitations
	vector<vector<SD_T>> list_ph_u, list_ph_d;
	
	create_list_ph(Nph_u_max, Nu, list_ph_u);
	if (Nd==Nu) list_ph_d=list_ph_u;
	else create_list_ph(Nph_d_max, Nd, list_ph_d);
	
	N_SD_l.resize(Nph_max+1);
	N_SD_l[0].resize(1);
	N_SD_l[0][0]=0;
	
//	cout<<setiosflags(ios::left);
	
	int nd_min, nd_max, n, nd;
	ind_T j, k, Nst_ph_u, Nst, m;
	SD_T st_tmp, SD_abs;
	int sgn;

//	cout<<0<<"\t\t";
//	print_SD(CI_ph_basis[0]);
//	cout<<"\t\t";
//	print_SD(CI_abs_basis[0]);
//	cout<<"\t\t"<<abs_to_ph_sign[0]<<endl;
	
	//for all the numbers of excitations
	for (n=1; n<=Nph_max; n++)
	{
		//	cout<<n<<" excitations total\n";
		//set the minimum and maximum numbers of spin-down excitations
		nd_min=std::max(0,n-Nph_u_max);
		nd_max=std::min(n,Nph_d_max);
		//	cout<<"nd_min: "<<nd_min<<endl;
		//	cout<<"nd_max: "<<nd_max<<endl;
		//set the size of N_SD_l[n] to the number of different numbers of spin-down excitations
		N_SD_l[n].resize(nd_max-nd_min+1);
		for (nd=nd_min; nd<=nd_max; nd++)
		{
			Nst_ph_u=list_ph_u[n-nd].size();
			Nst=list_ph_u[nd].size()*list_ph_u[n-nd].size();
			//j is the index of a string in the spin-down particle-hole strings list, which here is list_ph_u because N_down=N_up
			for (j=0; j<list_ph_d[nd].size(); j++)
			{
				//st_tmp is the spin-down part of the SD
				st_tmp=list_ph_d[nd][j];
				//shift st_tmp toward the left by L bits
				st_tmp=st_tmp<<L;
				//k is the index of a string in the spin-up particle-hole strings list
				for (k=0; k<Nst_ph_u; k++)
				{
					m=Nst_sum+k+j*Nst_ph_u;
					CI_ph_basis[m]=list_ph_u[n-nd][k]+st_tmp;
					it=CI_ph_basis_map.find(CI_ph_basis[m]);
					if (it!=CI_ph_basis_map.end()) cout<<"create_CI_ph_basis() error\n";
					else CI_ph_basis_map.insert(pair<SD_T,ind_T>(CI_ph_basis[m],m));
					convert_from_ph_SD(CI_ph_basis[m],CI_abs_basis[m],abs_to_ph_sign[m]);
				//	cout<<m<<"\t\t";
				//	print_SD(CI_ph_basis[m]);
				//	cout<<"\t\t";
				//	print_SD(CI_abs_basis[m]);
				//	cout<<"\t\t"<<abs_to_ph_sign[m]<<endl;
				}
			}
			N_SD_l[n][nd-nd_min]=Nst_sum;
			Nst_sum=Nst_sum+list_ph_d[nd].size()*list_ph_u[n-nd].size();
			
			//				cout<<"CI_ph_basis size: "<<Nst_sum<<endl;
		}
	}
	
/*
	for (j=0; j<NSD; j++)
	{
		k=get_SD_index(CI_ph_basis[j]);
		print_SD(CI_ph_basis[j]);
		cout<<"\t\t"<<setw(10)<<j<<setw(10)<<k<<j-k;
		cout<<endl;
	}
*/
	
	return true;
	
}

template<class ind_T, class val_T, class SD_T, class cfs_T>
bool Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::generate_states(vector<SD_T> &list, unsigned Nb, unsigned Nx)
{
	unsigned i;
	
	unsigned long Nst=binom_coef(Nb,Nx);
	
	list.resize(Nst);
	
	SD_T st_init=0;
	
	for (i=0; i<Nx; i++) setbit(st_init,i);
	
	ind_T j=0;
	list[j]=st_init;
	for (j=1; j<Nst; j++)
	{
		list[j]=next_state(list[j-1],Nb);
	}
/*
 	int Nmax_print=50;
 	if (Nmax_print>Nst) Nmax_print=Nst;
 
 	cout<<"number of states with "<<Nx<<" excitations on "<<Nb<<" sites: "<<Nst<<endl;
 
	for (j=0; j<Nmax_print; j++)
	{
		print_binary_string(list[j], Nb);
		cout<<endl;
	}
*/
	return true;
}

//starting from a string with n setbits, obtain the nex larger string with n setbits
template<class ind_T, class val_T, class SD_T, class cfs_T>
SD_T Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::next_state(SD_T s, unsigned Nb)
{
	unsigned i=0, j=0, k;
	
	SD_T s_out=s;
	
	//find the first pair 01 in the string starting at the least-significant bit and count the number of set bits on the right-hand side of that pair
	while (i<Nb-1 && !(tstbit(s,i) && !tstbit(s,i+1)))
	{
		if (tstbit(s,i)) j++;
		i++;
	}// if such a pair has been found, do 01->10, and put the j set bits at the right-hand side of the pair to the least signifiacnt positions
	if (i<Nb-1)
	{
		clrbit(s_out,i);
		setbit(s_out,i+1);
		for (k=0; k<j; k++) setbit(s_out,k);
		for (k=j; k<i; k++) clrbit(s_out,k);
	}
	
	return s_out;
}

//compute the reference SD diagonal energy
template<class ind_T, class val_T, class SD_T, class cfs_T>
void Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::compute_ref_diag_E()
{
//	cout<<"compute_ref_diag_E()\n";
	
	ind_T ind,i,j,k,l,m,n;
	auto V_it=SO_H_T::V_u_SO_basis.begin();
	
	E_ref=0;
	
	//one-body part
	for (k=0; k<Nu; k++)
	{
		l=SO_H_T::holes_ind[k];
		
		E_ref+=SO_H_T::t_SO_basis_u[l][l];
	}
	for (k=Nu; k<SO_H_T::N; k++)
	{
		l=SO_H_T::holes_ind[k];
		
		E_ref+=SO_H_T::t_SO_basis_d[l-L][l-L];
	}
	
	cout<<"E_ref kinetic:"<<E_ref<<endl;
	
	//up-up two-body part
	for (i=0; i<Nu-1; i++)
	{
		k=SO_H_T::holes_ind[i];
		for (j=i+1; j<Nu; j++)
		{
			l=SO_H_T::holes_ind[j];
			ind=k+k*s[1]+l*s[2]+l*s[3];
			V_it=SO_H_T::V_u_SO_basis.find(ind);
			if (V_it!=SO_H_T::V_u_SO_basis.end())
				E_ref+=V_it->second;
		}
	}
	
	//up-down two-body part
	for (i=0; i<Nu; i++)
	{
		k=SO_H_T::holes_ind[i];
		for (j=Nu; j<SO_H_T::N; j++)
		{
			l=SO_H_T::holes_ind[j];
			ind=k+k*s[1]+l*s[2]+l*s[3];
			V_it=SO_H_T::V_ud_SO_basis.find(ind);
			if (V_it!=SO_H_T::V_ud_SO_basis.end())
				E_ref+=V_it->second;
		}
	}
	
	//down-down two-body part
	for (i=Nu; i<SO_H_T::N-1; i++)
	{
		k=SO_H_T::holes_ind[i];
		for (j=i+1; j<SO_H_T::N; j++)
		{
			l=SO_H_T::holes_ind[j];
			ind=k+k*s[1]+l*s[2]+l*s[3];
			V_it=SO_H_T::V_d_SO_basis.find(ind);
			if (V_it!=SO_H_T::V_d_SO_basis.end())
				E_ref+=V_it->second;
		}
	}
	
	cout<<"E_ref total: "<<E_ref<<endl;
	
}

template<class ind_T, class val_T, class SD_T, class cfs_T>
void Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::create_pascal_matrix()
{
	unsigned n=L+1;
	int j,l;
	Pascal_mat.resize(n);
	Pascal_mat[0].resize(n,1);
	for (j=1; j<n; j++)
	{
		Pascal_mat[j].resize(n,0);
		Pascal_mat[j][0]=1;
	}
	for (j=1; j<n-1; j++)
		for (l=1; l<n-j; l++)
			Pascal_mat[j][l]=Pascal_mat[j][l-1]+Pascal_mat[j-1][l];
/*
	cout<<setiosflags(ios::left);
	for (j=0; j<n; j++)
	{
		for (l=0; l<n; l++)
			cout<<setw(10)<<Pascal_mat[j][l];
		cout<<endl;
	}
 */
}

template<class ind_T, class val_T, class SD_T, class cfs_T>
unsigned long Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T>::binom_coef(unsigned N, unsigned n)
{
	if (N<=L) return Pascal_mat[N-n][n];
	else
	{
		cout<<"binom_coef() warning: value outside precalculated Pascal triangle\n";
		return binom_cf(N,n);
	}
}


#endif /* Hamiltonian_diag_h */
