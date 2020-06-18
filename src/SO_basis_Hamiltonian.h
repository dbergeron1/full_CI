//
//  SO_basis_Hamiltonian.h
//  full_CI
//
//  Created by Dom2 on 20-03-28.
//

#ifndef SO_basis_Hamiltonian_h
#define SO_basis_Hamiltonian_h

#include "Hamiltonian.h"

template<class ind_T, class val_T, class SD_T>
class SO_basis_Hamiltonian : public Hamiltonian<ind_T, val_T>
{
	typedef Hamiltonian<ind_T, val_T> H_T;
public:
//	SO_basis_Hamiltonian(const Hamiltonian<ind_T, val_T> &H_par, array<unsigned,2> N_par, val_T tol_p);
	SO_basis_Hamiltonian(const array<unsigned, 2> &N_par, val_T tol_p, const coord_vec_T &lattice_par, const vector<val_T> &U_par, const vector<val_T> &t_par, const vector<val_T> &onsite_E_par);
	
	vector<vector<val_T>> SO_basis_u, SO_basis_d;
//	sparse_matrix_vec<ind_T, val_T> t_SO_basis_u, t_SO_basis_d, t_Fock_u, t_Fock_d;
	vector<vector<val_T>> t_SO_basis_u, t_SO_basis_d, t_Fock_u, t_Fock_d;
	
	map<ind_T, val_T> V_u_SO_basis, V_d_SO_basis, V_ud_SO_basis;
//	map<ind_T, val_T> V_u_copy, V_d_copy, V_ud_copy;
//	map<ind_T, val_T> V_u_tr, V_d_tr, V_ud_tr;
	
	ind_T &L=H_T::L;
	array<ind_T,4> &s=H_T::s;
	
	SD_T ref_SD;
	unsigned N, Nu, Nd;
	vector<unsigned> particles_ind;
	vector<unsigned> holes_ind;
	vector<unsigned> ph_indices;
	vector<unsigned> SO_ph_map; //map from spin-orbital indices to particle and hole indices
	val_T tol_t_min, tol_V_min;
	
	void set_SO_basis(const vector<vector<val_T>> &SO_basis_u_par, const vector<vector<val_T>> &SO_basis_d_par, SD_T ref_par, val_T r=1.0);
	void compute_Hartree_Fock_sol(int Niter_max=100, val_T r_init=0.01, val_T Dr=0.01, val_T tol=1.0e-8, double f_A=1.0, double g=0.5);
	void copy_V();
	void transform_V(vector<vector<val_T>> basis_u, vector<vector<val_T>> basis_d);
	void compare_V();
	void set_ph_indices();
	void define_SO_ph_map();
	void set_t_SO_basis();
	void set_V_SO_basis();
	void set_t_Fock(val_T r);
};

//SO_basis_Hamiltonian<ind_T, val_T, SD_T>::SO_basis_Hamiltonian(const Hamiltonian<ind_T, val_T> &H_par, array<unsigned,2> N_par, val_T tol_p): Hamiltonian<ind_T, val_T>(H_par)
template<class ind_T, class val_T, class SD_T>
SO_basis_Hamiltonian<ind_T, val_T, SD_T>::SO_basis_Hamiltonian(const array<unsigned,2> &N_par, val_T tol_p, const coord_vec_T &lattice_par, const vector<val_T> &U_par, const vector<val_T> &t_par, const vector<val_T> &onsite_E_par): Hamiltonian<ind_T, val_T>(lattice_par, U_par, t_par, onsite_E_par)
{
	Nu=N_par[0];
	Nd=N_par[1];
	N=Nu+Nd;
	tol_t_min=tol_p;
	tol_V_min=tol_p;
	
	cout<<"SO_basis_Hamiltonian, Nu: "<<Nu<<endl;
	cout<<"SO_basis_Hamiltonian, Nd: "<<Nd<<endl;
	
}

//compute the Hartree-Fock solution, r_init is the initial ratio of the interaction term in the iterative algorithm and the true one, and Dr is the step in r, f_A is used to define random orthogonal matrices for the initial spin-orbitals bases and can be seen as a deviation from the identity matrix 
template<class ind_T, class val_T, class SD_T>
void SO_basis_Hamiltonian<ind_T, val_T, SD_T>::compute_Hartree_Fock_sol(int Niter_max, val_T r_init, val_T Dr, val_T tol, double f_A, double g)
{
	int i,j,k;
	
	cout<<"computing Hartree-Fock orbitals...\n";
	
	mt19937 rnd_gen;
	uniform_int_distribution<unsigned> distr;
	uniform_real_distribution<double> distr_real(0,1.0);
	
	rnd_gen.seed(time(NULL));
	
	vector<vector<val_T>> Au(L), Aut, Au_as, Ad(L), Adt, Ad_as, Id(L);
	
	if (H_T::H_real)
	{
		for (i=0; i<L; i++)
		{
			Id[i].resize(L,0);
			Id[i][i]=1.0;
			Au[i].resize(L);
			Ad[i].resize(L);
			for (j=0; j<L; j++)
			{
				Au[i][j]=f_A*(2*distr_real(rnd_gen)-1.0);
				Ad[i][j]=f_A*(2*distr_real(rnd_gen)-1.0);
			}
		}
		
		val_T one=1.0, m_one=-1.0;
		
		transpose_matrix(Au, Aut);
		transpose_matrix(Ad, Adt);
		matrix_sum(one, Au, m_one, Aut, Au_as);
		matrix_sum(one, Ad, m_one, Adt, Ad_as);
		
		vector<vector<val_T>> IdpAu, IdpAd, invIdpAu, invIdpAd;
		
		//IdpA=Id+A_as;
		matrix_sum(one, Id, one, Au_as, IdpAu);
		matrix_sum(one, Id, one, Ad_as, IdpAd);
		
		invert_matrix(IdpAu, invIdpAu);
		invert_matrix(IdpAd, invIdpAd);
		
		vector<vector<val_T>>  IdmAu, IdmAd, OMu, OMd;
		
		//IdpA=Id-A_as;
		matrix_sum(one, Id, m_one, Au_as, IdmAu);
		matrix_sum(one, Id, m_one, Ad_as, IdmAd);
		
		matrix_prod(IdmAu,invIdpAu,OMu);
		matrix_prod(IdmAd,invIdpAd,OMd);
		
	//	cout<<"OMu:\n";
	//	print_matrix(OMu);
	//	cout<<"OMd:\n";
	//	print_matrix(OMd);
		
	/*
		//check orthogonality
		vector<vector<val_T>> Id_pr, OMut, OMdt;
		transpose_matrix(OMu,OMut);
		transpose_matrix(OMd,OMdt);
		matrix_prod(OMu,OMut,Id_pr);
		cout<<"OMu*OMut\n";
		print_matrix(Id_pr);
		matrix_prod(OMd,OMdt,Id_pr);
		cout<<"OMd*OMdt\n";
		print_matrix(Id_pr);
	*/
		
		SD_T ref_SD_tmp=0;
		for (i=0; i<Nu; i++) setbit(ref_SD_tmp, i);
		for (i=0; i<Nd; i++) setbit(ref_SD_tmp, i+L);
		
		set_SO_basis(OMu, OMd, ref_SD_tmp);
	//	set_SO_basis(Id, Id, ref_SD_tmp);
		
		vector<vector<val_T>> abs_basis_u_prec=OMu, abs_basis_d_prec=OMd;
		double sgn;
		
		char JOBZ[1];
		JOBZ[0]='V';
		char UPLO[1];
		UPLO[0]='U';
		int NA[1];
		NA[0]=L;
		int LDA[1];
		LDA[0]=L;
		int LWORK[1];
		LWORK[0]=3*L-1;
		int INFO[1];
		
		double *W=new double[L];
		double *A=new double[L*L];
		double *WORK=new double[LWORK[0]];
		
		vector<vector<val_T>> rel_basis_u(L), abs_basis_u, abs_basis_u_tmp1, abs_basis_u_tmp2=Id, rel_basis_d(L), abs_basis_d, abs_basis_d_tmp1, abs_basis_d_tmp2=Id;
		
		for (i=0; i<L; i++)
		{
			rel_basis_u[i].resize(L);
			rel_basis_d[i].resize(L);
		}
		
		val_T r, dev_u, dev_d, dev_t_Fock_u, dev_t_Fock_d;
		vector<val_T> eigs_u(L), eigs_u_prec(L), eigs_d(L), eigs_d_prec(L);
		
		int N_UT=L*(L-1)/2;
		val_T r_max=1.0+Dr/2.0;
		
		vector<val_T> SO_v1(L), SO_v2(L);
		val_T SO_v_max;
		int i_max;
		
		for (r=r_init; r<r_max; r+=Dr)
		{
			eigs_u_prec.assign(L,0);
			eigs_d_prec.assign(L,0);
			k=0;
			dev_u=1;
			dev_d=1;
			while (k<Niter_max && (dev_u>tol || dev_d>tol || dev_t_Fock_u>tol || dev_t_Fock_d>tol))
			{
			//	cout<<"r="<<r<<"\t"<<"k="<<k<<"\t"<<"U_max: "<<r*H_T::U[0]<<endl;
				//diagonalize t_Fock_u
				for (i=0; i<L; i++)
					for (j=0; j<L; j++)
						A[i+j*L]=dval(t_Fock_u[i][j]);
				
				dsyev_(JOBZ, UPLO, NA, A, LDA, W, WORK, LWORK, INFO);
				
				if (!INFO[0])
				{
					dev_u=0;
				//	cout<<"eigs(t_Fock_u):\n";
					for (i=0; i<L; i++)
					{
						eigs_u[i]=W[i];
						dev_u+=(eigs_u[i]-eigs_u_prec[i])*(eigs_u[i]-eigs_u_prec[i]);
				//		cout<<W[i]<<endl;
						for (j=0; j<L; j++)
							rel_basis_u[i][j]=A[i+j*L];
					}
					dev_u=sqrt(dev_u);
					eigs_u_prec=eigs_u;
				}
				else
				{
				//	r=2.0;
					cerr<<"compute_Hartree_Fock_sol() error: no solution found\n";
					exit(EXIT_FAILURE);
				}
				
				//diagonalize t_Fock_d
				for (i=0; i<L; i++)
					for (j=0; j<L; j++)
						A[i+j*L]=dval(t_Fock_d[i][j]);
				
				dsyev_(JOBZ, UPLO, NA, A, LDA, W, WORK, LWORK, INFO);
				
				if (!INFO[0])
				{
					dev_d=0;
				//	cout<<"eigs(t_Fock_d):\n";
					for (i=0; i<L; i++)
					{
						eigs_d[i]=W[i];
						dev_d+=(eigs_d[i]-eigs_d_prec[i])*(eigs_d[i]-eigs_d_prec[i]);
				//		cout<<W[i]<<endl;
						for (j=0; j<L; j++)
							rel_basis_d[i][j]=A[i+j*L];
					}
					dev_d=sqrt(dev_d);
					eigs_d_prec=eigs_d;
				}
				else
				{
				//	r=2.0;
					cerr<<"compute_Hartree_Fock_sol() error: no solution found\n";
					exit(EXIT_FAILURE);
				}
				
				matrix_prod(SO_basis_u,rel_basis_u,abs_basis_u_tmp1);
				matrix_prod(SO_basis_d,rel_basis_d,abs_basis_d_tmp1);
				
				
				for (j=0; j<L; j++)
				{
					for (i=0; i<L; i++)
					{
						SO_v1[i]=abs_basis_u_prec[i][j];
						SO_v2[i]=abs_basis_u_tmp1[i][j];
					}
					SO_v_max=abs(SO_v1[0]);
					i_max=0;
					for (i=1; i<L; i++)
						if (abs(SO_v1[i])>SO_v_max)
						{
							i_max=i;
							SO_v_max=SO_v1[i];
						}
					sgn=1;
					if (real(SO_v1[i_max])*real(SO_v2[i_max])<0) sgn=-1;
					for (i=0; i<L; i++)
					{
						abs_basis_u_tmp2[i][j]=g*SO_v1[i]+sgn*(1.0-g)*SO_v2[i];
					}
					
					for (i=0; i<L; i++)
					{
						SO_v1[i]=abs_basis_d_prec[i][j];
						SO_v2[i]=abs_basis_d_tmp1[i][j];
					}
					SO_v_max=abs(SO_v1[0]);
					i_max=0;
					for (i=1; i<L; i++)
						if (abs(SO_v1[i])>SO_v_max)
						{
							i_max=i;
							SO_v_max=SO_v1[i];
						}
					sgn=1;
					if (real(SO_v1[i_max])*real(SO_v2[i_max])<0) sgn=-1;
					for (i=0; i<L; i++)
					{
						abs_basis_d_tmp2[i][j]=g*SO_v1[i]+sgn*(1.0-g)*SO_v2[i];
					}
				}
				 
				Gram_Schmidt_orthog(abs_basis_u_tmp2, abs_basis_u);
				Gram_Schmidt_orthog(abs_basis_d_tmp2, abs_basis_d);
				abs_basis_u_prec=abs_basis_u;
				abs_basis_d_prec=abs_basis_d;
				
			//	Gram_Schmidt_orthog(abs_basis_u_tmp1, abs_basis_u);
			//	Gram_Schmidt_orthog(abs_basis_d_tmp1, abs_basis_d);
				
				set_SO_basis(abs_basis_u, abs_basis_d, ref_SD_tmp, r);
				
				dev_t_Fock_u=0;
				dev_t_Fock_d=0;
				for (i=0; i<L-1; i++)
				{
					for (j=i+1; j<L; j++)
					{
						dev_t_Fock_u+=t_Fock_u[i][j]*t_Fock_u[i][j];
						dev_t_Fock_d+=t_Fock_d[i][j]*t_Fock_d[i][j];
					}
				}
				dev_t_Fock_u=sqrt(dev_t_Fock_u/N_UT);
				dev_t_Fock_d=sqrt(dev_t_Fock_d/N_UT);
				
				k++;
				
			}
			
			if (k==Niter_max)
			{
				cerr<<"compute_Hartree_Fock_sol() error: Hartree-Fock orbitals cannot converge\n";
				exit(EXIT_FAILURE);
			}
			
		}
		 
		if (W) delete [] W;
		if (A) delete [] A;
		if (WORK) delete [] WORK;
		
		cout<<"t_Fock_u:\n";
		print_matrix(t_Fock_u);
		cout<<"t_Fock_d:\n";
		print_matrix(t_Fock_d);
		
		cout<<"SO_basis_u:\n";
		print_matrix(abs_basis_u);
		cout<<"SO_basis_d:\n";
		print_matrix(abs_basis_d);
		
		vector<vector<val_T>> sum_SO_bases, diff_SO_bases;
		matrix_sum(one,abs_basis_u,one,abs_basis_d,sum_SO_bases);
		matrix_sum(one,abs_basis_u,m_one,abs_basis_d,diff_SO_bases);
		cout<<"sum_SO_bases:\n";
		print_matrix(sum_SO_bases);
		cout<<"diff_SO_bases:\n";
		print_matrix(diff_SO_bases);
		
	//	val_T E_HF=0;
		cout<<"eigs(t_Fock_u):\n";
		for (i=0; i<Nu; i++)
		{
			cout<<eigs_u[i]<<endl;
		//	E_HF+=eigs_u[i];
		}
		cout<<"eigs(t_Fock_d):\n";
		for (i=0; i<Nd; i++)
		{
			cout<<eigs_d[i]<<endl;
		//	E_HF+=eigs_d[i];
		}
		
	//	cout<<"E_HF: "<<E_HF<<endl;
		
	}
}

template<class ind_T, class val_T, class SD_T>
void SO_basis_Hamiltonian<ind_T, val_T, SD_T>::set_SO_basis(const vector<vector<val_T>> &SO_basis_u_par, const vector<vector<val_T>> &SO_basis_d_par, SD_T ref_par, val_T r)
{
	SO_basis_u=SO_basis_u_par;
	SO_basis_d=SO_basis_d_par;
	
	set_t_SO_basis();
	set_V_SO_basis();
	
	ref_SD=ref_par;
	set_ph_indices();
	define_SO_ph_map();
	set_t_Fock(r);
}

template<class ind_T, class val_T, class SD_T>
void SO_basis_Hamiltonian<ind_T, val_T, SD_T>::define_SO_ph_map()
{
	int i;
	SO_ph_map.resize(2*L);
	
	//	print_binary_string(ref_SD,2*L);
	//	cout<<endl;
	
//	cout<<"SO_ph_map:\n";
	int hu_ind=0, pu_ind=Nu, hd_ind=L, pd_ind=L+Nd;
	for (i=0; i<L; i++)
	{
		if (tstbit(ref_SD,i))
		{
			SO_ph_map[i]=hu_ind;
			hu_ind++;
		}
		else
		{
			SO_ph_map[i]=pu_ind;
			pu_ind++;
		}
//		cout<<setw(10)<<i<<SO_ph_map[i]<<endl;
	}
	for (i=L; i<2*L; i++)
	{
		if (tstbit(ref_SD,i))
		{
			SO_ph_map[i]=hd_ind;
			hd_ind++;
		}
		else
		{
			SO_ph_map[i]=pd_ind;
			pd_ind++;
		}
//		cout<<setw(10)<<i<<SO_ph_map[i]<<endl;
	}
	
}

template<class ind_T, class val_T, class SD_T>
void SO_basis_Hamiltonian<ind_T, val_T, SD_T>::set_ph_indices()
{
	unsigned i,j,k;
	
//	cout<<"set_ph_indices()\n";
	
	holes_ind.resize(N);
	particles_ind.resize(2*L-N);
	ph_indices.resize(2*L);
	
	j=0;
	k=0;
	for (i=0; i<2*L; i++)
	{
		if (tstbit(ref_SD,i))
		{
			holes_ind[j]=i;
			j++;
		}
		else
		{
			particles_ind[k]=i;
			k++;
		}
	}
	
	if (j!=N || k!=(2*L-N))
	{
		cerr<<"set_ph_indices() error: reference SD does not match the number of particles in the system\n";
	}
	
	j=0;
	k=Nu;
	for (i=0; i<L; i++)
	{
		if (tstbit(ref_SD,i))
		{
			ph_indices[j]=i;
			j++;
		}
		else
		{
			ph_indices[k]=i;
			k++;
		}
	}
	j=L;
	k=L+Nd;
	for (i=L; i<2*L; i++)
	{
		if (tstbit(ref_SD,i))
		{
			ph_indices[j]=i;
			j++;
		}
		else
		{
			ph_indices[k]=i;
			k++;
		}
	}
	
//	cout<<"ref_SD:\n";
//	print_binary_string(ref_SD,2*L);
//	cout<<endl;
	
/*	cout<<"ph_indices:\n";
	for (i=0; i<2*L; i++)
	{
		cout<<setw(10)<<i<<ph_indices[i]<<endl;
	}
*/
/*
	cout<<"holes_ind:\n";
	for (i=0; i<N; i++)
	{
		cout<<holes_ind[i]<<endl;
	}
	cout<<"particles_ind:\n";
	for (i=0; i<N; i++)
	{
		cout<<particles_ind[i]<<endl;
	}
*/
	
}

//compute the Fock matrix
template<class ind_T, class val_T, class SD_T>
void SO_basis_Hamiltonian<ind_T, val_T, SD_T>::set_t_Fock(val_T r)
{
	ind_T ind, i,j,k,l,m,n;
	val_T V_tmp;
	auto V_it=V_u_SO_basis.begin();
	
	t_Fock_u=t_SO_basis_u;
	t_Fock_d=t_SO_basis_d;
	
	
	for (i=0; i<L; i++)
	{
		for (j=0; j<L; j++)
		{
			//up part of ref_SD
			for (k=0; k<Nu; k++)
			{
				l=holes_ind[k];
				
				//t_Fock_u with up-up interaction
				ind=i+j*s[1]+l*s[2]+l*s[3];
				V_it=V_u_SO_basis.find(ind);
				if (V_it!=V_u_SO_basis.end())
					t_Fock_u[i][j]+=r*V_it->second;
				
				//t_Fock_d with down-up interaction
				ind=i+L+(j+L)*s[1]+l*s[2]+l*s[3];
				V_it=V_ud_SO_basis.find(ind);
				if (V_it!=V_ud_SO_basis.end())
					t_Fock_d[i][j]+=r*V_it->second;
				
			}
			//down part of ref_SD
			for (k=Nu; k<N; k++)
			{
				l=holes_ind[k];
				
				//t_Fock_u with up-down interaction
				ind=i+j*s[1]+l*s[2]+l*s[3];
				V_it=V_ud_SO_basis.find(ind);
				if (V_it!=V_ud_SO_basis.end())
					t_Fock_u[i][j]+=r*V_it->second;
				
				//t_Fock_d with down-down interaction
				ind=i+L+(j+L)*s[1]+l*s[2]+l*s[3];
				V_it=V_d_SO_basis.find(ind);
				if (V_it!=V_d_SO_basis.end())
					t_Fock_d[i][j]+=r*V_it->second;
			}
		}
	}
	
//	cout<<"t_Fock_u:\n";
//	print_matrix(t_Fock_u);
	
//	cout<<"t_Fock_d:\n";
//	print_matrix(t_Fock_d);

}

template<class ind_T, class val_T, class SD_T>
void SO_basis_Hamiltonian<ind_T, val_T, SD_T>::set_t_SO_basis()
{
	t_SO_basis_u.clear();
	t_SO_basis_d.clear();
	
	vector<vector<val_T>> M_u, M_d;
	matrix_element<ind_T, val_T> m_elem;
	
	M_u.resize(L);
	M_d.resize(L);
	
	int i,j,k;
	
	//compute M_u=t*SO_basis_u and M_d=t*SO_basis_d
	for (i=0; i<L; i++)
	{
		M_u[i].resize(L,0);
		M_d[i].resize(L,0);
		for (j=0; j<L; j++)
		{
			for (k=0; k<H_T::t[i].size(); k++)
			{
				M_u[i][j]+=H_T::t[i][k].value*SO_basis_u[H_T::t[i][k].index][j];
				M_d[i][j]+=H_T::t[i][k].value*SO_basis_d[H_T::t[i][k].index][j];
			}
		}
	}
	
	// compute t_SO_basis_u=transpose(SO_basis_u)*M_u and t_SO_basis_d=transpose(SO_basis_d)*M_d
	t_SO_basis_u.resize(L);
	t_SO_basis_d.resize(L);
	for (i=0; i<L; i++)
	{
		t_SO_basis_u[i].resize(L,0);
		t_SO_basis_d[i].resize(L,0);
		for (j=0; j<L; j++)
		{
			for (k=0; k<L; k++)
			{
				t_SO_basis_u[i][j]+=conj(SO_basis_u[k][i])*M_u[k][j];
				t_SO_basis_d[i][j]+=conj(SO_basis_d[k][i])*M_d[k][j];
			}
		}
	}
	
}

/*
template<class ind_T, class val_T, class SD_T>
void SO_basis_Hamiltonian<ind_T, val_T, SD_T>::transform_V(vector<vector<val_T>> basis_u, vector<vector<val_T>> basis_d)
{
	V_u_tr.clear();
	V_d_tr.clear();
	V_ud_tr.clear();
	
	ind_T ind, i,j,k,l,m,n,q,r;
	val_T Vtmp, val_V, Um,Uq,Ur,Un;
	array<ind_T,4> indices;
	
	for (r=0; r<L; r++)
		for (q=0; q<L; q++)
			for (n=0; n<L; n++)
				for (m=0; m<L; m++)
				{
					Vtmp=0;
					for (auto &elem: V_u_SO_basis)
					{
						ind=elem.first;
						val_V=elem.second;
						
						H_T::get_indices(ind, indices);
						i=indices[0];
						j=indices[1];
						k=indices[2];
						l=indices[3];
						
						if (i<L && j<L && k<L && l<L)
						{
							Um=conj(basis_u[i][m]);
							Un=basis_u[j][n];
							Uq=conj(basis_u[k][q]);
							Ur=basis_u[l][r];
							Vtmp+=val_V*Um*Un*Uq*Ur;
						}
					}
					if (Vtmp)
					{
						ind=m + n*s[1] + q*s[2] + r*s[3];
						V_u_tr.insert(pair<ind_T,val_T>(ind,Vtmp));
					}
				}
	
	for (r=L; r<2*L; r++)
		for (q=L; q<2*L; q++)
			for (n=L; n<2*L; n++)
				for (m=L; m<2*L; m++)
				{
					Vtmp=0;
					for (auto &elem: V_d_SO_basis)
					{
						ind=elem.first;
						val_V=elem.second;
						
						H_T::get_indices(ind, indices);
						i=indices[0];
						j=indices[1];
						k=indices[2];
						l=indices[3];
						
						if (i>=L && j>=L && k>=L && l>=L)
						{
							Um=conj(basis_d[i%L][m%L]);
							Un=basis_d[j%L][n%L];
							Uq=conj(basis_d[k%L][q%L]);
							Ur=basis_d[l%L][r%L];
							Vtmp+=val_V*Um*Un*Uq*Ur;
						}
					}
					if (Vtmp)
					{
						ind=m + n*s[1] + q*s[2] + r*s[3];
						V_d_tr.insert(pair<ind_T,val_T>(ind,Vtmp));
					}
				}
	
	for (r=L; r<2*L; r++)
		for (q=L; q<2*L; q++)
			for (n=0; n<L; n++)
				for (m=0; m<L; m++)
				{
					Vtmp=0;
					for (auto &elem: V_ud_SO_basis)
					{
						ind=elem.first;
						val_V=elem.second;
						
						H_T::get_indices(ind, indices);
						i=indices[0];
						j=indices[1];
						k=indices[2];
						l=indices[3];
						
						if (i<L && j<L && k>=L && l>=L)
						{
							Um=conj(basis_u[i][m]);
							Un=basis_u[j][n];
							Uq=conj(basis_d[k%L][q%L]);
							Ur=basis_d[l%L][r%L];
							Vtmp+=val_V*Um*Un*Uq*Ur;
						}
					}
					if (Vtmp)
					{
						H_T::insert_elem_V(Vtmp, m, n, q, r, V_ud_tr);
						H_T::insert_elem_V(-Vtmp, m, r, q, n, V_ud_tr);
						H_T::insert_elem_V(-Vtmp, q, n, m, r, V_ud_tr);
						H_T::insert_elem_V(Vtmp, q, r, m, n, V_ud_tr);
					}
				}
}
*/
/*
template<class ind_T, class val_T, class SD_T>
void SO_basis_Hamiltonian<ind_T, val_T, SD_T>::compare_V()
{
	ind_T ind;
	auto it=V_u_copy.begin();

	cout<<"ind\t\t"<<"V_u_copy\t\t"<<"V_u_tr"<<endl;
	for (auto &elem: V_u_copy)
	{
		ind=elem.first;
		it=V_u_tr.find(ind);
		if (it!=V_u_tr.end())
		{
			cout<<ind<<"\t\t"<<elem.second<<"\t\t"<<it->second<<endl;
		}
		else
		{
			cout<<ind<<endl;
		}
	}
	cout<<"ind\t\t"<<"V_d_copy\t\t"<<"V_d_tr"<<endl;
	for (auto &elem: V_d_copy)
	{
		ind=elem.first;
		it=V_d_tr.find(ind);
		if (it!=V_d_tr.end())
		{
			cout<<ind<<"\t\t"<<elem.second<<"\t\t"<<it->second<<endl;
		}
		else
		{
			cout<<ind<<endl;
		}
	}
	cout<<"ind\t\t"<<"V_ud_copy\t\t"<<"V_ud_tr"<<endl;
	for (auto &elem: V_ud_copy)
	{
		ind=elem.first;
		it=V_ud_tr.find(ind);
		if (it!=V_ud_tr.end())
		{
			cout<<ind<<"\t\t"<<elem.second<<"\t\t"<<it->second<<endl;
		}
		else
		{
			cout<<ind<<endl;
		}
	}
 */
 /*
	cout<<"ind\t\t"<<"V_u_copy\t\t"<<"V_u_tr"<<endl;
	for (auto &elem: V_u_tr)
	{
		ind=elem.first;
		it=V_u_copy.find(ind);
		if (it!=V_u_copy.end())
		{
			cout<<ind<<"\t\t"<<elem.second<<"\t\t"<<it->second<<endl;
		}
		else
		{
			cout<<ind<<"\t\t"<<elem.second<<endl;
		}
	}
	cout<<"ind\t\t"<<"V_d_copy\t\t"<<"V_d_tr"<<endl;
	for (auto &elem: V_d_tr)
	{
		ind=elem.first;
		it=V_d_copy.find(ind);
		if (it!=V_d_copy.end())
		{
			cout<<ind<<"\t\t"<<elem.second<<"\t\t"<<it->second<<endl;
		}
		else
		{
			cout<<ind<<"\t\t"<<elem.second<<endl;
		}
	}
	cout<<"ind\t\t"<<"V_ud_copy\t\t"<<"V_ud_tr"<<endl;
	for (auto &elem: V_ud_tr)
	{
		ind=elem.first;
		it=V_ud_copy.find(ind);
		if (it!=V_ud_copy.end())
		{
			cout<<ind<<"\t\t"<<elem.second<<"\t\t"<<it->second<<endl;
		}
		else
		{
			cout<<ind<<"\t\t"<<elem.second<<endl;
		}
	}
  */
//}
/*
template<class ind_T, class val_T, class SD_T>
void SO_basis_Hamiltonian<ind_T, val_T, SD_T>::copy_V()
{
	V_u_copy.clear();
	V_d_copy.clear();
	V_ud_copy.clear();
	
	for (auto &elem: V_u_SO_basis)
	{
		V_u_copy.insert(pair<ind_T,val_T>(elem.first,elem.second));
	}
	for (auto &elem: V_d_SO_basis)
	{
		V_d_copy.insert(pair<ind_T,val_T>(elem.first,elem.second));
	}
	for (auto &elem: V_ud_SO_basis)
	{
		V_ud_copy.insert(pair<ind_T,val_T>(elem.first,elem.second));
	}
}
*/

template<class ind_T, class val_T, class SD_T>
void SO_basis_Hamiltonian<ind_T, val_T, SD_T>::set_V_SO_basis()
{
	V_u_SO_basis.clear();
	V_d_SO_basis.clear();
	V_ud_SO_basis.clear();
	
	ind_T ind, i,j,k,l,m,n,q,r;
	val_T Vtmp, val_V, Um,Uq,Ur,Un;
	array<ind_T,4> indices;
	
	for (r=0; r<L; r++)
		for (q=0; q<L; q++)
			for (n=0; n<L; n++)
				for (m=0; m<L; m++)
				{
					Vtmp=0;
					for (auto &elem: H_T::V_u)
					{
						ind=elem.first;
						val_V=elem.second;
						
						H_T::get_indices(ind, indices);
						i=indices[0];
						j=indices[1];
						k=indices[2];
						l=indices[3];
					
						if (i<L && j<L && k<L && l<L)
						{
							Um=conj(SO_basis_u[i][m]);
							Un=SO_basis_u[j][n];
							Uq=conj(SO_basis_u[k][q]);
							Ur=SO_basis_u[l][r];
							Vtmp+=val_V*Um*Un*Uq*Ur;
						}
					}
					//if (fabs(Vtmp)>tol_V_min*H_T::V_max)
					if (Vtmp)
					{
						ind=m + n*s[1] + q*s[2] + r*s[3];
						V_u_SO_basis.insert(pair<ind_T,val_T>(ind,Vtmp));
					}
				}
	
	for (r=L; r<2*L; r++)
		for (q=L; q<2*L; q++)
			for (n=L; n<2*L; n++)
				for (m=L; m<2*L; m++)
				{
					Vtmp=0;
					for (auto &elem: H_T::V_d)
					{
						ind=elem.first;
						val_V=elem.second;
						
						H_T::get_indices(ind, indices);
						i=indices[0];
						j=indices[1];
						k=indices[2];
						l=indices[3];
						
						if (i>=L && j>=L && k>=L && l>=L)
						{
							Um=conj(SO_basis_d[i%L][m%L]);
							Un=SO_basis_d[j%L][n%L];
							Uq=conj(SO_basis_d[k%L][q%L]);
							Ur=SO_basis_d[l%L][r%L];
							Vtmp+=val_V*Um*Un*Uq*Ur;
						}
					}
				//	if (fabs(Vtmp)>tol_V_min*H_T::V_max)
					if (Vtmp)
					{
						ind=m + n*s[1] + q*s[2] + r*s[3];
						V_d_SO_basis.insert(pair<ind_T,val_T>(ind,Vtmp));
					}
				}
	
	for (r=L; r<2*L; r++)
		for (q=L; q<2*L; q++)
			for (n=0; n<L; n++)
				for (m=0; m<L; m++)
				{
					Vtmp=0;
					for (auto &elem: H_T::V_ud)
					{
						ind=elem.first;
						val_V=elem.second;
						
						H_T::get_indices(ind, indices);
						i=indices[0];
						j=indices[1];
						k=indices[2];
						l=indices[3];
						
						if (i<L && j<L && k>=L && l>=L)
						{
							Um=conj(SO_basis_u[i][m]);
							Un=SO_basis_u[j][n];
							Uq=conj(SO_basis_d[k%L][q%L]);
							Ur=SO_basis_d[l%L][r%L];
							Vtmp+=val_V*Um*Un*Uq*Ur;
						}
					}
				//	if (fabs(Vtmp)>tol_V_min*H_T::V_max)
					if (Vtmp)
					{
						H_T::insert_elem_V(Vtmp, m, n, q, r, V_ud_SO_basis);
						H_T::insert_elem_V(-Vtmp, m, r, q, n, V_ud_SO_basis);
						H_T::insert_elem_V(-Vtmp, q, n, m, r, V_ud_SO_basis);
						H_T::insert_elem_V(Vtmp, q, r, m, n, V_ud_SO_basis);
				//		ind=m + n*s[1] + q*s[2] + r*s[3];
				//		V_ud_SO_basis.insert(pair<ind_T,val_T>(ind,Vtmp));
						
					}
				}
	
}

#endif /* SO_basis_Hamiltonian_h */
