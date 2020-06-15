//
//  full_CI_main.cpp
//
//  Created by Dom2 on 20-03-23.
//

#include "Hamiltonian_diag.h"
#include "misc.h"
//#include "SO_basis_Hamiltonian.h"

using namespace std;

int main(int arg_N, char *args[])
{
	mpf_set_default_prec(256);
	
	int i,j;
	
	int Nu, Nd;
	int L, L_max=20;
	using ind_T=int;
	using val_T=mpf_class;
	using SD_T=unsigned;
	using cfs_T=mpf_class;
	
	coord_vec_T lattice(L_max);
	
	ifstream params_file("../../full_CI.dat",ios::in);
	
	string str, displ_figs, iter_ref_SD, random_init_st;
	vector<val_T> U_par(3);
	vector<val_T> t_par(2);
	vector<unsigned> u;
	vector<unsigned> d;
	
	int Niter_Lanczos;
	int Niter_max_HF;
	val_T r_init_HF;
	val_T Dr_HF;
	val_T tol_HF;
	double f_A_HF;
	double g_HF;
	
	if (params_file)
	{
		params_file>>str>>U_par[0]>>U_par[1]>>U_par[2];
		params_file>>str>>t_par[0]>>t_par[1];
		params_file>>str>>Nu>>str>>Nd;
		cout<<"U: "<<U_par[0]<<"\t"<<U_par[1]<<"\t"<<U_par[2]<<endl;
		cout<<"t: "<<t_par[0]<<"\t"<<t_par[1]<<endl;
		cout<<"Nu: "<<Nu<<endl;
		cout<<"Nd: "<<Nd<<endl;
		u.resize(Nu);
		d.resize(Nd);
		getline(params_file,str);
		getline(params_file,str);
		L=0;
		i=0;
		while (L<=L_max && str.substr(0,6)!="pos_u:")
		{
			getline(params_file,str);
			for (j=0; j<str.size(); j++)
			{
				if (str[j]=='*')
				{
					lattice[L][0]=j;
					lattice[L][1]=i;
					L++;
			//		cout<<"x, y: "<<j<<"\t"<<i<<endl;
				}
			}
			i++;
		}
		cout<<"L: "<<L<<endl;
		if (str.substr(0,6)=="pos_u:")
		{
			params_file.seekg(-(str.size()+1),ios::cur);
			params_file>>str;
			cout<<str;
			for (i=0; i<Nu; i++)
			{
				params_file>>u[i];
				cout<<' '<<u[i];
			}
			cout<<endl;
			params_file>>str;
			cout<<str;
			for (i=0; i<Nd; i++)
			{
				params_file>>d[i];
				cout<<' '<<d[i];
			}
			cout<<endl;
		}
		params_file>>str>>Niter_Lanczos;
		cout<<str<<'\t'<<Niter_Lanczos<<endl;
		params_file>>str>>Niter_max_HF;
		cout<<str<<'\t'<<Niter_max_HF<<endl;
		params_file>>str>>r_init_HF;
		cout<<str<<'\t'<<r_init_HF<<endl;
		params_file>>str>>Dr_HF;
		cout<<str<<'\t'<<Dr_HF<<endl;
		params_file>>str>>tol_HF;
		cout<<str<<'\t'<<tol_HF<<endl;
		params_file>>str>>f_A_HF;
		cout<<str<<'\t'<<f_A_HF<<endl;
		params_file>>str>>g_HF;
		cout<<str<<'\t'<<g_HF<<endl;
		params_file>>str>>displ_figs;
		cout<<str<<'\t'<<displ_figs<<endl;
		params_file>>str>>random_init_st;
		cout<<str<<'\t'<<random_init_st<<endl;
		params_file>>str>>iter_ref_SD;
		cout<<str<<'\t'<<iter_ref_SD<<endl;
	}
	else
	{
		cerr<<"unable to open file\n";
		exit(EXIT_FAILURE);
	}
	lattice.resize(L);
	array<unsigned,2> N_u_d;
	N_u_d[0]=Nu;
	N_u_d[1]=Nd;
	int N=Nu+Nd;
	
	bool display_figures=false;
	if (displ_figs=="yes") display_figures=true;
	bool iterate_ref_SD=false;
	if (iter_ref_SD=="yes") iterate_ref_SD=true;
	bool rnd_init_st=false;
	if (random_init_st=="yes") rnd_init_st=true;
	
	
	//6 sites: U_c=2.69 => U_init_HF
	//8 sites U_c=1.9
	cout<<"U_init_HF: "<<r_init_HF*U_par[0]<<endl;
	cout<<"DU_HF: "<<Dr_HF*U_par[0]<<endl;
	

	vector<val_T> onsite_E;
	
	double tol=1e-8;
	Hamiltonian_diag<ind_T, val_T, SD_T, cfs_T> H_d(N_u_d, tol, lattice, U_par, t_par, onsite_E);
	
	vector<vector<val_T>> Id(L);
	
	for (i=0; i<L; i++)
	{
		Id[i].resize(L,0);
		Id[i][i]=1.0;
	}
	
	unsigned ref_SD=0;
	for (i=0; i<Nu; i++) setbit(ref_SD, u[i]);
	for (i=0; i<Nd; i++) setbit(ref_SD, d[i]+L);
	
	//compute the ground state using the local spin-orbital basis
	H_d.set_SO_basis(Id,Id,ref_SD);
	H_d.create_CI_ph_basis();
	H_d.compute_ground_state(Niter_Lanczos, iterate_ref_SD, display_figures, rnd_init_st);

	//compute the ground state using the natural spin-orbital computed with the preceding ground state
	H_d.compute_natural_spin_orbitals();
	H_d.create_CI_ph_basis();
	H_d.compute_ground_state(Niter_Lanczos, iterate_ref_SD, display_figures);

	//compute the ground state using the unrestricted Hartree-Fock orbitals
	H_d.compute_Hartree_Fock_sol(Niter_max_HF,r_init_HF,Dr_HF,tol_HF,f_A_HF,g_HF);
	H_d.create_CI_ph_basis();
	H_d.compute_ground_state(Niter_Lanczos, iterate_ref_SD, display_figures);

	return 0;
}
