//
//  Hamiltonian.h
//  full_CI
//
//  Copyright (C) Dominic Bergeron 2020
//

#ifndef Hamiltonian_h
#define Hamiltonian_h

#include "sparse_matrix.h"
#include <array>

typedef array<int, 2> site_T;
typedef vector<site_T> coord_vec_T;

template<class ind_T, class val_T>
class Hamiltonian
{
public:
	Hamiltonian(const coord_vec_T &lattice_par, const vector<val_T> &U_par, const vector<val_T> &t_par, const vector<val_T> &onsite_E_par);
	Hamiltonian(const Hamiltonian<ind_T, val_T> &H_par);
//	Hamiltonian & operator = (const Hamiltonian &H_par) {return *this;}
//	Hamiltonian(fstream &input_file);
	
	coord_vec_T lattice;
	ind_T L;
	array<ind_T,4> s;
	vector<val_T> U, t_rel, onsite_E;
	val_T t_max, V_max;
	int xmin, xmax, ymin, ymax, zmin, zmax;
	bool H_real;

	vector<vector<val_T>> t_m;
	sparse_matrix_vec<ind_T, val_T> t;
	//Vijkl is stored according to delta_{s_i s_j}delta_{s_k s_l}Vijkl c_i^+c_k^+c_lc_j, where s_i is the spin of spin-orbital i, and is antisymmetrized: Vijkl=-Vilkj
	map<ind_T, val_T> V_u, V_d, V_ud;
	
//	mat W_u,X_u, W_d,X_d, W_ud,X_ud;
//	vec S_u,S_d,S_ud;
	
	void get_indices(ind_T ind, array<ind_T,4> &indices);
	void insert_elem_V(val_T V_val, ind_T k1, ind_T k2, ind_T k3, ind_T k4, map<ind_T, val_T> &V_map);
	
//	void decompose_V();
	
//private:
	
	void set_t_Hubbard_2D();
	void set_V_Hubbard_2D();
	
	void insert_elems_Hubbard(ind_T i, ind_T j, val_T V_val, map<ind_T, val_T> &V_map);
	void insert_elems_ext_Hubbard(ind_T i, ind_T j, val_T V_val);
};

template<class ind_T, class val_T>
Hamiltonian<ind_T, val_T>::Hamiltonian(const Hamiltonian<ind_T, val_T> &H_par)
{
	lattice=H_par.lattice;
	L=lattice.size();
	U=H_par.U;
	t_rel=H_par.t_rel;
	onsite_E=H_par.onsite_E;
	t_max=H_par.t_max;
	V_max=H_par.V_max;
	xmin=H_par.xmin;
	xmax=H_par.xmax;
	ymin=H_par.ymin;
	ymax=H_par.ymax;
	zmin=H_par.zmin;
	zmax=H_par.zmax;
	
	t=H_par.t;
	V_u=H_par.V_u;
	V_d=H_par.V_d;
	V_ud=H_par.V_ud;
	
	s=H_par.s;
	
	H_real=H_par.H_real;
}


template<class ind_T, class val_T>
Hamiltonian<ind_T, val_T>::Hamiltonian(const coord_vec_T &lattice_par, const vector<val_T> &U_par, const vector<val_T> &t_par, const vector<val_T> &onsite_E_par)
{
	int i=0;
	
	U=U_par;
	t_rel=t_par;
	if (!U_par.size() || !t_par.size())
	{
		cerr<<"Hamiltonian(): empty U or t parameter\n";
		return;
	}
	onsite_E=onsite_E_par;
	lattice=lattice_par;
	L=lattice.size();
	if (L<=0)
	{
		cerr<<"Hamiltonian(): empty lattice\n";
		return;
	}
	s[0]=1;
	s[1]=2*L;
	s[2]=4*L*L;
	s[3]=8*L*L*L;
	
	t_max=fabs(t_rel[0]);
	for (i=1; i<t_rel.size(); i++)
		if (fabs(t_rel[i])>t_max) t_max=fabs(t_rel[i]);
	
	V_max=fabs(U[0]);
	for (i=1; i<U.size(); i++)
		if (fabs(U[i])>V_max) V_max=fabs(U[i]);
	
	xmin=lattice[0][0];
	ymin=lattice[0][1];
	xmax=lattice[L-1][0];
	ymax=lattice[L-1][1];

	cout<<setiosflags(ios::left);
	cout<<"lattice:\n";
	for (i=0; i<L; i++)
	{
		cout<<setw(10)<<lattice[i][0]<<lattice[i][1]<<endl;
		if (lattice[i][0]<xmin) xmin=lattice[i][0];
		if (lattice[i][1]<ymin) ymin=lattice[i][1];
		if (lattice[i][0]>xmax) xmax=lattice[i][0];
		if (lattice[i][1]>ymax) ymax=lattice[i][1];
	}
	
	cout<<"xmin, xmax: "<<xmin<<", "<<xmax<<endl;
	cout<<"ymin, ymax: "<<ymin<<", "<<ymax<<endl;

	set_t_Hubbard_2D();
	set_V_Hubbard_2D();
	
	H_real=true;

/*
	cout<<"t:\n";
	for (i=0; i<L; i++)
	{
		cout<<setw(10)<<i;
		for (auto const &item : t[i])
		{
			cout<<'('<<item.index<<','<<item.value<<setw(10)<<')';
		}
		cout<<endl;
	}
	cout<<endl;

	ind_T j,k,l,m, ind;
	ind_T s1=2*L, s2=4*L*L, s3=8*L*L*L;
	cout<<"V_u:\n";
	for (auto const &elem: V_u)
	{
		ind=elem.first;
		l=ind/s3;
		m=ind % s3;
		k=m/s2;
		m=m % s2;
		j=m/s1;
		i=m % s1;
		cout<<i<<"("<<i%L<<setw(10)<<")"<<j<<"("<<j%L<<setw(10)<<")"<<k<<"("<<k%L<<setw(10)<<")"<<l<<"("<<l%L<<setw(10)<<")"<<elem.second<<endl;
	}
	cout<<endl;
	cout<<"V_d:\n";
	for (auto const &elem: V_d)
	{
		ind=elem.first;
		l=ind/s3;
		m=ind % s3;
		k=m/s2;
		m=m % s2;
		j=m/s1;
		i=m % s1;
		cout<<i<<"("<<i%L<<setw(10)<<")"<<j<<"("<<j%L<<setw(10)<<")"<<k<<"("<<k%L<<setw(10)<<")"<<l<<"("<<l%L<<setw(10)<<")"<<elem.second<<endl;
	}
	cout<<endl;
	cout<<"V_ud:\n";
	for (auto const &elem: V_ud)
	{
		ind=elem.first;
		l=ind/s3;
		m=ind % s3;
		k=m/s2;
		m=m % s2;
		j=m/s1;
		i=m % s1;
		cout<<i<<"("<<i%L<<setw(10)<<")"<<j<<"("<<j%L<<setw(10)<<")"<<k<<"("<<k%L<<setw(10)<<")"<<l<<"("<<l%L<<setw(10)<<")"<<elem.second<<endl;
	}
	cout<<endl;
*/
}

template<class ind_T, class val_T>
void Hamiltonian<ind_T, val_T>::get_indices(ind_T ind, array<ind_T,4> &indices)
{
	int i;
	ind_T n=ind;
	for (i=s.size()-1; i>=0; i--)
	{
		indices[i]=n/s[i];
		n=n % s[i];
	}
}

template<class ind_T, class val_T>
void Hamiltonian<ind_T, val_T>::insert_elem_V(val_T V_val, ind_T k1, ind_T k2, ind_T k3, ind_T k4, map<ind_T, val_T> &V_map)
{
	ind_T ind=k1 + k2*s[1] + k3*s[2] + k4*s[3];
	V_map.insert(pair<ind_T,val_T>(ind,V_val));
}


//for a given pair of spin-orbital indices (j,k), call insert_elem_V() to insert all the possible finite elements in the antisymmetrized interaction V_map, which is one of V_u, V_d, V_ud
template<class ind_T, class val_T>
void Hamiltonian<ind_T, val_T>::insert_elems_Hubbard(ind_T j, ind_T k, val_T V_val, map<ind_T, val_T> &V_map)
{
	insert_elem_V(V_val, j,j,k,k, V_map);
	insert_elem_V(-V_val, j,k,k,j, V_map);
	insert_elem_V(-V_val, k,j,j,k, V_map);
	insert_elem_V(V_val, k,k,j,j, V_map);
}

//For a given pair of sites (i,j) connected by the interaction, generate all pairs of spin-orbitals (i_u, j_u), (i_u, j_d), (i_d, j_u), (i_d, j_d) and call insert_elems_Hubbard() with the corresponding map among V_u, V_d, V_ud.
template<class ind_T, class val_T>
void Hamiltonian<ind_T, val_T>::insert_elems_ext_Hubbard(ind_T i, ind_T j, val_T V_val)
{
	insert_elems_Hubbard(i, j, V_val, V_u);
	insert_elems_Hubbard(i, j+L, V_val, V_ud);
	insert_elems_Hubbard(i+L, j, V_val, V_ud);
	insert_elems_Hubbard(i+L, j+L, V_val, V_d);
}

//fill the antisymmetrized interaction tensor
template<class ind_T, class val_T>
void Hamiltonian<ind_T, val_T>::set_V_Hubbard_2D()
{
	ind_T i,m;
	
	int Nnn=U.size();
	
	int x, y;
	for (i=0; i<L; i++)
	{
		insert_elems_Hubbard(i, i+L, U[0], V_ud);
		
		if (Nnn>1 && U[1])
		{
			//Find all the distinct pairs of nearest neighbor sites (i,m) by looking for each nearest neighbor to site i at an index m>i in the lattice. All the nearest neighbors of i are searched, so the spatial distribution of site indices is not important.
			
			x=lattice[i][0];
			y=lattice[i][1]-1;
			if (y>=ymin)
			{
				m=i+1;
				while (m<L && (lattice[m][0]!=x || lattice[m][1]!=y)) m++;
				if (m<L) insert_elems_ext_Hubbard(i, m, U[1]);
			}
			
			x=lattice[i][0]-1;
			y=lattice[i][1];
			if (x>=xmin)
			{
				m=i+1;
				while (m<L && (lattice[m][0]!=x || lattice[m][1]!=y)) m++;
				if (m<L) insert_elems_ext_Hubbard(i, m, U[1]);
			}
			
			x=lattice[i][0]+1;
			y=lattice[i][1];
			if (x<=xmax)
			{
				m=i+1;
				while (m<L && (lattice[m][0]!=x || lattice[m][1]!=y)) m++;
				if (m<L) insert_elems_ext_Hubbard(i, m, U[1]);
			}
			
			x=lattice[i][0];
			y=lattice[i][1]+1;
			if (y<=ymax)
			{
				m=i+1;
				while (m<L && (lattice[m][0]!=x || lattice[m][1]!=y)) m++;
				if (m<L) insert_elems_ext_Hubbard(i, m, U[1]);
			}
			
			if (Nnn>2 && U[2])
			{
				//Find all the distinct pairs of next nearest neighbor sites (i,m) by looking for each next nearest neighbor to site i at an index m>i in the lattice. All the next nearest neighbors of i are searched, so the spatial distribution of site indices is not important.
				
				x=lattice[i][0]-1;
				y=lattice[i][1]-1;
				if (y>=ymin && x>=xmin)
				{
					m=i+1;
					while (m<L && (lattice[m][0]!=x || lattice[m][1]!=y)) m++;
					if (m<L) insert_elems_ext_Hubbard(i, m, U[2]);
				}
				
				x=lattice[i][0]+1;
				y=lattice[i][1]-1;
				if (y>=ymin && x<=xmax)
				{
					m=i+1;
					while (m<L && (lattice[m][0]!=x || lattice[m][1]!=y)) m++;
					if (m<L) insert_elems_ext_Hubbard(i, m, U[2]);
				}
				
				x=lattice[i][0]-1;
				y=lattice[i][1]+1;
				if (y<=ymax && x>=xmin)
				{
					m=i+1;
					while (m<L && (lattice[m][0]!=x || lattice[m][1]!=y)) m++;
					if (m<L) insert_elems_ext_Hubbard(i, m, U[2]);
				}
				
				x=lattice[i][0]+1;
				y=lattice[i][1]+1;
				if (y<=ymax && x<=xmax)
				{
					m=i+1;
					while (m<L && (lattice[m][0]!=x || lattice[m][1]!=y)) m++;
					if (m<L) insert_elems_ext_Hubbard(i, m, U[2]);
				}
			}
		}
	}
}

template<class ind_T, class val_T>
void Hamiltonian<ind_T, val_T>::set_t_Hubbard_2D()
{
	int i,j;
	matrix_element<ind_T, val_T> elem;
	
	t.resize(L);
	if (onsite_E.size()>=L)
	{
		for (i=0; i<L; i++)
		{
			elem.index=i;
			elem.value=onsite_E[i];
			t[i].push_back(elem);
		}
	}
	
	int x,y;
	for (i=0; i<L-1; i++)
	{
		x=lattice[i][0];
		y=lattice[i][1]-1;
		if (y>=ymin)
		{
			j=i+1;
			while (j<L && (lattice[j][0]!=x || lattice[j][1]!=y)) j++;
			if (j<L)
			{
				elem.index=j;
				elem.value=t_rel[0];
				t[i].push_back(elem);
				elem.index=i;
				elem.value=conj(t_rel[0]);
				t[j].push_back(elem);
			}
		}
		
		x=lattice[i][0]-1;
		y=lattice[i][1];
		if (x>=xmin)
		{
			j=i+1;
			while (j<L && (lattice[j][0]!=x || lattice[j][1]!=y)) j++;
			if (j<L)
			{
				elem.index=j;
				elem.value=t_rel[0];
				t[i].push_back(elem);
				elem.index=i;
				elem.value=conj(t_rel[0]);
				t[j].push_back(elem);
			}
		}
		
		x=lattice[i][0]+1;
		y=lattice[i][1];
		if (x<=xmax)
		{
			j=i+1;
			while (j<L && (lattice[j][0]!=x || lattice[j][1]!=y)) j++;
			if (j<L)
			{
				elem.index=j;
				elem.value=t_rel[0];
				t[i].push_back(elem);
				elem.index=i;
				elem.value=conj(t_rel[0]);
				t[j].push_back(elem);
			}
		}
		
		x=lattice[i][0];
		y=lattice[i][1]+1;
		if (y<=ymax)
		{
			j=i+1;
			while (j<L && (lattice[j][0]!=x || lattice[j][1]!=y)) j++;
			if (j<L)
			{
				elem.index=j;
				elem.value=t_rel[0];
				t[i].push_back(elem);
				elem.index=i;
				elem.value=conj(t_rel[0]);
				t[j].push_back(elem);
			}
		}
		
		if (t_rel.size()>1 && t_rel[1])
		{
			x=lattice[i][0]-1;
			y=lattice[i][1]-1;
			if (y>=ymin && x>=xmin)
			{
				j=i+1;
				while (j<L && (lattice[j][0]!=x || lattice[j][1]!=y)) j++;
				if (j<L)
				{
					elem.index=j;
					elem.value=t_rel[1];
					t[i].push_back(elem);
					elem.index=i;
					elem.value=conj(t_rel[1]);
					t[j].push_back(elem);
				}
			}
			
			x=lattice[i][0]+1;
			y=lattice[i][1]-1;
			if (y>=ymin && x<=xmax)
			{
				j=i+1;
				while (j<L && (lattice[j][0]!=x || lattice[j][1]!=y)) j++;
				if (j<L)
				{
					elem.index=j;
					elem.value=t_rel[1];
					t[i].push_back(elem);
					elem.index=i;
					elem.value=conj(t_rel[1]);
					t[j].push_back(elem);
				}
			}
			
			x=lattice[i][0]-1;
			y=lattice[i][1]+1;
			if (y<=ymax && x>=xmin)
			{
				j=i+1;
				while (j<L && (lattice[j][0]!=x || lattice[j][1]!=y)) j++;
				if (j<L)
				{
					elem.index=j;
					elem.value=t_rel[1];
					t[i].push_back(elem);
					elem.index=i;
					elem.value=conj(t_rel[1]);
					t[j].push_back(elem);
				}
			}
			
			x=lattice[i][0]+1;
			y=lattice[i][1]+1;
			if (y<=ymax && x<=xmax)
			{
				j=i+1;
				while (j<L && (lattice[j][0]!=x || lattice[j][1]!=y)) j++;
				if (j<L)
				{
					elem.index=j;
					elem.value=t_rel[1];
					t[i].push_back(elem);
					elem.index=i;
					elem.value=conj(t_rel[1]);
					t[j].push_back(elem);
				}
			}
		}
	}
/*
	cout<<setiosflags(ios::left);
	for (i=0; i<L; i++)
	{
		for (j=0; j<t[i].size(); j++)
		{
			cout<<i<<'\t'<<t[i][j].index<<'\t'<<t[i][j].value<<endl;
		}
	}
*/
}

/*
template<class ind_T, class val_T>
void Hamiltonian<ind_T, val_T>::set_t_Hubbard_2D()
{
	int i,j;
	matrix_element<ind_T, val_T> elem;

	cout<<setiosflags(ios::left);
	
	t.resize(L);
	if (onsite_E.size()>=L)
	{
		for (i=0; i<L; i++)
		{
			elem.index=i;
			elem.value=onsite_E[i];
			t[i].push_back(elem);
		}
	}
	
	int x,y;
	for (i=0; i<L; i++)
	{
		x=lattice[i][0];
		y=lattice[i][1]-1;
		if (y>=ymin)
		{
			j=i-1;
			while (j>=0 && (lattice[j][0]!=x || lattice[j][1]!=y)) j--;
			if (j<0)
			{
				j=i+1;
				while (j<L && (lattice[j][0]!=x || lattice[j][1]!=y)) j++;
			}
			if (j>=0 && j<L)
			{
				elem.index=j;
				elem.value=t_rel[0];
				t[i].push_back(elem);
			}
		}
		
		x=lattice[i][0]-1;
		y=lattice[i][1];
		if (x>=xmin)
		{
			j=i-1;
			while (j>=0 && (lattice[j][0]!=x || lattice[j][1]!=y)) j--;
			if (j<0)
			{
				j=i+1;
				while (j<L && (lattice[j][0]!=x || lattice[j][1]!=y)) j++;
			}
			if (j>=0 && j<L)
			{
				elem.index=j;
				elem.value=t_rel[0];
				t[i].push_back(elem);
			}
		}
		
		x=lattice[i][0]+1;
		y=lattice[i][1];
		if (x<=xmax)
		{
			j=i+1;
			while (j<L && (lattice[j][0]!=x || lattice[j][1]!=y)) j++;
			if (j==L)
			{
				j=i-1;
				while (j>=0 && (lattice[j][0]!=x || lattice[j][1]!=y)) j--;
				
			}
			if (j>=0 && j<L)
			{
				elem.index=j;
				elem.value=t_rel[0];
				t[i].push_back(elem);
			}
		}
		
		x=lattice[i][0];
		y=lattice[i][1]+1;
		if (y<=ymax)
		{
			j=i+1;
			while (j<L && (lattice[j][0]!=x || lattice[j][1]!=y)) j++;
			if (j==L)
			{
				j=i-1;
				while (j>=0 && (lattice[j][0]!=x || lattice[j][1]!=y)) j--;
				
			}
			if (j>=0 && j<L)
			{
				elem.index=j;
				elem.value=t_rel[0];
				t[i].push_back(elem);
			}
		}
		
		if (t_rel.size()>1)
		{
			x=lattice[i][0]-1;
			y=lattice[i][1]-1;
			if (y>=ymin && x>=xmin)
			{
				j=i-1;
				while (j>=0 && (lattice[j][0]!=x || lattice[j][1]!=y)) j--;
				if (j<0)
				{
					j=i+1;
					while (j<L && (lattice[j][0]!=x || lattice[j][1]!=y)) j++;
				}
				if (j>=0 && j<L)
				{
					elem.index=j;
					elem.value=t_rel[1];
					t[i].push_back(elem);
				}
			}

			x=lattice[i][0]+1;
			y=lattice[i][1]-1;
			if (y>=ymin && x<=xmax)
			{
				j=i-1;
				while (j>=0 && (lattice[j][0]!=x || lattice[j][1]!=y)) j--;
				if (j<0)
				{
					j=i+1;
					while (j<L && (lattice[j][0]!=x || lattice[j][1]!=y)) j++;
				}
				if (j>=0 && j<L)
				{
					elem.index=j;
					elem.value=t_rel[1];
					t[i].push_back(elem);
				}
			}

			x=lattice[i][0]-1;
			y=lattice[i][1]+1;
			if (y<=ymax && x>=xmin)
			{
				j=i+1;
				while (j<L && (lattice[j][0]!=x || lattice[j][1]!=y)) j++;
				if (j==L)
				{
					j=i-1;
					while (j>=0 && (lattice[j][0]!=x || lattice[j][1]!=y)) j--;
				}
				if (j>=0 && j<L)
				{
					elem.index=j;
					elem.value=t_rel[1];
					t[i].push_back(elem);
				}
			}

			x=lattice[i][0]+1;
			y=lattice[i][1]+1;
			if (y<=ymax && x<=xmax)
			{
				j=i+1;
				while (j<L && (lattice[j][0]!=x || lattice[j][1]!=y)) j++;
				if (j==L)
				{
					j=i-1;
					while (j>=0 && (lattice[j][0]!=x || lattice[j][1]!=y)) j--;
				}
				if (j>=0 && j<L)
				{
					elem.index=j;
					elem.value=t_rel[1];
					t[i].push_back(elem);
				}
			}	
		}	
	}
	
//	for (i=0; i<L; i++)
//	{
//		for (j=0; j<t[i].size(); j++)
//		{
//			cout<<i<<'\t'<<t[i][j].index<<'\t'<<t[i][j].value<<endl;
//		}
//	}
	
}
*/

/*
 //compute the SVD of V_u, V_d, V_ud
 template<class ind_T, class val_T>
 void Hamiltonian<ind_T, val_T>::decompose_V()
 {
 ind_T ind, i,j,k,l,m,n,p,q,r;
 mat Vu_tmp(L*L,L*L, fill::zeros), Vd_tmp(L*L,L*L, fill::zeros), Vud_tmp(L*L,L*L, fill::zeros);
 
 cout<<"decomposing V\n";
 array<ind_T, 4> indices;
 for (auto &elem: V_u)
 {
 ind=elem.first;
 get_indices(ind, indices);
 i=indices[0];
 j=indices[1];
 k=indices[2];
 l=indices[3];
 
 m=i+j*L;
 n=k+l*L;
 Vu_tmp(m,n)=elem.second;
 
 //		p=i+l*s[1]+k*s[2]+j*s[3];
 //		q=k+j*s[1]+i*s[2]+l*s[3];
 //		r=k+l*s[1]+i*s[2]+j*s[3];
 //		if (V[p]!=-elem.second) cout<<"AS problem: "<<setw(10)<<i<<setw(10)<<j<<setw(10)<<k<<l<<endl;
 //		if (V[q]!=-elem.second) cout<<"AS problem: "<<setw(10)<<i<<setw(10)<<j<<setw(10)<<k<<l<<endl;
 //		if (V[r]!=elem.second) cout<<"sym. problem: "<<setw(10)<<i<<setw(10)<<j<<setw(10)<<k<<l<<endl;
 }
 
 for (auto &elem: V_d)
 {
 ind=elem.first;
 get_indices(ind, indices);
 i=indices[0];
 j=indices[1];
 k=indices[2];
 l=indices[3];
 
 i=i%L;
 j=j%L;
 k=k%L;
 l=l%L;
 m=i+j*L;
 n=k+l*L;
 Vd_tmp(m,n)=elem.second;
 
 //		p=i+l*s[1]+k*s[2]+j*s[3];
 //		q=k+j*s[1]+i*s[2]+l*s[3];
 //		r=k+l*s[1]+i*s[2]+j*s[3];
 //		if (V[p]!=-elem.second) cout<<"AS problem: "<<setw(10)<<i<<setw(10)<<j<<setw(10)<<k<<l<<endl;
 //		if (V[q]!=-elem.second) cout<<"AS problem: "<<setw(10)<<i<<setw(10)<<j<<setw(10)<<k<<l<<endl;
 //		if (V[r]!=elem.second) cout<<"sym. problem: "<<setw(10)<<i<<setw(10)<<j<<setw(10)<<k<<l<<endl;
 
 }
 
 for (auto &elem: V_ud)
 {
 ind=elem.first;
 get_indices(ind, indices);
 i=indices[0];
 j=indices[1];
 k=indices[2];
 l=indices[3];
 
 if (i<L && j<L && k>=L && l>=L)
 {
 k=k%L;
 l=l%L;
 m=i+j*L;
 n=k+l*L;
 Vud_tmp(m,n)=elem.second;
 }
 
 //		p=i+l*s[1]+k*s[2]+j*s[3];
 //		q=k+j*s[1]+i*s[2]+l*s[3];
 //		r=k+l*s[1]+i*s[2]+j*s[3];
 //		if (V[p]!=-elem.second) cout<<"AS problem: "<<setw(10)<<i<<setw(10)<<j<<setw(10)<<k<<l<<endl;
 //		if (V[q]!=-elem.second) cout<<"AS problem: "<<setw(10)<<i<<setw(10)<<j<<setw(10)<<k<<l<<endl;
 //		if (V[r]!=elem.second) cout<<"sym. problem: "<<setw(10)<<i<<setw(10)<<j<<setw(10)<<k<<l<<endl;
 }
 */
/*
 mat Vt=Vu_tmp.t();
 mat dV=Vu_tmp-Vt;
 
 cout<<"max(Vd_tmp-Vt):\n";
 cout<<max(max(abs(dV)))<<endl;
 
 cout<<"Vu_tmp:\n";
 cout<<Vu_tmp.cols(0,12);
 */
//	svd(W_u,S_u,X_u,Vu_tmp);
//	svd(W_d,S_d,X_d,Vd_tmp);
//	svd(W_ud,S_ud,X_ud,Vud_tmp);

/*
 cout<<setiosflags(ios::left);
 cout<<"singular values:\n";
 cout<<"S_u\n";
 for (i=0; i<L*L; i++)
 {
 cout<<S_u(i)<<endl;
 }
 cout<<endl;
 
 cout<<"S_d\n";
 for (i=0; i<L*L; i++)
 {
 cout<<S_d(i)<<endl;
 }
 cout<<endl;
 
 cout<<"S_ud\n";
 for (i=0; i<L*L; i++)
 {
 cout<<S_ud(i)<<endl;
 }
 cout<<endl;
 */

/*
 double tol=1e-14;
 for (i=0; i<L*L; i++)
 for (j=0; j<L*L; j++)
 {
 if (abs(W(i,j))<tol) W(i,j)=0;
 if (abs(X(i,j))<tol) X(i,j)=0;
 }
 
 cout<<"W:\n";
 cout<<W.cols(0,12);
 
 cout<<"X:\n";
 cout<<X.cols(0,12);
 */
/*
 mat evecVu;
 vec evalVu;
 
 eig_sym(evalVu,evecVu, Vu,"std");
 
 cout<<setiosflags(ios::left);
 cout<<"eigenvalues:\n";
 for (i=0; i<L*L; i++)
 {
 cout<<evalVu(i)<<endl;
 }
 cout<<endl;
 */

//	cout<<"eigenvectors:\n";
//	cout<<evecVu.cols(0,12);

//}

#endif /* Hamiltonian_h */
