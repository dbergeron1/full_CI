//
//  sparse_matrix.h
//  project STTN
//
//  Copyright (C) Dominic Bergeron 2020
//

#ifndef sparse_matrix_h
#define sparse_matrix_h

#include <vector>
#include <map>
#include "misc.h"

//template<class ind_T, class val_T>
//using sparse_matrix_map = vector<map<ind_T, val_T>>;

template<class ind_T, class val_T>
class matrix_element
{
public:
	explicit matrix_element(ind_T indp=0, val_T vp=0){index=indp; value=vp;}
	matrix_element(const matrix_element &me){index=me.index; value=me.value;}
	
	ind_T index;
	val_T value;
};

template<class ind_T, class val_T>
using sparse_matrix_vec = vector<vector<matrix_element<ind_T, val_T>>>;

template<class ind_T, class val_T, class cfs_T>
void prod_mat_vec(vector<cfs_T> &M_v, const sparse_matrix_vec<ind_T, val_T> &M, const vector<cfs_T> &v)
{
	if (v.size()!=M.size())
	{
		cerr<<"prod_mat_vec() error\n";
		return;
	}
	M_v.clear();
	M_v.resize(v.size());
	
	ind_T i,j;
	cfs_T val;
	for (i=0; i<v.size(); i++)
	{
		val=0;
		for (j=0; j<M[i].size(); j++)
		{
			val+=M[i][j].value*v[M[i][j].index];
		}
		M_v[i]=val;
	}
	
}



#endif /* sparse_matrix_h */
