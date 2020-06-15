/*
file misc.h

copyright (C) 2020 Dominic Bergeron (dominic.bergeron@usherbrooke.ca)

miscellanous functions
*/

#include "includeDef.h"

#ifndef MISC_H
#define MISC_H

extern "C"
{
	//DGESV computes the solution to a real system of linear equations A * X = B where A is an N-by-N matrix and X and B are N-by-NRHS matrices
	void dgesv_(int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B, int *LDB, int *INFO);
	//computes all eigenvalues of a symmetric tridiagonal matrix using the Pal-Walker-Kahan variant of the QL or QR algorithm.
	// http://www.netlib.org/lapack/explore-html/d9/df2/dsterf_8f.html
	void dsterf_(int *N, double *D, double *E, int *INFO );
	//computes all eigenvalues and, optionally, eigenvectors of a symmetric tridiagonal matrix using the implicit QL or QR method.
	// http://www.netlib.org/lapack/explore-html/d9/d3f/dsteqr_8f.html
	void dsteqr_(char *COMPZ, int *N, double *D, double *E, double *Z, int *LDZ, double *WORK, int *INFO );
	//computes all eigenvalues and, optionally, eigenvectors of a (real) symmetric matrix
	// http://www.netlib.org/lapack/explore-html/dd/d4c/dsyev_8f_source.html
	void dsyev_(char *JOBZ, char *UPLO, int *N, double *A, int *LDA, double *W, double *WORK, int *LWORK, int *INFO );
	//computes the eigenvalues and, optionally, the left and/or right eigenvectors of a (complex) Hermitian matrix
	// http://www.netlib.org/lapack/explore-html/d6/dee/zheev_8f_source.html
	void zheev_(char *JOBZ, char *UPLO, int *N, dcomplex *A, int *LDA, double *W, dcomplex *WORK, int *LWORK, double *RWORK, int *INFO);
}

double dval(mpf_class v) {return v.get_d();}
double dval(double v) {return v;}
double dval(float v) {return (double)v;}

//solve the linear system AX=B for X, where B is a matrix
template<class elem_T>
void lin_solve(const vector<vector<elem_T>> &A, const vector<vector<elem_T>> &B, vector<vector<elem_T>> &X)
{
	int NRHS[1];
	NRHS[0]=B[0].size();
	int NA[1];
	NA[0]=A.size();
	int *IPIV=new int[NA[0]];
	int INFO[1];
	
	double *Av=new double[NA[0]*NA[0]];
	double *Bv=new double[NA[0]*NRHS[0]];
	int i,j;
	for (i=0; i<NA[0]; i++)
	{
		for (j=0; j<NA[0]; j++)	Av[i+j*NA[0]]=dval(A[i][j]);
		for (j=0; j<NRHS[0]; j++) Bv[i+j*NA[0]]=dval(B[i][j]);
	}
	
	dgesv_(NA, NRHS, Av, NA, IPIV, Bv, NA, INFO);
	
	if (!INFO[0])
	{
		X.clear();
		X.resize(NA[0]);
		for (i=0; i<NA[0]; i++)
		{
			X[i].resize(NRHS[0]);
			for (j=0; j<NRHS[0]; j++) X[i][j]=Bv[i+j*NA[0]];
		}
	}
	else
	{
		cerr<<"lin_solve() error: no solution found\n";
	}
	
	if (IPIV) delete [] IPIV;
	if (Av) delete [] Av;
	if (Bv) delete [] Bv;
}

//solve the linear system AX=B for X, where B is a vector
template<class elem_T>
void lin_solve(const vector<vector<elem_T>> &A, vector<elem_T> &B, vector<elem_T> &X)
{
	int NRHS[1];
	NRHS[0]=1;
	int NA[1];
	NA[0]=A.size();
	int *IPIV=new int[NA[0]];
	int INFO[1];
	
	double *Av=new double[NA[0]*NA[0]];
	double *Bd=new double[NA[0]];
	int i,j;
	for (i=0; i<NA[0]; i++)
	{
		Bd[i]=dval(B[i]);
		for (j=0; j<NA[0]; j++)	Av[i+j*NA[0]]=dval(A[i][j]);
	}
	
	dgesv_(NA, NRHS, Av, NA, IPIV, Bd, NA, INFO);
	
	if (!INFO[0])
	{
		X.clear();
		X.resize(NA[0]);
		for (i=0; i<NA[0]; i++)
		{
			X[i]=Bd[i];
		}
	}
	else
	{
		cerr<<"lin_solve() error: no solution found\n";
	}
	
	if (IPIV) delete [] IPIV;
	if (Av) delete [] Av;
	if (Bd) delete [] Bd;
}

float conj(float v)
{
	return v;
}

double conj(double v)
{
	return v;
}

mpf_class conj(mpf_class v)
{
	return v;
}

template<class ind_T, class cfs_T>
void scalar_prod_vec(cfs_T &pr, const vector<cfs_T> &v1, const vector<cfs_T> &v2)
{
	pr=0;
	
	for (ind_T i=0; i<v1.size(); i++) pr+=conj(v1[i])*v2[i];
}

template<class val_T, class cfs_T>
void prod_dense_mat_vec(vector<cfs_T> &M_v, const vector<vector<val_T>> &M, const vector<cfs_T> &v)
{
	if (v.size()!=M.size())
	{
		cerr<<"prod_mat_vec() error\n";
		return;
	}
	
	int s=v.size();
	M_v.resize(s);
	int i,j;
	val_T val;
	for (i=0; i<s; i++)
	{
		val=0;
		for (j=0; j<s; j++)
		{
			val+=M[i][j]*v[j];
		}
		M_v[i]=val;
	}
	
}

template<class ind_T, class cfs_T>
void vec_plus_equal(vector<cfs_T> &v1, const vector<cfs_T> &v2, cfs_T c=1)
{
	if (v1.size()!=v2.size())
	{
		cout<<"vec_plus_equal() error\n";
		return;
	}
	
	for (ind_T i=0; i<v1.size(); i++)
	{
		v1[i]+=c*v2[i];
	}
}

template<class ind_T, class cfs_T>
void prod_const_vec(vector<cfs_T> &c_v, const cfs_T &c, const vector<cfs_T> &v)
{
	c_v.resize(v.size());
	
	for (ind_T i=0; i<v.size(); i++)
	{
		c_v[i]=c*v[i];
	}
}

template<class ind_T, class cfs_T>
void vec_times_equal(vector<cfs_T> &v, const cfs_T &c)
{
	for (ind_T i=0; i<v.size(); i++)
	{
		v[i]=c*v[i];
	}
}

template<class ind_T, class cfs_T>
void vec_div_equal(vector<cfs_T> &v, const cfs_T &c)
{
	for (ind_T i=0; i<v.size(); i++)
	{
		v[i]=v[i]/c;
	}
}

template<class ind_T, class cfs_T>
void sum_vec(vector<cfs_T> &v_sum, const vector<cfs_T> &v1, const vector<cfs_T> &v2, cfs_T c1=1, cfs_T c2=1)
{
	if (v1.size()!=v2.size())
	{
		cout<<"sum_vec() error\n";
		return;
	}
	v_sum.resize(v1.size());
	
	for (ind_T i=0; i<v1.size(); i++)
	{
		v_sum[i]=c1*v1[i]+c2*v2[i];
	}
}

template<class elem_T>
void invert_matrix(const vector<vector<elem_T>> &A, vector<vector<elem_T>> &invA)
{
	int i;
	
	vector<vector<elem_T>> Id(A.size());
	for (i=0; i<A.size(); i++)
	{
		Id[i].resize(A.size(),0);
		Id[i][i]=1.0;
	}
	
	lin_solve(A, Id, invA);
}

//matrix prod C=A*B
template<class elem_T>
void matrix_prod(const vector<vector<elem_T>> &A, const vector<vector<elem_T>> &B, vector<vector<elem_T>> &C)
{
	int i,j,k;
	elem_T v;
	
	C.clear();
	C.resize(A.size());
	for (i=0; i<A.size(); i++)
	{
		C[i].resize(B[0].size(),0);
		for (j=0; j<B[0].size(); j++)
		{
			for (k=0; k<B.size(); k++)
			{
				C[i][j]+=A[i][k]*B[k][j];
			}
		}
	}
}

//matrix sum C=cA*A+cB*B
template<class elem_T>
void matrix_sum(elem_T cA, const vector<vector<elem_T>> &A, elem_T cB, const vector<vector<elem_T>> &B, vector<vector<elem_T>> &C)
{
	int i,j,k;
	elem_T v;
	
	C.clear();
	C.resize(A.size());
	for (i=0; i<A.size(); i++)
	{
		C[i].resize(A[0].size());
		for (j=0; j<B[0].size(); j++)
		{
			C[i][j]=cA*A[i][j]+cB*B[i][j];
		}
	}
}

template<class elem_T>
void transpose_matrix(const vector<vector<elem_T>> &A, vector<vector<elem_T>> &At)
{
	At.resize(A[0].size());
	int i,j;
	for (i=0; i<A[0].size(); i++)
	{
		At[i].resize(A.size());
		for (j=0; j<A.size(); j++)
		{
			At[i][j]=A[j][i];
		}
	}
}

template<class elem_T>
void print_matrix(const vector<vector<elem_T>> &A)
{
	int i, j;
	
	cout<<setiosflags(ios::left);
	cout<<setprecision(6);
	for (i=0; i<A.size(); i++)
	{
		for (j=0; j<A[i].size(); j++)
		{
			cout<<setw(14)<<A[i][j];
		}
		cout<<endl;
	}
}

template<class elem_T>
void Gram_Schmidt_orthog(const vector<vector<elem_T>> &A_in, vector<vector<elem_T>> &A_out)
{
	int i,j,k;
	
	int Nrows=A_in.size();
	int Ncols=A_in[0].size();
	
	vector<vector<elem_T>> U(Nrows);
	A_out.clear();
	A_out.resize(Nrows);
	for (i=0; i<Nrows; i++)
	{
		U[i].resize(Ncols,0);
		A_out[i].resize(Ncols,0);
	}
	
	elem_T sp_u, sp_u_v;
	vector<elem_T> u(Nrows), u_c(Nrows), v(Nrows), sp(Ncols);
	
	//u_c=A_in[:][0]
	//U[:][0]=u_c
	for (i=0; i<Nrows; i++)
	{
		u_c[i]=A_in[i][0];
	}
	//sp[0]=<u_c|u_c>
	scalar_prod_vec<int,elem_T>(sp[0], u_c, u_c);
	//A_out[:][0]=u_c/<u_c|u_c>;
	for (i=0; i<Nrows; i++)
	{
		U[i][0]=u_c[i];
		A_out[i][0]=U[i][0]/sqrt(sp[0]);
	}
	
	for (j=1; j<Ncols; j++)
	{
		//v=A_in[:][j]
		//u_c=v
		for (i=0; i<Nrows; i++)
		{
			v[i]=A_in[i][j];
			u_c[i]=v[i];
		}
		for (k=0; k<j; k++)
		{
			//u=U[:][k]
			for (i=0; i<Nrows; i++)
			{
				u[i]=U[i][k];
			}
			//sp_u=<u|u>
			sp_u=sp[k];
			//sp_u_v=<v|u>
			scalar_prod_vec<int,elem_T>(sp_u_v,v,u);
			//u_c=u_c-u*(<v|u>/<u|u>)
			vec_plus_equal<int,elem_T>(u_c, u, -sp_u_v/sp_u);
		}
		scalar_prod_vec<int,elem_T>(sp[j],u_c,u_c);
		for (i=0; i<Nrows; i++)
		{
			U[i][j]=u_c[i];
			A_out[i][j]=u_c[i]/sqrt(sp[j]);
		}
	}
}

float real(float v)
{
	return v;
}
double real(double v)
{
	return v;
}
mpf_class real(mpf_class v)
{
	return v;
}

mpf_class foor(mpf_class v)
{
	mpf_class tmp;
	mpf_floor(tmp.get_mpf_t(),v.get_mpf_t());
	return tmp;
}

mpf_class ceil(mpf_class v)
{
	mpf_class tmp;
	mpf_ceil(tmp.get_mpf_t(),v.get_mpf_t());
	return tmp;
}

mpf_class fabs(mpf_class v)
{
	if (v>=0) return v;
	else return -v;
}

template<class val_T>
val_T abs(val_T v)
{
	if (v>=0) return v;
	else return -v;
}

mpf_class round(mpf_class v)
{
	mpf_class tmp1=ceil(v);
	mpf_class tmp2=floor(v);
	if (fabs(tmp2-v)<fabs(tmp1-v)) return tmp2;
	else return tmp1;
}

mpf_class sqrt(mpf_class z)
{
	mpf_class tmp;
	mpf_sqrt(tmp.get_mpf_t(),z.get_mpf_t());
	
	return tmp;
}

mpz_class mp_factorial(unsigned long n)
{
	mpz_class fn;
	mpz_fac_ui(fn.get_mpz_t(), n);
	return fn;
}

mpz_class mp_binom(unsigned long N, unsigned long n)
{
	mpz_class b=N;
	
	for (unsigned i=1; i<n; i++)
	{
		b*=(N-i);
	}
	
	return b/mp_factorial(n);
}

unsigned long binom_cf(unsigned N, unsigned n)
{
	mpz_class b=N;
	
	for (unsigned i=1; i<n; i++)
	{
		b*=(N-i);
	}
	
	mpz_class mp_cf=b/mp_factorial(n);
	
	return mpz_get_ui(mp_cf.get_mpz_t());
}


void shift_left(unsigned &s, unsigned Nb) {s=s<<Nb;}
void shift_left(uL &s, unsigned Nb) {s=s<<Nb;}
void shift_left(uLL &s, unsigned Nb) {s=s<<Nb;}
void shift_left(mpz_class &s, unsigned Nb) {mpz_mul_2exp(s.get_mpz_t(), s.get_mpz_t(), Nb);}
void shift_right(unsigned &s, unsigned Nb) {s=s>>Nb;}
void shift_right(uL &s, unsigned Nb) {s=s>>Nb;}
void shift_right(uLL &s, unsigned Nb) {s=s>>Nb;}
void shift_right(mpz_class &s, unsigned Nb) {mpz_tdiv_q_2exp(s.get_mpz_t(), s.get_mpz_t(), Nb);}

void setbit(unsigned &s, unsigned l)
{
	s=s|((unsigned)1<<l);
}

void setbit(uL &s, unsigned l)
{
	s=s|((uL)1<<l);
}

void setbit(uLL &s, unsigned l)
{
	s=s|((uLL)1<<l);
}

void setbit(mpz_class &s, unsigned l)
{
	mpz_setbit(s.get_mpz_t(), l);
}

void clrbit(unsigned &s, unsigned l)
{
	uL stmp=(unsigned)1<<l;
	s=s&(~stmp);
}

void clrbit(uL &s, unsigned l)
{
	uL stmp=(uL)1<<l;
	s=s&(~stmp);
}

void clrbit(uLL &s, unsigned l)
{
	uLL stmp=(uLL)1<<l;
	s=s&(~stmp);
}

void clrbit(mpz_class &s, unsigned l)
{
	mpz_clrbit(s.get_mpz_t(), l);
}

int tstbit(unsigned s, unsigned l)
{
	return (s>>l)&1;
}

int tstbit(uL s, unsigned l)
{
	return (s>>l)&1;
}

int tstbit(uLL s, unsigned l)
{
	return (s>>l)&1;
}

int tstbit(mpz_class s, unsigned l)
{
	return mpz_tstbit(s.get_mpz_t(), l);
}

uL int_div(uL s1, uL s2)
{
	return s1/s2;
}

uLL int_div(uLL s1, uLL s2)
{
	return s1/s2;
}

mpz_class int_div(mpz_class s1, mpz_class s2)
{
	mpz_class s;
	mpz_tdiv_q(s.get_mpz_t(), s1.get_mpz_t(), s2.get_mpz_t());
	return s;
}

unsigned mod(unsigned s1, unsigned s2)
{
	return s1%s2;
}

uL mod(uL s1, uL s2)
{
	return s1%s2;
}

uLL mod(uLL s1, uLL s2)
{
	return s1%s2;
}

mpz_class mod(mpz_class s1, mpz_class s2)
{
	mpz_class s;
	mpz_tdiv_r(s.get_mpz_t(), s1.get_mpz_t(), s2.get_mpz_t());
	return s;
}

unsigned mod_p2(unsigned s1, unsigned p)
{
	if (p<8*sizeof(unsigned))
		return s1&(((unsigned)1<<p)-1);
	else
		return s1;
}

uL mod_p2(uL s1, unsigned p)
{
	if (p<8*sizeof(uL))
		return s1&(((uL)1<<p)-1);
	else
		return s1;
}

uLL mod_p2(uLL s1, unsigned p)
{
	if (p<8*sizeof(uLL))
		return s1&(((uLL)1<<p)-1);
	else
		return s1;
}

mpz_class mod_p2(mpz_class s1, unsigned p)
{
	mpz_class z;
	mpz_tdiv_r_2exp(z.get_mpz_t(), s1.get_mpz_t(), p);
	return z;
}

uL int_div_p2(uL s1, unsigned p)
{
	if (p<8*sizeof(uL))
		return s1/((uL)1<<p);
	else
		return 0;
}

uLL int_div_p2(uLL s1, unsigned p)
{
	if (p<8*sizeof(uLL))
		return s1/((uLL)1<<p);
	else
		return 0;
}

mpz_class int_div_p2(mpz_class s1, unsigned p)
{
	mpz_class z;
	mpz_tdiv_q_2exp(z.get_mpz_t(), s1.get_mpz_t(), p);
	return z;
}

uL bin_AND(uL s1, uL s2)
{
	return s1&s2;
}

uLL bin_AND(uLL s1, uLL s2)
{
	return s1&s2;
}

mpz_class bin_AND(mpz_class s1, mpz_class s2)
{
	mpz_class s;
	mpz_and(s.get_mpz_t(), s1.get_mpz_t(), s2.get_mpz_t());
	return s;
}

uL bin_OR(uL s1, uL s2)
{
	return s1|s2;
}

uLL bin_OR(uLL s1, uLL s2)
{
	return s1|s2;
}

mpz_class bin_OR(mpz_class s1, mpz_class s2)
{
	mpz_class s;
	mpz_ior(s.get_mpz_t(), s1.get_mpz_t(), s2.get_mpz_t());
	return s;
}

uL bin_XOR(uL s1, uL s2)
{
	return s1^s2;
}

uLL bin_XOR(uLL s1, uLL s2)
{
	return s1^s2;
}

mpz_class bin_XOR(mpz_class s1, mpz_class s2)
{
	mpz_class s;
	mpz_xor(s.get_mpz_t(), s1.get_mpz_t(), s2.get_mpz_t());
	return s;
}

template<class st_T>
void print_binary_string(const st_T &s, int lgth)
{	
//	cout<<setiosflags(ios::left);
	for (int j=lgth-1; j>=0; j--)
	{
		cout<<tstbit(s,j);
	}
//	cout<<endl;
}

uL bin_comp(uL s, unsigned n)
{
	return mod_p2(~s,n);
}

uLL bin_comp(uLL s, unsigned n)
{
	return mod_p2(~s,n);
}

mpz_class bin_comp(mpz_class s, unsigned n)
{
	mpz_class p2n=1;
	shift_left(p2n,n);
	return p2n-1-s;
//	mpz_com(ns.get_mpz_t(), s.get_mpz_t());
//	return ns;
}

double get_double(double v){return v;}
double get_double(mpf_class v){return v.get_d();}
double get_double(mpq_class v){return v.get_d();}

double conv_to_double(unsigned v){return (double)v;}
double conv_to_double(uL v){return (double)v;}
double conv_to_double(uLL v){return (double)v;}
double conv_to_double(mpz_class v){return v.get_d();}

unsigned get_ui(unsigned v){return v;}
unsigned get_ui(uL v){return (unsigned)v;}
unsigned get_ui(uLL v){return (unsigned)v;}
unsigned get_ui(mpz_class v){return v.get_ui();}


unsigned add_bits(const unsigned &s, unsigned startbit, unsigned endbit)
{
	//	if (endbit<startbit) return 0;
	
	unsigned s1=mod_p2(s,endbit+1);
	
	shift_right(s1,startbit);
	
	return __builtin_popcount(s1);
}

unsigned add_bits(const uL &s, unsigned startbit, unsigned endbit)
{
//	if (endbit<startbit) return 0;
		
	uL s1=mod_p2(s,endbit+1);
	
	shift_right(s1,startbit);
	
	return __builtin_popcountl(s1);
}

unsigned add_bits(const uLL &s, unsigned startbit, unsigned endbit)
{
//	if (endbit<startbit) return 0;
	
	uLL s1=mod_p2(s,endbit+1);
	
	shift_right(s1,startbit);
	
	return __builtin_popcountll(s1);
}

unsigned add_bits(const mpz_class &s, unsigned startbit, unsigned endbit)
{
//	if (endbit<startbit) return 0;
	
	mpz_class s1=mod_p2(s,endbit+1);
	
	shift_right(s1,startbit);
	
	return mpz_popcount(s1.get_mpz_t());
}

void remove_spaces_front(string &str)
{
	int j=0;
	while (str[j]==' ' || str[j]=='\t') j++;
	str=str.substr(j);
}

void remove_spaces_back(string &str)
{
	int j=str.size()-1;
	while (str[j]==' ' || str[j]=='\t') j--;
	str.resize(j+1);
}

void remove_spaces_ends(string &str)
{
	remove_spaces_front(str);
	remove_spaces_back(str);
}

#endif
