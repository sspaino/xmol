#include <iostream>
#include <cmath>
#include <cstring>

extern "C" 
{
/*
lapack and blas from fortran version.
*/
int dsyev_ (char *jobz, char *uplo, int *n, double *a, int *lda, double *w,double *work, int *lwork, int *info);
int dsymm_ (char *side, char *uplo, int *m, int *n, double *alpha, double *a,int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);
int dgemm_ (char *transa, char *transb, int *m, int * n, int *k, double *alpha,double *a, int *lda, double *b, int *ldb, double *beta, double *c,int *ldc);
int dgemv_ (char *trans, int *m, int *n, double * alpha, double *a, int *lda,double *x, int *incx, double *beta, double *y, int *incy);
int dgesv_ (int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b,int* ldb, int * info);
double ddot_ (int *n, double *dx, int *incx, double *dy, int *incy);
}


#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))

int solve_AxB (double * a, double *b, double * d, int n);
void mmult(double *a,double *b, double *c ,int n);
void mmult_sym (double *a, double *b, double * c, int n);
double *matrix_add(double *a, double *b,int nbf);
int eigen(double *a,double *eval,int n);
void transform (double *a, double *b, int n);
double *orthogonalize(double *ovlap,int nbf);
void contribution_4ind (int i, int n, double t, double* Fock, double* P);
void contribution_3ind (int i, int l, int n, double t, double* Fock, double* P);
void contribution_2indpairA (int i, int k, int n, double t, double* Fock,double* P);
void contribution_2indpairB (int i, int j, int n, double t, double* Fock,double* P);
void contribution_2indA (int i, int k, int l, int n, double t, double* Fock,double* P);
void contribution_2indB (int i, int j, int l, int n, double t, double* Fock,double* P);
void contribution_noind (int i, int j, int k, int l, int n, double t, double* Fock,double* P);
int packed_index(int i ,int j);
void bulidfock(double *F, double *Hcore,double *P,double *EE,int nbf);
void bulidP(double *P,double *C,int *occ,int nbf);
double calc_energy(double *P,double *Hcore,double *F,int nbf);
