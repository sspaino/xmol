#include "utils.h"





int solve_AxB (double * a, double *b, double * d, int n)
{
    int *ipiv = new int[n]();
    int nrhs = 1;
    int info;
    double *c = new double[n * n]();
    for (int i = 0; i < n * n; i++)
        c[i] = a[i];
    for (int i = 0; i < n; i++)
        d[i] = b[i];
    dgesv_ (&n, &nrhs, c, &n, ipiv, d, &n, &info);
    delete[] ipiv;
    delete[] c;
    return info;

}


double norm(double *a,double *b,int n)
{

    double *temp = new double[n]();
    char transa = 'N';
    double alpha = 1.,beta = 0;
    int inc = 1;
    dgemv_ (&transa, &n, &n, &alpha, a, &n, b, &inc, &beta, temp, &inc); //temp=a*b
    double morm_val = ddot_(&n,a,&inc,temp,&inc);
    delete[] temp;
    return morm_val;   
}


// matrix-matrix multiplication c=a b

void mmult(double *a,double *b, double *c ,int n)
{
    double alpha = 1.,beta = 0.;
    char transa = 'N';
    char transb = 'N';
    dgemm_ (&transa, &transb, &n, &n, &n, &alpha, a, &n, b, &n, &beta, c, &n); //c=a*b
}

// matrix-matrix multiplication c=a b , with symmetric matrix a
void mmult_sym (double *a, double *b, double * c, int n) 
{
    char transa = 'L';
    char transb = 'U';
    double alpha = 1, beta = 0.;
    dsymm_ (&transa,&transb, &n, &n, &alpha, a, &n, b, &n, &beta, c, &n);
}


//matrix addition 
double *matrix_add(double *a, double *b,int nbf)
{
    double* ret = new double[nbf * nbf];
    for (int i = 0; i < nbf; i++)
    {
        for (int j = 0; j < nbf; j++)
        {
           ret[j+i*nbf] = a[j+i*nbf] + b[j+i*nbf];
        }
    }
    return ret;
}

//eigenvalue and eigenvector
int eigen(double *a,double *eval,int n)
{
    int lwork, info;
    double* work;
    double opt_lwork = 0.0;
    lwork = -1;
    char transa = 'V';
    char transb = 'L';
    dsyev_ (&transa, &transb, &n, a, &n, eval, (double *) &opt_lwork, &lwork,&info);
    lwork = (int) opt_lwork;
    work = new double[lwork]();
    dsyev_ (&transa, &transb, &n, a, &n, eval, work, &lwork, &info);
    delete[] work;
    if (info < 0)
    {
        printf("argument.");
    }
    if (info > 0)
    {
        printf("convergence.");
    }
    return ((int) info);
}


void transform (double *a, double *b, int n)
{
    double *temp = new double[n*n]();
    char transa = 'N';
    char transb = 'N';
    double alpha = 1.,beta = 0.;
    dgemm_ (&transa, &transb, &n, &n, &n, &alpha, a, &n, b, &n, &beta, temp, &n); //temp=a*b
    transa = 'T';
    transb = 'N';
    dgemm_ (&transa, &transb, &n, &n, &n, &alpha, b, &n, temp, &n, &beta, a, &n); //a=b^t*temp
    delete[] temp;
}



/*
there are two method to get X form overlap matrix szabo'book p143~144.
1.symmetric orthogonalization:X = S^-1/2
2.canonical orthogonalization:X = Us^-1/2
here we use first method.
*/
double *orthogonalize(double *ovlap,int nbf)
{
    //copy ovlap
    
    double* eigvec = new double[nbf * nbf]();
    double* eval = new double[nbf]();
    memcpy(eigvec, ovlap, nbf * nbf * sizeof(double));

    //eigenvalue and eigenvector
    int state = eigen(eigvec,eval,nbf);
    //X
    double *olap_inv_sqrt = new double [nbf * nbf]();  
    for (int i = 0; i < nbf; i++)
    {
        for (int j = 0; j < nbf; j++)
        {
            for (int k = 0; k < nbf; k++)
            {
                olap_inv_sqrt[i * nbf + j] += eigvec[k * nbf + i] * eigvec[k * nbf + j] / sqrt (eval[k]);
            }
        }
    }
    delete[] eval;
    return olap_inv_sqrt;
}

void contribution_4ind (int i, int n, double t, double* Fock, double* P) 
{
  // 4 indices equal i,j,k,l
  Fock[i * n + i] += P[i * n + i] * t ;
}

void contribution_3ind (int i, int l, int n, double t, double* Fock, double* P) 
{
  // 3 indices equal i=j=k
  Fock[i * n + i] += P[i * n + l] * t * 2.0;
  Fock[i * n + l] += P[i * n + i] * t * 1.0;
  Fock[l * n + i] += P[i * n + i] * t * 1.0;
}

void contribution_2indpairA (int i, int k, int n, double t, double* Fock,double* P) 
{
  // 2 index pairs equal i=j k=l
  Fock[i * n + i] +=  P[k * n + k] * t * 2.0;
  Fock[k * n + k] +=  P[i * n + i] * t * 2.0;
  Fock[i * n + k] += -P[i * n + k] * t;
  Fock[k * n + i] += -P[k * n + i] * t;
}

void contribution_2indpairB (int i, int j, int n, double t, double* Fock,double* P) 
{
  // 2 index pairs equal i=k j=l
  Fock[i * n + j] +=  P[i * n + j] * t *3.0;
  Fock[j * n + i] +=  P[i * n + j] * t *3.0;
  Fock[i * n + i] += -P[j * n + j] * t;
  Fock[j * n + j] += -P[i * n + i] * t;
}

void contribution_2indA (int i, int k, int l, int n, double t, double* Fock,double* P) 
{
  // 2 indices equal i=j
  Fock[i * n + i] +=  P[k * n + l] * t * 4.0;
  Fock[k * n + l] +=  P[i * n + i] * t * 2.0;
  Fock[l * n + k] +=  P[i * n + i] * t * 2.0;
  Fock[i * n + k] += -P[i * n + l] * t;
  Fock[i * n + l] += -P[i * n + k] * t;
  Fock[k * n + i] += -P[l * n + i] * t;
  Fock[l * n + i] += -P[k * n + i] * t;
}

void contribution_2indB (int i, int j, int l, int n, double t, double* Fock,double* P) 
{
  // 2 indices equal i=k
  Fock[i * n + j] +=  P[i * n + l] * t * 3.0;
  Fock[j * n + i] +=  P[i * n + l] * t * 3.0;
  Fock[i * n + l] +=  P[i * n + j] * t * 3.0;
  Fock[l * n + i] +=  P[i * n + j] * t * 3.0;
  Fock[i * n + i] += -P[j * n + l] * t * 2.0;
  Fock[j * n + l] += -P[i * n + i] * t;
  Fock[l * n + j] += -P[i * n + i] * t;
}

void contribution_noind (int i, int j, int k, int l, int n, double t, double* Fock,double* P) 
{
  // no indices equal i!=k (j!=l)
  Fock[i * n + j] +=  P[k * n + l] * t * 4.0;
  Fock[j * n + i] +=  P[k * n + l] * t * 4.0;
  Fock[k * n + l] +=  P[i * n + j] * t * 4.0;
  Fock[l * n + k] +=  P[i * n + j] * t * 4.0;
  Fock[i * n + k] += -P[j * n + l] * t;
  Fock[j * n + k] += -P[i * n + l] * t;
  Fock[i * n + l] += -P[j * n + k] * t;
  Fock[j * n + l] += -P[i * n + k] * t;
  Fock[k * n + i] += -P[l * n + j] * t;
  Fock[k * n + j] += -P[l * n + i] * t;
  Fock[l * n + i] += -P[k * n + j] * t;
  Fock[l * n + j] += -P[k * n + i] * t;
}


int packed_index(int i ,int j)
{
      int min = MIN0(i, j);
      int max = MAX0(i, j);
      return (min + max * (max + 1) / 2);
}

void bulidfock(double *F, double *Hcore,double *P,double *EE,int nbf)
{

    for (int i = 0; i < nbf*nbf; i++)
    {
        F[i] = Hcore[i];
    }
    int ij, kl;
    int ijkl;
    for (int i = 0; i < nbf ; i++)
    {
      for (int j = 0; j <= i; j++)
        {
            ij = packed_index(i,j);
            for (int k = 0; k < nbf; k++)
            {
                for (int l = 0; l <= k; l++)
                {
                    kl = packed_index(k,l);
                    if (ij < kl)
                        continue;
                    ijkl = packed_index(ij,kl);
                    double t = EE[ijkl];
                    if (i == j && k == l && i == k) 
                    {
                      // 4 indices equal
                      contribution_4ind(i, nbf, t, F, P);
                    } 
                    else if (i == j && j == k) 
                    {
                      // 3 indices equal i,j,k
                      contribution_3ind(i, l, nbf, t, F, P);
                    } 
                    else if (i == j && j == l ) 
                    {
                      // 3 indices equal i,j,l
                      contribution_3ind(i, k, nbf, t, F, P);
                    } 
                    else if (i == k && k == l) 
                    {
                      // 3 indices equal i,k,l
                      contribution_3ind(i, j, nbf, t, F, P);
                    } 
                    else if (i == j && k == l) {
                      // two indices pairwise equal ij and kl
                      contribution_2indpairA(i, k, nbf, t, F, P);
                    } 
                    else if (i == k && j == l) 
                    {
                      // two indices pairwise equal ik and jl
                      contribution_2indpairB(i, j, nbf, t, F, P);
                    } 
                    else if (i == j) 
                    {
                      // two indices  equal i=j
                      contribution_2indA(i, k, l, nbf, t, F, P);
                    } 
                    else if (k == l) 
                    {
                      // two indices  equal k=l
                      contribution_2indA(k, i, j, nbf, t, F, P);
                    } 
                    else if (i == k) 
                    {
                      // two indices  equal i=k
                      contribution_2indB(i, j, l, nbf, t, F, P);
                    } 
                    else if (j == l) 
                    {
                      // two indices  equal j=l
                      contribution_2indB(j, i, k, nbf, t, F, P);
                    } 
                    else 
                    {
                      contribution_noind(i, j, k, l, nbf, t, F, P);
                    }//if                       
                }   
            }   
        }
    }
}


/*
szabo's book p139->3.145
*/
void bulidP(double *P,double *C,int *occ,int nbf)
{
    memset(P, 0, nbf * nbf * sizeof(double));
    for (int i = 0; i < nbf; i++)
    {
        for (int j = i; j < nbf; j++)
        {
            for (int k = 0; k < nbf; k++)
            {
                P[j+i*nbf] += 0.5 * occ[k] * C[k*nbf + i] * C[k*nbf + j];
            }
            P[i+j*nbf] = P[j+i*nbf];
        }
        
    }

}

double calc_energy(double *P,double *Hcore,double *F,int nbf)
{
    double* fp = new double[nbf * nbf]();
    double* hp = new double[nbf * nbf]();
    mmult_sym(F,P,fp,nbf);
    mmult_sym(Hcore,P,hp,nbf);
    double E = 0;
    for (int i = 0; i < nbf; i++)
    {
      E += fp[i * nbf + i] + hp[i * nbf + i];
    }
    delete[] fp;
    delete[] hp;
    return(E);
}