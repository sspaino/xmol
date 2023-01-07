#include "scf.h"






SCF::SCF(const MOLECULE &mol):_mol(mol)
{
    _S = _mol._S;
    _T = _mol._T;
    _V = _mol._V;
    _EE = _mol._EE;
    _nbf = mol._nbf;
    _occ = mol._occ;
    _enuc = mol._enuc;



    _PO = new double[_nbf * _nbf]();
    _P  = new double[_nbf * _nbf]();
    _FO = new double[_nbf * _nbf]();
    _F0 = new double[_nbf * _nbf]();
    _F  = new double[_nbf * _nbf]();
    _CO = new double[_nbf * _nbf]();
    _C0 = new double[_nbf * _nbf]();
    _C = new double[_nbf * _nbf]();
    _eo = new double[_nbf * _nbf]();
    _e  = new double[_nbf * _nbf]();
    _X = orthogonalize(_S,_nbf);
    _Hcore = new double[_nbf * _nbf]();
    for (int i = 0; i < _nbf * _nbf; i++){_Hcore[i] = _T[i] + _V[i];}


    for (int step = 0; step < 100; step++)
    {
        for (int i = 0; i < _nbf * _nbf; i++){_FO[i] = _F[i];}
        bulidfock(_F,_Hcore,_P,_EE,_nbf);//add G to the core hanmiltonian to bulid fock matrix.
        if ((step > 0 && _ndiis > 0) && _isdiis)
        {
            diis();
        }
        for (int i = 0; i < _nbf * _nbf; i++){_F0[i] = _F[i];}
        transform(_F0,_X,_nbf); //F' = XFX;
        for (int i = 0; i < _nbf; i++){_eo[i] = _e[i];}
        for (int i = 0; i < _nbf * _nbf; i++){_CO[i] = _C[i];}
        for (int i = 0; i < _nbf * _nbf; i++){_C0[i] = _F0[i];}
        eigen(_C0,_e,_nbf); //diagonlize F'
        mmult(_X,_C0,_C,_nbf);
        for (int i = 0; i < _nbf * _nbf; i++){_PO[i] = _P[i];}
        bulidP(_P,_C,_occ,_nbf); // form a new density mmatrix P
        _EO = _E;
        _E = calc_energy(_P,_Hcore,_F,_nbf) + _enuc;
        _dE = _E - _EO;
        if (step > 1 && fabs(_dE ) < _conver)
        {
            printf("SCF Done After cycle %d.\n",step);
            printf("E(RHF)  = %20.8lf a.u\n",_E);
            printf("dE(RHF) = %20.8lf a.u\n",_dE);
            break;
        }
        else if (step == (_maxiter - 1))
        {
            printf ("CONVERGENCE FAILED\n");
        }

    }




}





void SCF::diis()
{


    _f   = new double[_nbf * _nbf](); // copy F to f.
    _err = new double[_nbf * _nbf]();
    _FP  = new double[_nbf * _nbf]();
    _FPS = new double[_nbf * _nbf]();
    _SP  = new double[_nbf * _nbf]();
    _SPF = new double[_nbf * _nbf]();
    memcpy(_f,_F,_nbf * _nbf * sizeof(double));
    mmult_sym(_F,_P,_FP,_nbf);
    mmult(_FP,_S,_FPS,_nbf);
    mmult_sym(_S,_P,_SP,_nbf);
    mmult(_SP,_F,_SPF,_nbf);


    double t = 0.;
    for (int i = 0; i < _nbf*_nbf; i++)
    {
        _err[i] = _SPF[i] - _FPS[i];
        if (fabs(_SPF[i] - _FPS[i] > t))
        {
            t = fabs(_SPF[i] - _FPS[i]);
        }
    }

    if (_errset.size() >= _ndiis)
    {
        _errset.erase(_errset.begin());
        _fockset.erase(_fockset.begin());
    }
    _errset.push_back(_err);
    _fockset.push_back(_f);
    
    int N = _errset.size();
    double *bmatrix = new double[(N+1) * (N+1)];
    double *rhs =  new double[(N+1)];
    memset(bmatrix, 0, (N+1) * (N+1) * sizeof(double));
    memset(rhs, 0, (N+1) * sizeof(double));

    for (int i = 0; i < N; i++)
    {
        for (int j = i; j < N; j++)
        {
            double val = 0.;
            for (int k = 0; k < _nbf * _nbf; k++)
            {
                val += _errset[i][k] * _errset[j][k];
            }
            bmatrix[i * (N + 1) + j] = val;
            bmatrix[j * (N + 1) + i] = val;
        }
    }




    for (int i = 0; i < N; i++)
    {
        bmatrix[i * (N + 1) + N] = -1.;
        bmatrix[N * (N + 1) + i] = -1.;
        rhs[i] = 0.;
    }
    bmatrix[N * (N + 1) + N] = 0.;
    rhs[N] = -1.;

    if (bmatrix[0] > 1.) 
    {
        double factor = 1. / bmatrix[0];
        for (int i = 0; i < N; i++) 
        {
            for (int j = 0; j < N; j++) 
            {
                bmatrix[i * (N + 1) + j] *= factor;
            }
        }
    }



    double* coeff = new double[(N+1)]();
    int info = solve_AxB (bmatrix, rhs, coeff, N + 1);

    if (info == 0) 
    {
        memset(_F, 0, _nbf * _nbf * sizeof(double));
        if (N == 1){coeff[0] = 1.;}
        
        for (int i = 0; i < N; i++) 
        {
            for (int k = 0; k < _nbf * _nbf; k++)
            {   

                _F[k] += coeff[i] * _fockset[i][k];   
            }
        }
    }
    else
    {
        for (int k = 0; k <  _nbf * _nbf; k++)
        {
            _F[k] = _f[k];
        }
    }




}







SCF::~SCF()
{

    delete[] _PO; 
    delete[] _P; 
    delete[] _FO; 
    delete[] _F0; 
    delete[] _F; 
    delete[] _CO; 
    delete[] _C0;
    delete[] _C; 
    delete[] _eo; 
    delete[] _e; 
    delete[] _f;
    delete[] _err;
    delete[] _FP;
    delete[] _FPS;
    delete[] _SP;
    delete[] _SPF;

}







int main()
{
    const std::string filename1 = "../../example/H2O.in";
    const std::string filename2 = "../../example/N2.in";
    const std::string filename3 = "../../example/BZ.in";
    //const std::string filename = "test1.in";

    MOLECULE mol1(filename1);
    MOLECULE mol2(filename2);
    MOLECULE mol3(filename3);

    std::cout << std::endl;
    std::cout << "test1*********************************************" << std::endl;
    mol1.read();
    mol1.bulid();
    SCF hf1(mol1);
    std::cout << std::endl;
    std::cout << "test2*********************************************" << std::endl;
    mol2.read();
    mol2.bulid();
    SCF hf2(mol2);
    std::cout << std::endl;
    std::cout << "test3*********************************************" << std::endl;
    mol3.read();
    mol3.bulid();
    SCF hf3(mol3);
    return 0;
}