#include <vector>
#include "utils.h"
#include "molecule.h"


class SCF
{


    public:
        SCF(const MOLECULE &mol);
        ~SCF();



    private:
        const MOLECULE &_mol;
        //scf============
        const int _maxiter = 100;
        const double _conver = 1e-8;
        const double _P_RMS_conver = 1e-8;
        bool _isdiis = true; // open diis;
        int _ndiis = 8; //diis space;
        int _nbf;
        //bool isdiis = false;

        
        double _EO = 0;//old energy 
        double _E = 0; // new energy ;
        double _dE; //energy gap;
        double _P_RMS = 0.;//  

        

        double* _S;
        double* _T;
        double* _V;
        double* _EE;
        int *_occ;
        double _enuc;

        double* _PO; //old density matrix;
        double* _P; //new density matrix

        double* _FO; //old Fock matrix
        double* _F0;
        double* _F; //new Fock matrix
        
        double* _CO;
        double* _C0;
        double* _C;
        
        double* _eo;
        double* _e;

        double *_X;//diagonlize the overlap matrix
        double *_Hcore;
        //scf============



        //diis============
        std::vector<double*> _errset; // error set
        std::vector<double*> _fockset; //fock set 
        double *_f;
        double *_err;
        double *_FP;
        double *_FPS;
        double *_SP;
        double *_SPF;
        void diis();
        //diis============















};