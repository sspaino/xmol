#include <fstream>
#include <iostream>
#include <vector>
#include <regex>
#include <cstring>
#include <cmath>
#include <ctime>
#include <algorithm>
#include "chemutils.h"



extern "C"
{
#include <cint.h>
int cint1e_ovlp_sph (double *buf, int *shls, int *atm, int natm, int *bas,int nshl, double *env);
int cint1e_kin_sph (double *buf, int *shls, int *atm, int natm, int *bas,int nshl, double *env);
int cint1e_nuc_sph (double *buf, int *shls, int *atm, int natm, int *bas,int nshl, double *env);
FINT cint2e_sph(double *opijkl, FINT *shls,FINT *atm, FINT natm, FINT *bas, FINT nshl, double *env,CINTOpt *opt);
}





#define MAX_ENV 10000 
#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))




typedef struct PRIM
{

    int symbol;
    int shell;
    double exp;
    double coeff;
} PRIM;


typedef struct CTRG
{
    int shell;
    int nprim;
    int nctr;
} CTRG;


typedef struct ATOM
{
    std::string symbol;
    int charge;
    double mass;
    std::vector<double> coor;
    std::vector <CTRG> ctrgs;
    std::vector <PRIM> prims;
} ATOM;


class MOLECULE
{
    public:
        int _charge;
        double _spin;
        int _nelec;
        int _nocc;
        int _nbasis;
        int _natm;
        int _nshell;
        int _nbf;
        std::string _input;
        std::string _basis_set = "sto-3g";
        std::string _basis_filename = "../../basis/" + _basis_set + ".0.gbs";
        std::vector<ATOM*> _atoms;
        int* _atm;
        int* _bas;
        double* _env;
        int* _shell_label;
        double* _S;
        double* _T;
        double* _V;
        double* _EE;
        int *_occ;
        double _enuc;
        void read();
        void bulid();
        MOLECULE(const std::string& input);
        ~MOLECULE();


    private:
        void basisload();
        void cut(const std::string & symbol,int *start , int *end);
        std::vector<std::string> split(const std::string &str, const std::string &delim);
        int qm(char shellname);
        static bool cmp_prim(const PRIM &a, const PRIM &b);
        static bool cmp_ctrg(const CTRG &a, const CTRG &b);
        int packed_index(int i ,int j);
        void normalize();
        void one_electron_integrals();
        void two_electron_integrals();



};


