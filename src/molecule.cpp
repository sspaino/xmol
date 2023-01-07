#include "molecule.h"





//ATOM::ATOM(const std::string& symbol,const int& charge,const double& mass,const std::vector<double>& coor)
//:_symbol(symbol),_charge(charge),_mass(mass),_coor(coor){}


MOLECULE::MOLECULE(const std::string& input)
:_input(input){}

MOLECULE::~MOLECULE()
{

    delete[] _atm;
    delete[] _bas;   
    delete[] _env;
    delete[] _shell_label;
    
    //delete[] _S;
    //delete[] _T;
    //delete[] _V;
    //delete[] _EE;

    //free ATOM
    for (auto it = _atoms.begin(); it != _atoms.end(); it++)
    {
	    if (*it != NULL)
	    {
		    delete *it;
		    *it = NULL;
		}
	}
}




void MOLECULE::cut(const std::string & symbol,int *start , int *end)
{
    std::ifstream file1(_basis_filename); 
    if(!file1){std::cout <<"open input file error!"<< std::endl;exit(0);}
    std::string line;
    std::regex regex_symbol("^"+symbol+"\\s+0$");
    std::smatch match_symbol;
    int i = 0;
    while (getline(file1,line))
    {

        if (regex_search(line,match_symbol,regex_symbol))
        {
            break;
        }
        i++;
    }
    file1.close();

    std::ifstream file2(_basis_filename); 
    if(!file2){std::cout <<"open input file error!"<< std::endl;exit(0);}
    std::regex regex_split("\\*\\*\\*\\*");
    std::smatch match_split;
    int j = 0;
    while (getline(file2,line))
    {
        //std::cout << line << std::endl;
        if (regex_search(line,match_split,regex_split))
        {
            if (j > i)
            {
                break;
            }
        }
        j++;
    }
    file2.close();
    *start = i;
    *end = j;
}



std::vector<std::string> MOLECULE::split(const std::string &str, const std::string &delim)
{

    std::vector<std::string> res;
	char * strs = new char[str.length() + 1]; 
	strcpy(strs, str.c_str()); 
 
	char * de = new char[delim.length() + 1];
	strcpy(de, delim.c_str());
 
	char *token = strtok(strs, de);
	while(token) {
		std::string temp = token; 
		res.push_back(temp); 
		token = strtok(NULL, de);
	}
    delete[] strs;
    delete[] de;
    return res;
}






int MOLECULE::qm(char shellname)
{
    int ret; 
    switch (tolower(shellname))
    {
        case 's':
            ret = 0;
        break;
        case 'p':
            ret = 1;
            break;
        case 'd':
            ret = 2;
            break;
        case 'f':
            ret = 3;
            break;
    }
    return ret;
}



bool MOLECULE::cmp_prim(const PRIM &a, const PRIM &b){return a.shell < b.shell;}
bool MOLECULE::cmp_ctrg(const CTRG &a, const CTRG &b){return a.shell < b.shell;}





void MOLECULE::basisload()
{
    int nshell = 0;
    for (ATOM* atom:_atoms)
    {
        int start,end;
        cut(atom->symbol,&start,&end);
        std::ifstream file(_basis_filename); 
        std::string line;
        std::regex regex_shell("^([A-Za-z]{1,2})\\s+(\\d+)\\s+1.00$");   
        std::regex regex_prim("(\\s+-?\\d+\\.\\d+[D|E|e|d]?[-|+]?\\d+){2,3}");
        std::regex sci_notation("[D|E|e|d]");
        std::smatch match_shell;
        std::smatch match_prim;
        int iline = 0;
        std::vector <CTRG> ctrgs;
        std::vector <PRIM> prims;
        std::vector<std::string> shell_info(2);
        while (getline(file,line))
        {   
            if ( (start <= iline) && (iline <= end) )
            {    
                if (std::regex_search(line,match_shell,regex_shell))
                {
                    shell_info[0] = match_shell[1].str();
                    shell_info[1] = match_shell[2].str();
                    if (shell_info[0].size()==1)// single shell
                    {
                        CTRG ctrg;
                        ctrg.shell = qm(shell_info[0][0]);
                        ctrg.nprim = stoi(shell_info[1]);
                        ctrg.nctr = 1;
                        ctrgs.push_back(ctrg);
                        nshell++;
                    }

                    if (shell_info[0].size()==2)// sp shell
                    {
                        CTRG ctrg1,ctrg2;
                        ctrg1.shell = qm(shell_info[0][0]);
                        ctrg1.nprim = stoi(shell_info[1]);
                        ctrg1.nctr = 1;
                        nshell++;
                        ctrg2.shell = qm(shell_info[0][1]);
                        ctrg2.nprim = stoi(shell_info[1]);
                        ctrg2.nctr = 1;
                        nshell++;
                        ctrgs.push_back(ctrg1);
                        ctrgs.push_back(ctrg2);
                    }
                }

                if (std::regex_search(line,match_prim,regex_prim))
                {
                    std::vector<std::string> res = split(match_prim[0].str()," ");

                    if (res.size() == 2) // single shell
                    {
                        PRIM prim;
                        prim.symbol = shell_info[0][0];
                        prim.shell = qm(shell_info[0][0]);
                        prim.exp = stod(std::regex_replace(res[0],sci_notation,"e"));
                        prim.coeff = stod(std::regex_replace(res[1],sci_notation,"e"));
                        prims.push_back(prim);

                    }
                    if (res.size() == 3) //sp shell 
                    {
                        PRIM prim1,prim2;
                        prim1.symbol = shell_info[0][0];
                        prim1.shell = qm(shell_info[0][0]);
                        prim1.exp = stod(std::regex_replace(res[0],sci_notation,"e"));
                        prim1.coeff = stod(std::regex_replace(res[1],sci_notation,"e"));
                        prim2.symbol = shell_info[0][1];
                        prim2.shell = qm(shell_info[0][1]);
                        prim2.exp = stod(std::regex_replace(res[0],sci_notation,"e"));
                        prim2.coeff = stod(std::regex_replace(res[2],sci_notation,"e"));
                        prims.push_back(prim1);
                        prims.push_back(prim2);
                    }
                }
            }   
            iline++;
        }
        file.close();
        std::sort(prims.begin(),prims.end(),cmp_prim);
        std::sort(ctrgs.begin(),ctrgs.end(),cmp_ctrg);
        atom->prims = prims;
        atom->ctrgs = ctrgs;
        _nshell = nshell;

    }
}



void MOLECULE::read()
{
    std::ifstream file(_input); 
    if(!file){std::cout <<"open input file error!"<< std::endl;exit(0);}
    std::string line;
    std::regex charge_multipy("^(\\d)\\s(\\d)$");
    std::regex coordinate("([a-zA-z]{1,2})\\s+(-?\\d+\\.\\d+)\\s+(-?\\d+\\.\\d+)\\s+(-?\\d+\\.\\d+)");
    std::regex basisset("basis=(.+)");
    std::smatch match_charge_multipy;
    std::smatch match_coordinate;
    std::smatch match_basisset;

    while (getline(file,line))
    {
        std::cout << line << std::endl;
        if (std::regex_search(line,match_basisset,basisset))
        {
            _basis_set = match_basisset[1].str();
            _basis_filename = "../../basis/" + match_basisset[1].str()  + ".0.gbs";        
        }
        if (regex_search(line,match_charge_multipy,charge_multipy))
        {
            _charge = stoi(match_charge_multipy[1].str());
            _spin = stoi(match_charge_multipy[2].str());
        }
        if (regex_search(line,match_coordinate,coordinate))
        {
            const std::string symbol = match_coordinate[1].str();
            const int charge = getcharge(symbol);
            const double mass = getmass(symbol);
            const double x = stod(match_coordinate[2].str())/BOHR2ANG;//convert to bhor.
            const double y = stod(match_coordinate[3].str())/BOHR2ANG;
            const double z = stod(match_coordinate[4].str())/BOHR2ANG;
            ATOM *atom = new(ATOM);
            atom->symbol = symbol;
            atom->charge = charge;
            atom->mass = mass;
            atom->coor = std::vector<double> {x,y,z};
            _atoms.push_back(atom);

        }
    }
    _natm = _atoms.size();
    file.close();
    basisload();

    _atm = new int[_natm * ATM_SLOTS]();
    _bas = new int[_nshell * BAS_SLOTS]();
    _env = new double[MAX_ENV]();
    _shell_label = new int[_nshell]();

    /*

    for (const ATOM* atom:_atoms)
    {
        std::cout << atom->symbol << std::endl;
        for (const CTRG ctrg:atom->ctrgs)
        {
            printf("%d%10d%10d\n",ctrg.shell,ctrg.nprim,ctrg.nctr);
        }
        for (const PRIM prim:atom->prims)
        {
            printf("%c%10d%10.3lf%10.3lf\n",prim.symbol,prim.shell,prim.exp,prim.coeff);
        }
    }
    */



    int off = PTR_ENV_START;
    int total_charge=0;
    for (int iatm = 0; iatm < _natm; iatm++)
    {
        _atm[CHARGE_OF + ATM_SLOTS * iatm] = _atoms[iatm]->charge;//charge for atom
        _atm[PTR_COORD + ATM_SLOTS * iatm] = off;//
        _atm[NUC_MOD_OF + ATM_SLOTS * iatm] = 1;
        total_charge += _atoms[iatm]->charge; //total charge for each atom.
        _env[off + 0] = _atoms[iatm]->coor[0];//x angstrom to bhor
        _env[off + 1] = _atoms[iatm]->coor[1];//y
        _env[off + 2] = _atoms[iatm]->coor[2];//z
        off +=3;
        _atm[PTR_ZETA + ATM_SLOTS * iatm] = off;
        _env[off + 3] = 0;
        off += 1;
    }
    

    _nbf = 0;
    int shl_off=0;
    for (int iatm = 0; iatm < _natm; iatm++)
    {   
        int ioff = 0;
        int joff = 0;
        for (int ishl = shl_off; ishl < shl_off + _atoms[iatm]->ctrgs.size(); ishl++)
        {
            _bas[ATOM_OF  + BAS_SLOTS * ishl] = iatm;
            _bas[ANG_OF   + BAS_SLOTS * ishl] = _atoms[iatm]->ctrgs[ioff].shell;
            int l = _atoms[iatm]->ctrgs[ioff].shell;
            _shell_label[ishl] = _nbf;
            _nbf += 2*l+1;
            _bas[NPRIM_OF + BAS_SLOTS * ishl] = _atoms[iatm]->ctrgs[ioff].nprim;
            _bas[NCTR_OF  + BAS_SLOTS * ishl] = 1;
            _bas[PTR_EXP  + BAS_SLOTS * ishl]  = off;
            int koff = 0;
            for (int iprim = joff; iprim < joff+_atoms[iatm]->ctrgs[ioff].nprim; iprim++)
            {
                _env[off + koff] = _atoms[iatm]->prims[iprim].exp;
                koff++;
            }
  
            off += _atoms[iatm]->ctrgs[ioff].nprim;
            _bas[PTR_COEFF+ BAS_SLOTS * ishl] = off;

            int loff = 0;
            for (int iprim = joff; iprim < joff + _atoms[iatm]->ctrgs[ioff].nprim; iprim++)
            {
                double exp = _atoms[iatm]->prims[iprim].exp;
                double coeff = _atoms[iatm]->prims[iprim].coeff;
                int ang = _atoms[iatm]->ctrgs[ioff].shell;
                _env[off + loff] = coeff * CINTgto_norm(ang, exp);
                loff++;
            }
            off += _atoms[iatm]->ctrgs[ioff].nprim;
            joff += _atoms[iatm]->ctrgs[ioff].nprim;
            ioff++;      
        }
        shl_off += _atoms[iatm]->ctrgs.size();
    }
    normalize();

/*

    printf("_natm\n");
    for (int i = 0; i < _natm * ATM_SLOTS; i++)
    {
        printf("%8d",_atm[i]);
        if ((i+1)%ATM_SLOTS==0)
        {
            printf("\n");
        }
    }
    printf("\n");

    printf("_bas\n");
    for (int i = 0; i < _nshell * BAS_SLOTS; i++)
    {
        printf("%8d",_bas[i]);
        if ((i+1)%BAS_SLOTS==0)
        {
            printf("\n");
        }
    }
    printf("\n");

    printf("_env\n");
    for (int i = 0; i < 200; i++)
    {
        printf("%10.3lf",_env[i]);
        if ((i+1)%4==0)
        {
            printf("\n");
        }

    }
    printf("\n");
*/


    // calculate occ;
    _occ = new int[_nbf]();

    int nocc = (total_charge -_charge)/2;
    for (int i = 0; i < nocc; i++)
    {
        _occ[i] = 2;
    }

    //calculate enuc
    _enuc = 0.;
    for (int i = 0; i < _natm; i++)
    {
        for (int j = i+1; j < _natm; j++)

        {
            double dx = _atoms[i]->coor[0] - _atoms[j]->coor[0];
            double dy = _atoms[i]->coor[1] - _atoms[j]->coor[1];
            double dz = _atoms[i]->coor[2] - _atoms[j]->coor[2];
            double dist = sqrt (dx * dx + dy * dy + dz * dz);
            _enuc += (double)_atoms[i]->charge * (double)_atoms[j]->charge / dist;  
        }
    }

    //printf("\n");
    //printf("TOTAL CHARGE:%d\n",total_charge);
    //printf("OCC= ");
    //for (int i = 0; i < _nbf; i++)
    //{
    //    printf("%d ",_occ[i]);
    //}
    //printf("\n");
    //printf("ENUC:%15.10lf a.u\n",_enuc);


    //printf("%d\n",_nbf);
    //for (int i = 0; i < _nshell; i++)
    //{
    //    std::cout << _shell_label[i] << " ";
    //}
    //std::cout << std::endl;

    /*
    printf("atm\n");
    for (int i = 0; i < _natm * ATM_SLOTS; i++)
    {
          printf("%8d",_atm[i]);
          if ((i+1)%ATM_SLOTS == 0)
          {
                printf("\n");
          }
    }
    printf("\n");
    printf("bas\n");
    for (int i = 0; i < _nshell * BAS_SLOTS; i++)
    {
          printf("%8d",_bas[i]);
          if ((i+1)%BAS_SLOTS == 0)
          {
                printf("\n");
          }
    }
    printf("\n");
    printf("env\n");
    for (int i = 20; i < 73; i++)
    {
          printf("%15.3lf",_env[i]);
          if ((i+1)%4 == 0)
          {
                printf("\n");
          }
    }
    printf("\n");
    */
}

int MOLECULE::packed_index(int i ,int j)
{
      int min = MIN0(i, j);
      int max = MAX0(i, j);
      return (min + max * (max + 1) / 2);
}


void MOLECULE::normalize()
{
    int shls[2]={0},di,dj;
    for (int ishl = 0; ishl < _nshell; ishl++)
    {
        shls[0] = ishl;
        shls[1] = ishl;
        di = CINTcgto_spheric(ishl,_bas);
        dj = CINTcgto_spheric(ishl,_bas);
        double *buf = new double[di * dj]();
        cint1e_ovlp_sph(buf,shls,_atm,_natm,_bas,_nshell,_env);
        for (int j = 0; j < _bas[NPRIM_OF + BAS_SLOTS * ishl]; j++)
        {
            _env[_bas[PTR_COEFF + BAS_SLOTS * ishl] + j] /= sqrt(buf[0]);
        }
        delete buf;
    }
}






void MOLECULE::one_electron_integrals()
{
   /*
      para:
      atm:  the information of atom in system.
      bas:  the records of basis in system.
      env:  the information of basis in system.
      natm: the number of atom in system.
      nshl: the number of shell for nolecule.
      nbf: the numbr of orbital form basis function.
    */

    int shls[2]={0},di,dj;
    _S = new double[_nbf * _nbf]();
    _T = new double[_nbf * _nbf]();
    _V = new double[_nbf * _nbf]();
    for (int i = 0; i < _nshell; i++)
    {
        shls[0] = i;
        di = CINTcgto_spheric(i, _bas);
        for (int j = i; j < _nshell; j++)
        {
            shls[1] = j;
            dj = CINTcgto_spheric(j, _bas);
            double *bufovp = new double[di * dj]();
            double *bufkin = new double[di * dj]();
            double *bufnuc = new double[di * dj]();
            int retovp = cint1e_ovlp_sph(bufovp,shls,_atm,_natm,_bas,_nshell,_env);
            int retkin = cint1e_kin_sph(bufkin,shls,_atm,_natm,_bas,_nshell,_env);
            int retnuc = cint1e_nuc_sph(bufnuc,shls,_atm,_natm,_bas,_nshell,_env);
            if (retovp == 0)
                  continue;
            for (int i2 = _shell_label[j], k = 0; k < dj; i2++, k++) 
            {
                  for (int i1 = _shell_label[i],l = 0; l < di; i1++, l++) 
                  {
                        _S[i1 * _nbf + i2] = bufovp[di * k + l];
                        _S[i2 * _nbf + i1] = bufovp[di * k + l];
                  }
            }

            if (retkin == 0)
                  continue;
            for (int i2 = _shell_label[j], k = 0; k < dj; i2++, k++) 
            {
                  for (int i1 = _shell_label[i],l = 0; l < di; i1++, l++) 
                  {
                        _T[i1 * _nbf + i2] = bufkin[di * k + l];
                        _T[i2 * _nbf + i1] = bufkin[di * k + l];
                  }
            }

            if (retnuc == 0)
                  continue;
            for (int i2 = _shell_label[j], k = 0; k < dj; i2++, k++) 
            {
                  for (int i1 = _shell_label[i],l = 0; l < di; i1++, l++) 
                  {
                        _V[i1 * _nbf + i2] = bufnuc[di * k + l];
                        _V[i2 * _nbf + i1] = bufnuc[di * k + l];
                  }
            }
            delete[] bufovp;
            delete[] bufkin;
            delete[] bufnuc;
        }
    
    }


}





void MOLECULE::two_electron_integrals()
{
    int shls[4]={0},di,dj,dk,dl;
    unsigned long long int nij = _nbf * (_nbf + 1) / 2; //important
    unsigned long long int nijkl = nij * (nij + 1) / 2; //important
    _EE = new double[nijkl]();
    CINTOpt *opt = NULL;
    cint2e_sph_optimizer (&opt,_atm,_natm,_bas,_nshell,_env);
    for (int i = 0; i < _nshell; i++)
    {
        shls[0] = i;
        di = CINTcgto_spheric(i,_bas);
        for (int j = 0; j <= i; j++)
        {
            shls[1] = j;
            dj = CINTcgto_spheric(j,_bas);
            for (int k = 0; k < _nshell; k++)
            {
                shls[2] = k;
                dk = CINTcgto_spheric(k,_bas);
                for (int l = 0; l <= k; l++)
                {
                    shls[3] = l;
                    dl = CINTcgto_spheric(l,_bas);
                    int imax = _shell_label[i] + di;
                    int jmax = _shell_label[j] + dj;
                    int ijmax = jmax + imax * (imax + 1) / 2;
                    int kmin = _shell_label[k];
                    int lmin = _shell_label[l];
                    int klmin = lmin + kmin * (kmin + 1) / 2;
                    if (ijmax < klmin)
                        continue;
                    double *buferi = new double[di*dj*dk*dl]();
                    int ret2e = cint2e_sph(buferi,shls,_atm,_natm,_bas,_nshell,_env,opt);
                    if (ret2e == 0)
                        continue;
                    for (int i4 = _shell_label[l], j4 = 0; j4 < dl; i4++, j4++)
                    {
                        for (int i3 = _shell_label[k], j3 = 0; j3 < dk; i3++, j3++)
                        {
                            if (i3 < i4)
                                continue;
                            int kl = packed_index (i3, i4);
                            for (int i2 = _shell_label[j], j2 = 0; j2 < dj; i2++, j2++)
                            {
                                for (int i1 = _shell_label[i], j1 = 0; j1 < di; i1++, j1++)
                                {
                                    if (i1 < i2)
                                        continue;
                                    int ij = packed_index (i1, i2);
                                    if (ij < kl)
                                        continue;
                                    int ijkl = packed_index (ij, kl);
                                    _EE[ijkl] = buferi[j1 + j2 * di + j3 * di * dj+ j4 * di * dj * dk];
                                }
                            }
                        }

                    }
                    delete[] buferi;
                }
            }
        }
    }
    CINTdel_optimizer (&opt);
}




void MOLECULE::bulid()
{
    printf("\n");
    printf("call libcint integral libary...\n");
    clock_t integral_start = clock();
    one_electron_integrals();
    two_electron_integrals();
    clock_t integral_end = clock();
    printf("libcint calculation time = %.3lf SEC.\n",(double)(integral_end-integral_start)/CLOCKS_PER_SEC);


/*
  
    printf("S\n");
    for (int i = 0; i < _nbf * _nbf; i++)
    {
          printf("%8.3lf",_S[i]);
          if ((i+1)%_nbf == 0)
          {
                printf("\n");
          }
    }
    printf("\n");

    printf("T\n");
    for (int i = 0; i < _nbf * _nbf; i++)
    {
          printf("%8.3lf",_T[i]);
          if ((i+1)%_nbf == 0)
          {
                printf("\n");
          }
    }
    printf("\n");


    printf("V\n");
    for (int i = 0; i < _nbf * _nbf; i++)
    {
          printf("%8.3lf",_V[i]);
          if ((i+1)%_nbf == 0)
          {
                printf("\n");
          }
    }
    printf("\n");
    */  


}




