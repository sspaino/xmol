

from pyscf import gto
import numpy as np
from pyscf import scf
np.set_printoptions(threshold=10000,linewidth=1000,precision=3,suppress=True)

mol = gto.Mole()
mol.atom = """
C   1.21773989,   -0.70306245,   0.00000000;   
H   2.17299146,   -1.25457720,   0.00000000;   
C   1.21773989,    0.70306245,   0.00000000;   
H   2.17299146,    1.25457720,   0.00000000;   
C   0.00000000,    1.40612490,   0.00000000;   
H   0.00000000,    2.50915441,   0.00000000;   
C  -1.21773989,    0.70306245,   0.00000000;   
H  -2.17299146,    1.25457720,   0.00000000;   
C  -1.21773989,   -0.70306245,   0.00000000;   
H  -2.17299146,   -1.25457720,   0.00000000;   
C   0.00000000,   -1.40612490,   0.00000000;   
H   0.00000000,   -2.50915441,   0.00000000;
"""
mol.basis = '6-31g'
#mol.basis = 'def2-svp'
#mol.basis = 'aug-cc-pvtz'
mol.build()
rhf_h2o = scf.RHF(mol)
rhf_h2o.kernel()
#S = mol.intor("int1e_ovlp")
#V = mol.intor("int1e_nuc")
#T = mol.intor("int1e_kin")
#eri = mol.intor("int2e")
#print(mol._atm)
#print(mol._bas)
#print(mol.atom_coords())
#print(mol._env)

#print("Overlap matrix")
#print(S)
#print(S.shape)
#lam_s, l_s = np.linalg.eigh(S)
#lam_s = lam_s * np.eye(len(lam_s))
#lam_sqrt_inv = np.sqrt(np.linalg.inv(lam_s))
#symm_orthog = np.dot(l_s, np.dot(lam_sqrt_inv, l_s.T))
#print("Symmetric orthogonalization matrix")
#print(symm_orthog)
#print(S.shape)
#print(mol._basis)
#print(len(mol._env))
#for key,value in mol._basis.items():
#        print(key)
#        for i in value:
#                for j in i:
#                        print(j)
#

#converged SCF energy = -230.618848335009