

from pyscf import gto
import numpy as np
from pyscf import scf
np.set_printoptions(threshold=10000,linewidth=1000,precision=3,suppress=True)

mol = gto.Mole()
mol.atom = """
N 0.00000000,  0.00000000,  0.00000000;
N 0.00000000,  0.00000000,  1.10000000;
"""
mol.basis = '6-31g'
mol.basis = 'def2-svp'
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

#converged SCF energy = -108.851780982863