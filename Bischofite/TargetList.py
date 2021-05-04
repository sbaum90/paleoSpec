from __future__ import division
import numpy as np

# define target list
mmol=203.30 # g/mol
Th_length=520. # [Aa]
# target masses in GeV and number of nucleons
#Hydrogen
m_H=0.9389
A_H=1.
Z_H=1.
n_H=12.
# Oxygen
m_O=14.903
A_O=(99.76*16+.04*17+.2*18)/100
Z_O=8.
n_O=6.
# Magnesium
m_Mg=22.64
A_Mg=(24.*79.+25.*10.+26.*11.)/100.   
Z_Mg=12.
n_Mg=1.
# Chlorine
m_Cl=33.02
A_Cl=(76*35+24*37)/100
Z_Cl=17.
n_Cl=2.


# lists
mT_list=[m_H,m_O,m_Mg,m_Cl]
AT_list=[A_H,A_O,A_Mg,A_Cl]
ZT_list=[Z_H,Z_O,Z_Mg,Z_Cl]
nameT_list=['H','O','Mg','Cl']
massFrac_list=np.array([m_H*n_H,m_O*n_O,m_Mg*n_Mg,m_Cl*n_Cl])/(m_H*n_H+m_O*n_O+m_Mg*n_Mg+m_Cl*n_Cl)
