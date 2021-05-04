from __future__ import division
import numpy as np

# define target list
mmol=237.69 # g/mol
Th_length=460. # [Aa]
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
# Chlorine
m_Cl=33.02
A_Cl=(76*35+24*37)/100
Z_Cl=17.
n_Cl=2.
# Nickel
m_Ni=54.6726
A_Ni=(58.*68.077+60.*26.223+61.*1.140+62.*3.635+64.*0.926)/100.   
Z_Ni=28.
n_Ni=1.

# lists
mT_list=[m_H,m_O,m_Cl,m_Ni]
AT_list=[A_H,A_O,A_Cl,A_Ni]
ZT_list=[Z_H,Z_O,Z_Cl,Z_Ni]
nameT_list=['H','O','Cl','Ni']
massFrac_list=np.array([m_H*n_H,m_O*n_O,m_Cl*n_Cl,m_Ni*n_Ni])/(m_H*n_H+m_O*n_O+m_Cl*n_Cl+m_Ni*n_Ni)
