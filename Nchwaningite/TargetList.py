from __future__ import division
import numpy as np

# define target list
mmol=237.99 # g/mol
Th_length=298. # [Aa]
# target masses in GeV and number of nucleons
#Hydrogen
m_H=0.9389
A_H=1.
Z_H=1.
n_H=4.
# Oxygen
m_O=14.903
A_O=(99.76*16+.04*17+.2*18)/100
Z_O=8.
n_O=6.
# Silicon
m_Si=26.161
A_Si=(92.2*28+4.7*29+3.1*30)/100
Z_Si=14.
n_Si=1.
# Manganese
m_Mn=51.17446
A_Mn=55.
Z_Mn=25.
n_Mn=2.

# lists
mT_list=[m_H,m_O,m_Si,m_Mn]
AT_list=[A_H,A_O,A_Si,A_Mn]
ZT_list=[Z_H,Z_O,Z_Si,Z_Mn]
nameT_list=['H','O','Si','Mn']
massFrac_list=np.array([m_H*n_H,m_O*n_O,m_Si*n_Si,m_Mn*n_Mn])/(m_H*n_H+m_O*n_O+m_Si*n_Si+m_Mn*n_Mn)