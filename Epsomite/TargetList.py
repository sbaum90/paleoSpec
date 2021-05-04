from __future__ import division
import numpy as np

# define target list
mmol=246.48 # g/mol
Th_length=457. # [Aa]
# target masses in GeV and number of nucleons
#Hydrogen
m_H=0.9389
A_H=1.
Z_H=1.
n_H=14.
# Oxygen
m_O=14.903
A_O=(99.76*16+.04*17+.2*18)/100
Z_O=8.
n_O=11.
#Magnesium
m_Mg=22.64
A_Mg=(24.*79.+25.*10.+26.*11.)/100.
Z_Mg=12.
n_Mg=1.
# Sulfur
m_S=29.86
A_S=(94.99*32+.75*33+4.25*34+.01*36)/100
Z_S=16.
n_S=1.

# lists
mT_list=[m_H,m_O,m_Mg,m_S]
AT_list=[A_H,A_O,A_Mg,A_S]
ZT_list=[Z_H,Z_O,Z_Mg,Z_S]
nameT_list=['H','O','Mg','S']
massFrac_list=np.array([m_H*n_H,m_O*n_O,m_Mg*n_Mg,m_S*n_S])/(m_H*n_H+m_O*n_O+m_Mg*n_Mg+m_S*n_S)
