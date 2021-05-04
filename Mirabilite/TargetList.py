from __future__ import division
import numpy as np

# define target list
mmol=322.20 # g/mol
Th_length=505. # [Aa]
# target masses in GeV and number of nucleons
#Hydrogen
m_H=0.9389
A_H=1.
Z_H=1.
n_H=20.
# Oxygen
m_O=14.903
A_O=(99.76*16+.04*17+.2*18)/100
Z_O=8.
n_O=14.
#Sodium
m_Na=21.41483
A_Na=23.
Z_Na=11.
n_Na=2.
# Sulfur
m_S=29.86
A_S=(94.99*32+.75*33+4.25*34+.01*36)/100
Z_S=16.
n_S=1.

# lists
mT_list=[m_H,m_O,m_Na,m_S]
AT_list=[A_H,A_O,A_Na,A_S]
ZT_list=[Z_H,Z_O,Z_Na,Z_S]
nameT_list=['H','O','Na','S']
massFrac_list=np.array([m_H*n_H,m_O*n_O,m_Na*n_Na,m_S*n_S])/(m_H*n_H+m_O*n_O+m_Na*n_Na+m_S*n_S)
