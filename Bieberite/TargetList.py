from __future__ import division
import numpy as np

# define target list
mmol=281.10 # g/mol
Th_length=440. # [Aa]
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
# Sulfur
m_S=29.86
A_S=(94.99*32+.75*33+4.25*34+.01*36)/100
Z_S=16.
n_S=1.
#Cobalt
m_Co=54.89592
A_Co=59.
Z_Co=27.
n_Co=1.

# lists
mT_list=[m_H,m_O,m_S,m_Co]
AT_list=[A_H,A_O,A_S,A_Co]
ZT_list=[Z_H,Z_O,Z_S,Z_Co]
nameT_list=['H','O','S','Co']
massFrac_list=np.array([m_H*n_H,m_O*n_O,m_S*n_S,m_Co*n_Co])/(m_H*n_H+m_O*n_O+m_S*n_S+m_Co*n_Co)
