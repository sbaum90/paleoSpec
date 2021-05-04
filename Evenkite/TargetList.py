from __future__ import division
import numpy as np

# define target list
mmol=338.66 # g/mol
Th_length=687. # [Aa]
# target masses in GeV and number of nucleons
#Hydrogen
m_H=0.9389
A_H=1.
Z_H=1.
n_H=6.+44.
# Carbon
m_C=11.188
A_C=(12.*98.9+13.*1.1)/100.
Z_C=6.
n_C=2.+22.

# lists
mT_list=[m_H,m_C]
AT_list=[A_H,A_C]
ZT_list=[Z_H,Z_C]
nameT_list=['H','C']
massFrac_list=np.array([m_H*n_H,m_C*n_C])/(m_H*n_H+m_C*n_C)
