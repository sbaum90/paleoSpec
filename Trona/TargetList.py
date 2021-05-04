from __future__ import division
import numpy as np

# define target list
mmol=226.03 # g/mol
Th_length=367. # [Aa]
# target masses in GeV and number of nucleons
#Hydrogen
m_H=0.9389
A_H=1.
Z_H=1.
n_H=5.
#Carbon
m_C=11.188
A_C=(12.*98.9+13.*1.1)/100.
Z_C=6.
n_C=2.
# Oxygen
m_O=14.903
A_O=(99.76*16+.04*17+.2*18)/100
Z_O=8.
n_O=8.
# Sodium
m_Na=21.41483
A_Na=23.
Z_Na=11.
n_Na=3.
# lists
mT_list=[m_H,m_C,m_O,m_Na]
AT_list=[A_H,A_C,A_O,A_Na]
ZT_list=[Z_H,Z_C,Z_O,Z_Na]
nameT_list=['H','C','O','Na']
massFrac_list=np.array([m_H*n_H,m_C*n_C,m_O*n_O,m_Na*n_Na])/(m_H*n_H,m_C*n_C,m_O*n_O,m_Na*n_Na)
