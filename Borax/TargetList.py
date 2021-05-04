from __future__ import division
import numpy as np

# define target list
mmol=381.37 # g/mol
Th_length=438. # [Aa]
# target masses in GeV and number of nucleons
# Hydrogen
m_H=0.9389
A_H=1.
Z_H=1.
n_H=20.
# Boron
m_B=10.07
A_B=(10.*20.+11.*80.)/100.
Z_B=5.
n_B=4.
# Oxygen
m_O=14.903
A_O=(99.76*16+.04*17+.2*18)/100
Z_O=8.
n_O=17.
#Sodium
m_Na=21.41483
A_Na=23.
Z_Na=11.
n_Na=2.
# lists
mT_list=[m_H,m_B,m_O,m_Na]
AT_list=[A_H,A_B,A_O,A_Na]
ZT_list=[Z_H,Z_B,Z_O,Z_Na]
nameT_list=['H','B','O','Na']
massFrac_list=np.array([m_H*n_H,m_B*n_B,m_O*n_O,m_Na*n_Na])/(m_H*n_H+m_B*n_B+m_O*n_O+m_Na*n_Na)
