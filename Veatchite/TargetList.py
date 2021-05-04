from __future__ import division
import numpy as np

# define target list
mmol=653.19 # g/mol
Th_length=341. # [Aa]
# target masses in GeV and number of nucleons
# Hydrogen
m_H=0.9389
A_H=1.
Z_H=1.
n_H=7.
# Boron
m_B=10.07
A_B=(10.*20.+11.*80.)/100.
Z_B=5.
n_B=11.
# Oxygen
m_O=14.903
A_O=(99.76*16+.04*17+.2*18)/100
Z_O=8.
n_O=22.
# Strontium
m_Sr=81.62
A_Sr=(84.*0.56+86.*9.86+87.*7.00+88.*82.58)/100.
Z_Sr=38.
n_Sr=2.
# lists
mT_list=[m_H,m_B,m_O,m_Sr]
AT_list=[A_H,A_B,A_O,A_Sr]
ZT_list=[Z_H,Z_B,Z_O,Z_Sr]
nameT_list=['H','B','O','Sr']
massFrac_list=np.array([m_H*n_H,m_B*n_B,m_O*n_O,m_Sr*n_Sr])/(m_H*n_H+m_B*n_B+m_O*n_O+m_Sr*n_Sr)
