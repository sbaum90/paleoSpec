from __future__ import division
import numpy as np

# define target list
mmol=153.31 # g/mol
Th_length=302. # [Aa]
# target masses in GeV and number of nucleons
# Hydrogen
m_H=0.9389
A_H=1.
Z_H=1.
n_H=1.
# Oxygen
m_O=14.903
A_O=(99.76*16+.04*17+.2*18)/100
Z_O=8.
n_O=11.
# Fluorine
m_F=17.6969
A_F=19.
Z_F=9.
n_F=1.
# Magnesium
m_Mg=22.64
A_Mg=(24.*79.0+25.*10.0+26.*11.0)/100.
Z_Mg=12.
n_Mg=3.
# Aluminium
m_Al=25.13314
A_Al=27.
Z_Al=13.
n_Al=1.
# Silicon
m_Si=26.161
A_Si=(92.2*28+4.7*29+3.1*30)/100
Z_Si=14.
n_Si=3.
# Potassium
m_K=36.4198
A_K=(39.*93.258+40.*0.012+41.*6.730)/100.
Z_K=19.
n_K=1.


# lists
mT_list=[m_H,m_O,m_F,m_Mg,m_Al,m_Si,m_K]
AT_list=[A_H,A_O,A_F,A_Mg,A_Al,A_Si,A_K]
ZT_list=[Z_H,Z_O,Z_F,Z_Mg,Z_Al,Z_Si,Z_K]
nameT_list=['H','O','F','Mg','Al','Si','K']
massFrac_list=np.array([m_H*n_H,m_O*n_O,m_F*n_F,m_Mg*n_Mg,m_Al*n_Al,m_Si*n_Si,m_K*n_K])/(m_H*n_H+m_O*n_O+m_F*n_F+m_Mg*n_Mg+m_Al*n_Al+m_Si*n_Si+m_K*n_K)
