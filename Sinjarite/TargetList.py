from __future__ import division
import numpy as np

# define target list
mmol=174.01 # g/mol
Th_length=507. # [Aa]
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
n_O=2.
# Chlorine
m_Cl=33.02
A_Cl=(76*35+24*37)/100
Z_Cl=17.
n_Cl=2.
#Calcium
m_Ca=37.332
A_Ca=(40.*96.941+42.*0.647+43.*0.135+44.*2.086+48.*0.187)/100.
Z_Ca=20.
n_Ca=1.

# lists
mT_list=[m_H,m_O,m_Cl,m_Ca]
AT_list=[A_H,A_O,A_Cl,A_Ca]
ZT_list=[Z_H,Z_O,Z_Cl,Z_Ca]
nameT_list=['H','O','Cl','Ca']
massFrac_list=np.array([m_H*n_H,m_O*n_O,m_Cl*n_Cl,m_Ca*n_Ca])/(m_H*n_H+m_O*n_O+m_Cl*n_Cl+m_Ca*n_Ca)
