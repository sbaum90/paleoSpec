from __future__ import division
import numpy as np

# define target list
mmol=659.66 # g/mol
Th_length=449. # [Aa]
# target masses in GeV and number of nucleons
#Hydrogen
m_H=0.9389
A_H=1.
Z_H=1.
n_H=2*22.05
# Oxygen
m_O=14.903
A_O=(99.76*16+.04*17+.2*18)/100
Z_O=8.
n_O=22.05+4*2.01
# Magnesium
m_Mg=22.64
A_Mg=(24.*79.0+25.*10.0+26.*11.0)/100.
Z_Mg=12.
n_Mg=2.92
# Phosphorus
m_P=28.85188
A_P=31.
Z_P=15.
n_P=2.01
#Iron
m_Fe=52.019
A_Fe=(54.*5.85+56.*91.75+57.*2.12+58*0.28)/100.
Z_Fe=26.
n_Fe=0.01

# lists
mT_list=[m_H,m_O,m_Mg,m_P,m_Fe]
AT_list=[A_H,A_O,A_Mg,A_P,A_Fe]
ZT_list=[Z_H,Z_O,Z_Mg,Z_P,Z_Fe]
nameT_list=['H','O','Mg','P','Fe']
massFrac_list=np.array([m_H*n_H,m_O*n_O,m_Mg*n_Mg,m_P*n_P,m_Fe*n_Fe])/(m_H*n_H+m_O*n_O+m_Mg*n_Mg+m_P*n_P+m_Fe*n_Fe)
