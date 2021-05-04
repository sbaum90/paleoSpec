from __future__ import division
import numpy as np

# define target list
mmol=153.31 # g/mol
Th_length=261. # [Aa]
# target masses in GeV and number of nucleons
# Oxygen
m_O=14.903
A_O=(99.76*16+.04*17+.2*18)/100
Z_O=8.
n_O=4.
# Magnesium
m_Mg=22.64
A_Mg=(24.*79.0+25.*10.0+26.*11.0)/100.
Z_Mg=12.
n_Mg=1.6
# Silicon
m_Si=26.161
A_Si=(92.2*28+4.7*29+3.1*30)/100
Z_Si=14.
n_Si=1.
#Iron
m_Fe=52.019
A_Fe=(54.*5.85+56.*91.75+57.*2.12+58*0.28)/100.
Z_Fe=26.
n_Fe=0.4

# lists
mT_list=[m_O,m_Mg,m_Si,m_Fe]
AT_list=[A_O,A_Mg,A_Si,A_Fe]
ZT_list=[Z_O,Z_Mg,Z_Si,Z_Fe]
nameT_list=['O','Mg','Si','Fe']
massFrac_list=np.array([m_O*n_O,m_Mg*n_Mg,m_Si*n_Si,m_Fe*n_Fe])/(m_O*n_O+m_Mg*n_Mg+m_Si*n_Si+m_Fe*n_Fe)
