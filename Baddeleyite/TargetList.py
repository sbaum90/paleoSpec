from __future__ import division
import numpy as np

# define target list
mmol=123.22 # g/mol
Th_length=219. # [Aa]
# target masses in GeV and number of nucleons
# Oxygen
m_O=14.903
A_O=(99.76*16+.04*17+.2*18)/100
Z_O=8.
n_O=2.
# Zirkonium
m_Zr=84.975
A_Zr=(90.*51.45+91.*11.22+92.*17.15+94.*17.38+96.*2.80)/100
Z_Zr=40.
n_Zr=1.


# lists
mT_list=[m_O,m_Zr]
AT_list=[A_O,A_Zr]
ZT_list=[Z_O,Z_Zr]
nameT_list=['O','Zr']
massFrac_list=np.array([m_O*n_O,m_Zr*n_Zr])/(m_O*n_O+m_Zr*n_Zr)
