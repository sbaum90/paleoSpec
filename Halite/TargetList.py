from __future__ import division
import numpy as np

# define target list
mmol=58.44 # g/mol
Th_length=453. # [Aa]
# Sodium
m_Na=21.41483
A_Na=23.
Z_Na=11.
n_Na=1.
# Chlorine
m_Cl=33.02
A_Cl=(76*35+24*37)/100
Z_Cl=17.
n_Cl=1.
# lists
mT_list=[m_Na,m_Cl]
AT_list=[A_Na,A_Cl]
ZT_list=[Z_Na,Z_Cl]
nameT_list=['Na','Cl']
massFrac_list=np.array([m_Na*n_Na,m_Cl*n_Cl])/(m_Na*n_Na+m_Cl*n_Cl)
