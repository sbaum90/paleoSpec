from __future__ import division
import numpy as np

# define target list
mmol=12.01 # g/mol
Th_length=196. # [Aa]
# target masses in GeV and number of nucleons
#Carbon
m_C=11.188
A_C=(12.*98.9+13.*1.1)/100.
Z_C=6.
n_C=1.
# lists
mT_list=[m_C]
AT_list=[A_C]
ZT_list=[Z_C]
nameT_list=['C']
massFrac_list=np.array([m_C*n_C])/(m_C*n_C)
