from __future__ import division
import numpy as np

# define target list
mmol=136.14 # g/mol
Th_length=304. # [Aa]
# target masses in GeV and number of nucleons
# Oxygen
m_O=14.903
A_O=(99.76*16+.04*17+.2*18)/100
Z_O=8.
n_O=4.
# Sulfur
m_S=29.86
A_S=(94.99*32+.75*33+4.25*34+.01*36)/100
Z_S=16.
n_S=1.
#Calcium
m_Ca=37.332
A_Ca=(40.*96.941+42.*0.647+43.*0.135+44.*2.086+48.*0.187)/100.
Z_Ca=20.
n_Ca=1.

# lists
mT_list=[m_O,m_S,m_Ca]
AT_list=[A_O,A_S,A_Ca]
ZT_list=[Z_O,Z_S,Z_Ca]
nameT_list=['O','S','Ca']
massFrac_list=np.array([m_O*n_O,m_S*n_S,m_Ca*n_Ca])/(m_O*n_O+m_S*n_S+m_Ca*n_Ca)
