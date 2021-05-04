import numpy as np


Flux_pp = 5.98e10  
Flux_hep = 7.98e3 
Flux_8B = 5.16e6 
Flux_13N = 2.78e8 
Flux_15O = 2.05e8  
Flux_17F = 5.29e6  
Flux_Atm = 10.54e0  

solar_pp_shape = np.loadtxt('CAJOH_pp_shape.dat')
solar_hep_shape = np.loadtxt('CAJOH_hep_shape.dat')
solar_8B_shape = np.loadtxt('CAJOH_8B_shape.dat')
solar_13N_shape = np.loadtxt('CAJOH_13N_shape.dat')
solar_15O_shape = np.loadtxt('CAJOH_15O_shape.dat')
solar_17F_shape = np.loadtxt('CAJOH_17F_shape.dat')
atm_shape = np.loadtxt('CAJOH_Atm_shape.dat')

fo = open('CAJOH_pp.dat','w')
fo.write('# (E [MeV])  (Flux [1/MeV/cm^2/s]\n')
for i in range(len(solar_pp_shape[:,0])):
   fo.write('{:3E}  '.format(solar_pp_shape[i,0]))
   fo.write('{:3E}\n'.format(Flux_pp*solar_pp_shape[i,1]))

fo.close()

fo = open('CAJOH_hep.dat','w')
fo.write('# (E [MeV])  (Flux [1/MeV/cm^2/s]\n')
for i in range(len(solar_hep_shape[:,0])):
   fo.write('{:3E}  '.format(solar_hep_shape[i,0]))
   fo.write('{:3E}\n'.format(Flux_hep*solar_hep_shape[i,1]))

fo.close()

fo = open('CAJOH_8B.dat','w')
fo.write('# (E [MeV])  (Flux [1/MeV/cm^2/s]\n')
for i in range(len(solar_8B_shape[:,0])):
   fo.write('{:3E}  '.format(solar_8B_shape[i,0]))
   fo.write('{:3E}\n'.format(Flux_8B*solar_8B_shape[i,1]))

fo.close()

fo = open('CAJOH_13N.dat','w')
fo.write('# (E [MeV])  (Flux [1/MeV/cm^2/s]\n')
for i in range(len(solar_13N_shape[:,0])):
   fo.write('{:3E}  '.format(solar_13N_shape[i,0]))
   fo.write('{:3E}\n'.format(Flux_13N*solar_13N_shape[i,1]))

fo.close()

fo = open('CAJOH_15O.dat','w')
fo.write('# (E [MeV])  (Flux [1/MeV/cm^2/s]\n')
for i in range(len(solar_15O_shape[:,0])):
   fo.write('{:3E}  '.format(solar_15O_shape[i,0]))
   fo.write('{:3E}\n'.format(Flux_15O*solar_15O_shape[i,1]))

fo.close()

fo = open('CAJOH_17F.dat','w')
fo.write('# (E [MeV])  (Flux [1/MeV/cm^2/s]\n')
for i in range(len(solar_17F_shape[:,0])):
   fo.write('{:3E}  '.format(solar_17F_shape[i,0]))
   fo.write('{:3E}\n'.format(Flux_17F*solar_17F_shape[i,1]))

fo.close()

fo = open('CAJOH_Atm.dat','w')
fo.write('# (E [MeV])  (Flux [1/MeV/cm^2/s]\n')
for i in range(len(atm_shape[:,0])):
   fo.write('{:3E}  '.format(atm_shape[i,0]))
   fo.write('{:3E}\n'.format(Flux_Atm*atm_shape[i,1]))

fo.close()


