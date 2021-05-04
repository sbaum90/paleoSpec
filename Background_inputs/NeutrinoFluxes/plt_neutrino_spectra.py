import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import MaxNLocator

import numpy as np

# fix some plotting settings
fs = 14
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'legend.fontsize': fs,
       'axes.labelsize': fs,
       'axes.titlesize': fs,
       'xtick.labelsize': fs,
       'ytick.labelsize': fs})
plt.rcParams["figure.figsize"] = (7, 6)

plt_colors = ['k',
   '#1b9e77',
   '#d95f02',
   '#7570b3',
   '#e7298a',
   '#66a61e',
   '#e6ab02']

# import spectra
solar_pp = np.loadtxt('CAJOH_pp.dat')
solar_hep = np.loadtxt('CAJOH_hep.dat')
solar_8B = np.loadtxt('CAJOH_8B.dat')
solar_13N = np.loadtxt('CAJOH_13N.dat')
solar_15O = np.loadtxt('CAJOH_15O.dat')
solar_17F = np.loadtxt('CAJOH_17F.dat')
solar_lines = np.loadtxt('CAJOH_lines.dat')

DSNB = np.loadtxt('DSNB_spec.dat')
GSNB = np.loadtxt('Galactic_SN_spec.dat')

atm = np.loadtxt('CAJOH_Atm.dat')

# function for point spectra
def plt_lines(E, Flux, color, linestyle, label):
   Evec = np.logspace(-1, 3, 4001)
   Fvec = np.zeros(Evec.shape)
   Eind = np.argmin(np.abs(Evec-E))
   Evec[Eind] = E
   Fvec[Eind] = Flux
   plt.plot(Evec[:Eind+1], Fvec[:Eind+1], c=color, linestyle=linestyle, label=label)


plt.close('all')
plt.plot(solar_pp[:,0], solar_pp[:,1], c=plt_colors[1], linestyle='-', label=r'pp')
plt_lines(solar_lines[2,0], solar_lines[2,1], plt_colors[1], '--', r'pep')
plt.plot(solar_hep[:,0], solar_hep[:,1], c=plt_colors[1], linestyle='-.', label=r'hep')
plt_lines(solar_lines[0,0], solar_lines[0,1], plt_colors[2], '-', r'$^7$Be [384\,keV]')
plt_lines(solar_lines[1,0], solar_lines[1,1], plt_colors[2], '--', r'$^7$Be [861\,keV]')
plt.plot(solar_8B[:,0], solar_8B[:,1], c=plt_colors[2], linestyle='-.', label=r'$^{8}$B')
plt.plot(solar_13N[:,0], solar_13N[:,1], c=plt_colors[3], linestyle='-', label=r'$^{13}$N')
plt.plot(solar_15O[:,0], solar_15O[:,1], c=plt_colors[3], linestyle='--', label=r'$^{15}$O')
plt.plot(solar_17F[:,0], solar_17F[:,1], c=plt_colors[3], linestyle='-.', label=r'$^{17}$F')

plt.plot(DSNB[:,0], DSNB[:,1], c=plt_colors[4], linestyle='-', label=r'DSNB')
plt.plot(GSNB[:,0], GSNB[:,1], c=plt_colors[4], linestyle='--', label=r'GSNB')

plt.plot(atm[:,0], atm[:,1], c=plt_colors[5], linestyle='-', label=r'atm')

plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'Neutrino Energy [MeV]')
plt.ylabel(r'Neutrino Flux [cm$^{-2}$s$^{-1}$MeV$^{-1}$]')
plt.xlim(0.1, 1e3)
plt.ylim(1e-4, 1e12)
plt.legend()
plt.tick_params(right=True,top=True)
plt.tick_params(which='minor',right=True,top=True)
plt.tight_layout()
plt.savefig('Neutrino_Fluxes.pdf')