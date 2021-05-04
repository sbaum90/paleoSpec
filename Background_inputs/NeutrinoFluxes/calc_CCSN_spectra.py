from __future__ import division
import numpy as np
import os
import subprocess
import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator
from scipy.interpolate import interp1d as interp1d
from scipy.integrate import simps as intsimps
import scipy.special as ss_funs

#######################################################
# functions for DSNB spectra
#######################################################

def SFR(z):
   """ returns the cosmic star formation rate in M_\odot/yr/Mpc^3 following Strolger+15"""
   A = 0.015 #M_\odot/yr/Mpc^3
   B = 1.5
   C = 5.0
   D = 6.1
   return A*(1+z)**C/( ((1+z)/B)**D + 1 )

def RCC(z):
   """ returns the core collapse SN rate in 1/yr/Mpc^3 followng Strolger+15"""
   k = 0.0091 # 1/M_\odot
   return k*SFR(z)

def neu_spec(E,Etot,Eavg,alpha):
   """ returns neutrino spectrum in 1/MeV
       pinched Fermi-Dirac following Horiuchi+15 """
   return Etot*(1+alpha)**(1+alpha)/ss_funs.gamma(1+alpha)*E**alpha/Eavg**(2+alpha)*np.exp(-(1+alpha)*E/Eavg)

def neu_spec_antie(E):
   """ parameters for anti-electron nu spectrum averaged over IMF from Horiuchi+15 """
   Etot = 4.3e52*6.2415e5 # total neutrino energy in MeV
   Eavg = 14.6 # average neutrino energy in MeV
   alpha = 3.3 # shape parameter
   return neu_spec(E,Etot,Eavg,alpha)

def neu_spec_e(E):
   """ parameters for electrun nu spectrum estimated from Fig. 3 in Horiuchi+15 """
   Etot = 6e52*6.2415e5 # total neutrino energy in MeV
   Eavg = 13.3 # average neutrino energy in MeV
   alpha = 3 # shape parameter
   return neu_spec(E,Etot,Eavg,alpha)

def neu_spec_x(E):
   """ parameters for other nu spectrum estimated from Fig. 3 in Horiuchi+15 """
   Etot = 2e52*6.2415e5 # total neutrino energy in MeV
   Eavg = 15 # average neutrino energy in MeV
   alpha = 3 # shape parameter
   return neu_spec(E,Etot,Eavg,alpha)

def DSNB_flux(E):
   """ returns the DSNB flux summed over all neutrino species in 1/cm^2/s/MeV """
   c = 2.997e5 # speed of light in km/s
   H0 = 67.32 # Hubble constant in km/s/Mpc
   Om = 0.3158 # matter energy density
   OLam = 1.-Om # DE energy density
   unit_fac = 1/(3.154e7*3.086e24**2) # unit conversion factor
   zvec=np.linspace(0.,5.,int(5e3+1))
   int_vec=(
      (neu_spec_e((1+zvec)*E)+neu_spec_antie((1+zvec)*E)+4*neu_spec_x((1+zvec)*E))
      *RCC(zvec)
      /( H0*np.sqrt(OLam + Om*(1+zvec)**3) )
      )
   return unit_fac*c*intsimps(int_vec,zvec)

#######################################################
# functions for galactic CC SN background
#######################################################

# get pdf for local SN distance
# start with MC
NMC=int(5e7)
Nruns=20
R_SN=np.linspace(0.,100.,int(1e4+1)) # [kpc]
R_SN_pdf=np.zeros(len(R_SN)-1)
# get SN positions in galactocentric cylindrical coordiantes [R,\phi,\theta]
Rd = 2.9 # disc scale in kpc
H = 95e-3 # disc heigth in kpc
rv = np.linspace(0.,100.,int(1e6+1))
pv = rv*np.exp(-rv/Rd)
for n in range(Nruns):
   Rvec=np.random.choice(rv,size=NMC,p=pv/np.sum(pv))
   zvec=np.random.exponential(scale=H,size=NMC)*np.random.choice([1,-1],size=NMC)
   phivec=2*np.pi*np.random.rand(NMC)
   # compute distance
   Rsol = 8.7 # distance of the Sun from the galactic center in kpc
   Hsol = 24e-3 # height of the Sun over mid-plane in kpc
   RSN=np.sqrt( (Rvec*np.cos(phivec)-Rsol)**2 + (Rvec*np.sin(phivec))**2 + (zvec-Hsol)**2 )
   # get the pdf
   R_SN_pdf+=np.histogram(RSN,bins=R_SN,density=True)[0]
   print('pdf run '+str(int(n+1))+' of '+str(int(Nruns))+' done')

R_SN_pdf=R_SN_pdf/Nruns
del rv,pv,Rvec,zvec,phivec,RSN

def Gal_flux(E,pdf=[(R_SN[1:]+R_SN[:-1])/2,R_SN_pdf]):
   """ returns the galactic flux summed over all neutrino species in 1/cm^2/s/MeV """
   RCC_gal = 2.3e-2 # [1/yr] Galactic CC SN rate from Li+ 1006.4613
   unit_fac = 1/(3.154e7*3.086e21**2) # unit conversion factor
   return (
      unit_fac/(4*np.pi)
      *(neu_spec_e(E)+neu_spec_antie(E)+4*neu_spec_x(E))
      *RCC_gal
      *intsimps(pdf[1]/pdf[0]**2,pdf[0])
      )

#######################################################
# compute spectra and write to disk
#######################################################

# export spectra
Evec = np.logspace(-1,2,301)
DSNB_vec=[DSNB_flux(E) for E in Evec]
fo=open('DSNB_spec.dat','w')
fo.write('# (E [MeV])  (Flux [1/MeV/cm^2/s]\n')
for i in range(len(Evec)):
   fo.write('{:3E}  {:3E}\n'.format(Evec[i],DSNB_vec[i]))

fo.close()

Gal_vec=[Gal_flux(E) for E in Evec]
fo=open('Galactic_SN_spec.dat','w')
fo.write('# (E [MeV])  (Flux [1/MeV/cm^2/s]\n')
for i in range(len(Evec)):
   fo.write('{:3E}  {:3E}\n'.format(Evec[i],Gal_vec[i]))

fo.close()
