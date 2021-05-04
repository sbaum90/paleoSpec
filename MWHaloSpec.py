import numpy as np
import matplotlib.pyplot as plt

from scipy.special import erf

from paleoSpec.helper_functions import eta, F2helm

def dRdE(mDM, SIDD, mN, AN, xiN, rhoDM=0.3, vrel=248.0, sigv=166.0, vesc=550.0):
    """
    returns a tuple
    [recoil energies, differential rate of recoil events
    per unit recoil energy and unit target mass]
    for the target nucleus N from the MW halo, assuming a standard
    truncate Maxwell-Boltzmann velocity distribution
    inputs:
       mDM - DM mass in [GeV]
       SIDD - SI cross section in [cm^2]
       ---
       mN - mass of N in [GeV]
       AN - atomic number of N
       xiN - mass fraction N comprises of the target
       ----
       rhoDM - local DM density in [GeV/cm^3]
       vrel - relative speed wrt the galactic rest frame in [km/s]
       sigv - velocity dispersion in [km/s]
       vesc - escape velocity in [km/s]
    output:
       E [keV]
       (dR/dE)_N in [dru = 1/keV/kg/day]
    """
    # constants
    mp = 0.938  # proton mass in GeV
    c = 2.997e5  # speed of light in km/s
    # unit conversion
    kmPERcm = 1e-5
    GeVPERkg = 5.61e26
    sPERday = 8.64e4
    GeVPERkeV = 1e-6
    # let's go
    Evec = np.logspace(-3, 3, 1201)  # grid of recoil energies in keV
    dRdE = np.zeros(len(Evec))  # create variable for result
    # calculate
    vmin = c * np.sqrt(1e-6 * Evec * (mN + mDM) ** 2 / (2.0 * mN * mDM ** 2))
    prefac = xiN * AN ** 2 * SIDD * rhoDM / (2.0 * mDM ** 3 * mp ** 2 / (mDM + mp) ** 2)
    unitconv = c ** 2 * GeVPERkeV * GeVPERkg * sPERday / kmPERcm
    dRdE = F2helm(2.0 * 1e6 * mN * Evec, AN) * eta(vmin, vrel, sigv, vesc)
    dRdE = unitconv * prefac * dRdE
    return Evec, dRdE
