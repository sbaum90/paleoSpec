import numpy as np
import matplotlib.pyplot as plt

from scipy.special import erf

# unit conversions
kmPERpc = 3.086e13
GeVPERMsol = 1.115e57

# functions
from paleoSpec.helper_functions import eta, F2helm


def eta_av_delta(vmin, theta, vsun, vearth=29.8):
    """
    returns the inverse mean speed averaged
    over on sidereal year
    For now, the velocity dispersion of the dark disk
    is assumed to be a delta function (in the rest
    frame of the dark disk) - can be replaced by
    replacing etav
    input:
       vmin - minima velocity for eta integral
       theta - angle between the motion of the solar system relative to
          the dark disk and the orbital plane in [rad]
       vsun - relative speed of the solar system with respect
          to the dark disk
       vearth - relative speed of the Earth around the Sun
    units:
       input theta in [rad]
       input all speeds in same units (default vearth is in km/s)
       returns averaged eta in inverse units of the speed
    """
    phiv = np.linspace(0, 2.0 * np.pi, 1000)
    vrel = np.sqrt(
        2.0 * vsun * vearth * np.sin(phiv) * np.cos(theta) + vearth ** 2 + vsun ** 2
    )
    etav = 1.0 / vrel * np.heaviside(vrel / vmin - 1.0, 0.0)
    return 1.0 / (2.0 * np.pi) * np.trapz(etav, x=phiv)


def eta_av_MB(vmin, theta, vsun, sigvDD, vearth=29.8, vescDD=100.0):
    """
    returns the inverse mean speed averaged
    over on sidereal year
    assumes the velocity distribution of the dark
    disk to be a Maxwell-Boltzmann distribution
    with inputs sigvDD and vescDD
    input:
       vmin - minima velocity for eta integral
       theta - angle between the motion of the solar system relative to
          the dark disk and the orbital plane in [rad]
       vsun - relative speed of the solar system with respect
          to the dark disk
       sigvDD - velocity dispersion of the dark disk in [km/s]
       vearth - relative speed of the Earth around the Sun
       vescDD - escape velocity of the dark disk in [km/s]
    units:
       input theta in [rad]
       input all speeds in same units (default vearth is in km/s)
       returns averaged eta in inverse units of the speed
    """
    phiv = np.linspace(0, 2.0 * np.pi, 1000)
    vrel = np.sqrt(
        2.0 * vsun * vearth * np.sin(phiv) * np.cos(theta) + vearth ** 2 + vsun ** 2
    )
    etav = eta(vmin, vrel, sigvDD, vescDD)
    return 1.0 / (2.0 * np.pi) * np.trapz(etav, x=phiv)


def dndE(mDM, SIDD, mN, AN, xiN, SigDD, vvDD, thetavDD, thetaorbitDD, sigvDD):
    """
    returns a tuple
    [recoil energies, differential number of recoil events
    per unit recoil energy and unit target mass]
    for the target nucleus N from crossing through the dark disk
    inputs:
       mDM - DM mass in [GeV]
       SIDD - SI cross section in [cm^2]
       ---
       mN - mass of N in [GeV]
       AN - atomic number of N
       xiN - mass fraction N comprises of the target
       ----
       SigDD - surface density of the dark disk in [Msol/pc^2]
       vvDD - vertical speed relative to the dark disk in [km/s]
       thetavDD - crossing angle of the solar system relative
          to the dark disk in [rad]
       thetaorbitDD - angle between the orbital plane of the Earth
          around the Sun and the velocity of the solar system
          relative to the dark disk in [rad]
       sigvDD - velocity dispersion of the dark disk in [km/s]
    output:
       E [keV]
       (dn/dE)_N in [1/keV/kg]
    """
    # constants
    mp = 0.938  # proton mass in GeV
    c = 2.997e5  # speed of light in km/s
    # units
    kmPERcm = 1e-5
    GeVPERkeV = 1e-6
    GeVPERkg = 5.61e26
    # let's go
    Evec = np.logspace(-3, 3, 1201)  # grid of recoil energies in keV
    dndE = np.zeros(len(Evec))  # create variable for result
    vmin = c * np.sqrt(1e-6 * Evec * (mN + mDM) ** 2 / (2.0 * mN * mDM ** 2))
    # compute the recoil spectrum
    prefac = (
        xiN
        * AN ** 2
        * SIDD
        / (2.0 * mDM ** 3 * mp ** 2 / (mDM + mp) ** 2)
        * SigDD
        / vvDD
    )
    unitfac = kmPERcm ** 2 * c ** 2 / kmPERpc ** 2 * GeVPERMsol * GeVPERkeV * GeVPERkg
    for i in range(len(Evec)):
        dndE[i] = F2helm(2.0 * 1e6 * mN * Evec[i], AN) * eta_av_MB(
            vmin[i], thetaorbitDD, vvDD / np.cos(thetavDD), sigvDD
        )
    dndE = unitfac * prefac * dndE
    return Evec, dndE
