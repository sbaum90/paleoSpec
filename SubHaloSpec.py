import numpy as np
import matplotlib.pyplot as plt

from scipy.special import erf

# constants
G_N_pcMsolkms = 4.3009e-3  # [pc/Msol.km^2/s^2]
# cosmology
cosmo_H0 = 67.32  # [km/s/Mpc]
cosmo_Omega_mat = 0.3158
cosmo_Omega_lam = 1.0 - cosmo_Omega_mat

# unit conversions
secPERyr = 3.154e7
kmPERpc = 3.086e13
GeVPERMsol = 1.115e57

##############################################
from paleoSpec.helper_functions import eta, F2helm

# functions for NFW profile
def rho_NFW(r, rhos, rs):
    """
    returns density of a NFW profile with
    characteristic density rhos and scale radius rs
    note that the units for r are set by rs,
    and the density is returned in units of rhos
    """
    return rhos / (r / rs * (1.0 + r / rs) ** 2)


def M_r_NFW(r, rhos, rs):
    """
    returns the mass inside radius r for a NFW profile
    with characteristic density rhos and scale radius rs
    note that the units for r are set by rs,
    and the enclosed mass is returned in units of rhos*r^3
    """
    return np.where(r / rs > 1e-4,
                    4.0 * np.pi * rhos * rs ** 3 * (np.log((r + rs) / rs) - r / (r + rs)),
                    2.0 * np.pi * rhos * rs * r ** 2)


def pot_r_NFW(r, rhos, rs):
    """
    returns the gravitational potential at radius r for a
    NFW profile with characteristic density rhos and
    scale radius rs
    units: returns potential in [km^2/s^2]
       enter r, rs in [pc]
       enter rhos in [Msol/pc^3]
    """
    return (
        -4.0
        * np.pi
        * G_N_pcMsolkms
        * rhos
        * rs ** 2
        * (rs / r * (np.log((r + rs) / rs) - r / (r + rs)) + rs / (r + rs))
    )


def sigv_NFW(r, rhos, rs):
    """
    returns the 1d-velocity dispersion at radius r
    in a NFW halo with characteristic density rhos
    and scale radius rs, assuming that the halo is
    virialized at all radii
    units: returns velocity dispersion [km/s]
       enter r, rs in [pc]
       enter rhos in [Msol/pc^3]
    """
    return np.sqrt(G_N_pcMsolkms * M_r_NFW(r, rhos, rs) / r)


def vesc_NFW(r, rhos, rs):
    """
    returns the escape velocity at radius r
    in a NFW halo with characteristic density rhos
    and scale radius rs, assuming that the halo is
    virialized at all radii
    units: returns escape velocity [km/s]
       enter r, rs in [pc]
       enter rhos in [Msol/pc^3]
    """
    return np.sqrt(2 * np.abs(pot_r_NFW(r, rhos, rs)))


def rho_crit(z):
    """
    returns critical density at redshift z in [Msol/pc^3]
    """
    return (
        3.0
        * 1e-12
        * cosmo_H0 ** 2
        / (8.0 * np.pi * G_N_pcMsolkms)
        * (cosmo_Omega_mat * (1.0 + z) ** 3 + cosmo_Omega_lam)
    )


def Dvir(z):
    """
    returns critical overdensity for a halo to decouple
    from cosmic expansion at redshift z
    """
    d = (
        cosmo_Omega_mat
        * (1.0 + z) ** 3
        / (cosmo_Omega_mat * (1.0 + z) ** 3 + cosmo_Omega_lam)
        - 1.0
    )
    return 18.0 * np.pi ** 2 + 82.0 * d - 39.0 * d ** 2


def rvir_NFW(Mvir, z):
    """
    returns the virial radius in [pc] for a NFW halo
    with virial mass Mvir [Msol] forming at redshift z
    """
    return (3.0 / (4.0 * np.pi) * Mvir / (Dvir(z) * rho_crit(z))) ** (1.0 / 3.0)


def rs_NFW(Mvir, cvir, z):
    """
    returns the scale radius in [pc] for a NFW halo
    with virial mass Mvir [Msol], concentration parameters
    cvir, forming at redshift z
    """
    return rvir_NFW(Mvir, z) / cvir


def f(c):
    return np.log(1.0 + c) - c / (1.0 + c)


def rhos_NFW(Mvir, cvir, z):
    """
    returns the characteristic denisty [Msol/pc^3]
    for a NFW halo with virial mass Mvir [Msol],
    concentration parameters cvir, forming at redshift z
    """
    return cvir ** 3 / f(cvir) * Dvir(z) * rho_crit(z)


def mk_xvec(x0, x1, bsh, xN):
    """
    a very ugly function to return a 
    ~evenly spaced thing on a log-scale.
    this is used to prepare spatial grid 
    for the numerical integral over the 
    DM distribution
    """
    xN = int(xN)
    if x0 > 0 and x1 > 0:
        x = np.logspace(np.log10(x0), np.log10(x1), xN)
    elif x0 < 0 and x1 < 0:
        x = -np.logspace(np.log10(-x0), np.log10(-x1), xN)
    elif x0 < 0 and x1 > 0:
        if np.abs(x0) < np.abs(x1):
            xN1 = int(np.log10(-x0) / np.log10(x1) * xN / 2.0)
            xN2 = xN - xN1
        else:
            xN2 = int(np.log10(x1) / np.log10(-x0) * xN / 2.0)
            xN1 = xN - xN2
        xN1 = np.max([xN1, 2])
        xN2 = np.max([xN2, 2])
        xmin = bsh/xN
        x = np.concatenate(
            [
                -np.logspace(np.log10(-x0), np.log10(xmin), xN1),
                np.zeros(1),
                np.logspace(np.log10(xmin), np.log10(x1), xN2),
            ]
        )
    elif x0 > 0 and x1 < 0:
        if np.abs(x0) < np.abs(x1):
            xN1 = int(np.log10(x0) / np.log10(-x1) * xN / 2.0)
            xN2 = xN - xN1
        else:
            xN2 = int(np.log10(-x1) / np.log10(x0) * xN / 2.0)
            xN1 = xN - xN2
        xN1 = np.max([xN1, 2])
        xN2 = np.max([xN2, 2])
        xmin = bsh/xN
        x = np.concatenate(
            [
                np.logspace(np.log10(x0), np.log10(xmin), xN1),
                np.zeros(1),
                -np.logspace(np.log10(xmin), np.log10(-x1), xN2),
            ]
        )
    return x


def dndE(mDM, SIDD, mN, AN, xiN, Msh, csh, bsh, vsh, x0, x1, zsh=0.0):
    """
    returns a tuple
    [recoil energies, differential number of recoil events
    per unit recoil energy and unit target mass]
    for the target nucleus N from flying through a mini-halo
    inputs:
       mDM - DM mass in [GeV]
       SIDD - SI cross section in [cm^2]
       ---
       mN - mass of N in [GeV]
       AN - atomic number of N
       xiN - mass fraction N comprises of the target
       ----
       Msh - mass of minihalo in [Msol]
       csh - concentration parameter of minihalo
       bsh - impact parameter for minihalo [pc]
       vsh - relative speed wrt minihalo [km/s]
       x0 - distance to closest approach at start [pc]
       x1 - distance to closest approach at end [pc]
    output:
       E [keV]
       (dn/dE)_N in [1/keV/kg]
    """
    # constants
    mp = 0.938  # proton mass in GeV
    c = 2.997e5  # speed of light in km/s
    # unit conversion
    kmPERcm = 1e-5
    hbarc_GeVpc = 6.384e-33
    hbarc_keVs = 6.582e-19
    hbarc_kgs = 1.173e-51
    # let's go
    Evec = np.logspace(-3, 3, 1201)  # grid of recoil energies in keV
    dndE = np.zeros(len(Evec))  # create variable for result
    vmin = c * np.sqrt(1e-6 * Evec * (mN + mDM) ** 2 / (2.0 * mN * mDM ** 2))
    x = mk_xvec(x0, x1, bsh, 1e3)
    # get all relevant quantities from the halo
    r = np.sqrt(x ** 2 + bsh ** 2)
    rhos = rhos_NFW(Msh, csh, zsh)
    rs = rs_NFW(Msh, csh, zsh)
    rhox = rho_NFW(r, rhos, rs)
    sigvx = sigv_NFW(r, rhos, rs)
    vescx = vesc_NFW(r, rhos, rs)
    # compute the recoil spectrum
    prefac = xiN * AN ** 2 * SIDD / (2.0 * mDM ** 3 * mp ** 2 / (mDM + mp) ** 2 * vsh)
    unitconv = kmPERcm ** 2 * hbarc_GeVpc ** 2 * GeVPERMsol / hbarc_keVs / hbarc_kgs
    for i in range(len(Evec)):
        dndE[i] = F2helm(2.0 * 1e6 * mN * Evec[i], AN) * np.trapz(
            rhox * eta(vmin[i], vsh, sigvx, vescx), x=x
        )
    dndE = unitconv * prefac * dndE
    return Evec, dndE
