import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.special import erf

###########################################
# load neutrino spectra
###########################################
PATH_abs = os.path.dirname(os.path.abspath(__file__))
flist_neuflux = [
    PATH_abs + "/Background_inputs/NeutrinoFluxes/CAJOH_pp.dat",
    PATH_abs + "/Background_inputs/NeutrinoFluxes/CAJOH_hep.dat",
    PATH_abs + "/Background_inputs/NeutrinoFluxes/CAJOH_8B.dat",
    PATH_abs + "/Background_inputs/NeutrinoFluxes/CAJOH_13N.dat",
    PATH_abs + "/Background_inputs/NeutrinoFluxes/CAJOH_15O.dat",
    PATH_abs + "/Background_inputs/NeutrinoFluxes/CAJOH_17F.dat",
    PATH_abs + "/Background_inputs/NeutrinoFluxes/CAJOH_lines.dat",
    PATH_abs + "/Background_inputs/NeutrinoFluxes/DSNB_spec.dat",
    PATH_abs + "/Background_inputs/NeutrinoFluxes/Galactic_SN_spec.dat",
    PATH_abs + "/Background_inputs/NeutrinoFluxes/CAJOH_Atm.dat",
]

Enu_vec = np.logspace(-3, 6, 9001)
Fnu_vec_tot = np.zeros(Enu_vec.shape)
# interpolate and sum solar fluxes
Fnu_vec_solar = np.zeros(Enu_vec.shape)
for i in range(0, 6):
    in_flux = np.loadtxt(flist_neuflux[i])
    Fnu_vec_solar += np.exp(
        np.interp(
            np.log(Enu_vec),
            np.log(in_flux[:, 0]),
            np.log(in_flux[:, 1]),
            left=-999.0,
            right=-999.0,
        )
    )

# add monochromatic spectra
in_lines = np.loadtxt(flist_neuflux[6])
for i in range(len(in_lines[:, 0])):
    Fnu_vec_temp = np.zeros(Enu_vec.shape)
    Eind = np.argmin(np.abs(Enu_vec - in_lines[i, 0]))
    Fnu_vec_temp[Eind] = 1.0
    Fnu_vec_solar += in_lines[i, 1] / np.trapz(Fnu_vec_temp, x=Enu_vec) * Fnu_vec_temp

Fnu_vec_tot += Fnu_vec_solar

# DSNB flux
in_flux = np.loadtxt(flist_neuflux[7])
Fnu_vec_DSNB = np.exp(
    np.interp(
        np.log(Enu_vec),
        np.log(in_flux[:, 0]),
        np.log(in_flux[:, 1]),
        left=-999.0,
        right=-999.0,
    )
)
Fnu_vec_tot += Fnu_vec_DSNB

# Galactic SNB flux
in_flux = np.loadtxt(flist_neuflux[8])
Fnu_vec_GSNB = np.exp(
    np.interp(
        np.log(Enu_vec),
        np.log(in_flux[:, 0]),
        np.log(in_flux[:, 1]),
        left=-999.0,
        right=-999.0,
    )
)
Fnu_vec_tot += Fnu_vec_GSNB

# atmospheric flux
in_flux = np.loadtxt(flist_neuflux[9])
Fnu_vec_Atm = np.exp(
    np.interp(
        np.log(Enu_vec),
        np.log(in_flux[:, 0]),
        np.log(in_flux[:, 1]),
        left=-999.0,
        right=-999.0,
    )
)
Fnu_vec_tot += Fnu_vec_Atm


###########################################
# functions
###########################################
from paleoSpec.helper_functions import F2helm


def diff_xsec(ER, Enu, mN, AN, ZN):
    """
    returns the differential scattering cross section
    per nuclear recoil energy ER for a neutrino with
    energy Enu scattering elastically (CEvNS) off a
    nucleus with N
    input:
       ER - nuclear recoil energy in [keV]
       Enu - neutrino energy in [MeV]
       mN - mass of N in [GeV]
       AN - atomic number of N
       ZN - number of protons in N
    output:
       differential cross section in [cm^2/keV]
    """
    # constants
    GF = 1.1663787e-5  # Fermi constant in [1/GeV^2]
    sin2W = 0.23121  # weak mixing angle
    hbarc_GeVcm = 1.97e-14  # in [GeV cm]
    # unit conversion
    GeVPERkeV = 1e-6
    #
    unitconv = hbarc_GeVcm ** 2 * GeVPERkeV
    prefac = GF ** 2 / (4.0 * np.pi)
    QN2 = (AN - ZN - (1.0 - 4.0 * sin2W) * ZN) ** 2
    xsec = (
        QN2 * mN * (1.0 - mN * ER / (2.0 * Enu ** 2)) * F2helm(2.0 * 1e6 * mN * ER, AN)
    )
    xsec = np.clip(xsec, 0, 1e30)
    return unitconv * prefac * xsec


def dRdE(mN, AN, ZN, xiN, Fnu):
    """
    returns a tuple
    [recoil energies, differential rate of recoil events
    per unit recoil energy and unit target mass]
    for the target nucleus N from the neutrino background
    inputs:
       mN - mass of N in [GeV]
       AN - atomic number of N
       ZN - number of protons in N
       xiN - mass fraction N comprises of the target
       ---- optional
       Fnu - spectral neutrino flux,
          default is the total flux computed above
    output:
       E [keV]
       (dR/dE)_N in [dru = 1/keV/kg/day]
    """
    # unit conversion
    GeVPERkg = 5.61e26
    sPERday = 8.64e4
    # let's go
    unitconv = GeVPERkg * sPERday
    prefac = xiN / mN
    Evec = np.logspace(-3, 3, 1201)  # grid of recoil energies in keV
    dRdE = np.zeros(len(Evec))  # create variable for result
    for i in range(len(dRdE)):
        Enu_min = np.sqrt(mN * Evec[i] / 2.0)
        inds = np.where(Enu_vec > Enu_min)[0]
        integrand = diff_xsec(Evec[i], Enu_vec[inds], mN, AN, ZN) * Fnu[inds]
        dRdE[i] = np.trapz(integrand, x=Enu_vec[inds])
    dRdE = unitconv * prefac * dRdE
    return Evec, dRdE


def dRdE_solar(mN, AN, ZN, xiN):
    return dRdE(mN, AN, ZN, xiN, Fnu_vec_solar)


def dRdE_DSNB(mN, AN, ZN, xiN):
    return dRdE(mN, AN, ZN, xiN, Fnu_vec_DSNB)


def dRdE_GSNB(mN, AN, ZN, xiN):
    return dRdE(mN, AN, ZN, xiN, Fnu_vec_GSNB)


def dRdE_atm(mN, AN, ZN, xiN):
    return dRdE(mN, AN, ZN, xiN, Fnu_vec_Atm)
