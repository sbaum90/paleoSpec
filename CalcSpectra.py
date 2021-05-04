import numpy as np
import os
import sys
from scipy.special import erf
import importlib

from paleoSpec import MWHaloSpec
from paleoSpec import BkgNeuSpec

# some constants
secPERyr = 3.154e7
dayPERMyr = 3.65e8
kmPERpc = 3.086e13

PATH_abs = os.path.dirname(os.path.abspath(__file__))

good_mineral_list = [
    "Anhydrite",
    "Baddeleyite",
    "Bieberite",
    "Bischofite",
    "Borax",
    "Cattiite",
    "Diamond",
    "Epsomite",
    "Evenkite",
    "Gypsum",
    "Halite",
    "Mirabilite",
    "Natron",
    "Nchwaningite",
    "Nickelbischofite",
    "Olivine",
    "Phlogopite",
    "Sinjarite",
    "Trona",
    "Veatchite",
]


class CalcSpectra:
    def __init__(self, mineral, switch_keep_H=False):

        self.name_mineral = mineral

        if self.name_mineral in good_mineral_list:
             self.TargetList = importlib.import_module("paleoSpec."+self.name_mineral+".TargetList")
             importlib.reload(self.TargetList) # this step is neccessary to ensure that TargetList is loaded correctly if switch_keep_H has been changed
        else:
            print("Oh no, we don't know this mineral")

        # remove Hydrogen from target list unless switch_keep_H == True
        if switch_keep_H == False and "H" in self.TargetList.nameT_list:
            indH = self.TargetList.nameT_list.index("H")
            del self.TargetList.mT_list[indH]
            del self.TargetList.AT_list[indH]
            del self.TargetList.ZT_list[indH]
            del self.TargetList.nameT_list[indH]
            self.TargetList.massFrac_list = np.delete(
                self.TargetList.massFrac_list, indH
            )

        # load ranges
        self.Trange = []
        for name in self.TargetList.nameT_list:
            self.Trange.append(
                np.loadtxt(
                    PATH_abs
                    + "/"
                    + self.name_mineral
                    + "/Ranges/"
                    + name
                    + self.name_mineral[:3]
                    + "_SRIM.txt"
                )
            )

    # -----------------------------------------------
    # smooth Milky Way halo signal calculation
    def calc_dRdx_MW(self, mDM, SIDD, rhoDM=0.3, vrel=248.0, sigv=166.0, vesc=550.0):
        """
        returns a tuple
        [track lengths, differential Rate(per unit time) of tracks
        per unit track length and unit target mass]
        summed over the mineral from the MW, assuming a standard
        truncated Maxwell-Boltzmann velocity distribution
        inputs:
        mDM - DM mass in [GeV]
        SIDD - SI cross section in [cm^2]
        ----
        rhoDM - local DM density in [GeV/cm^3]
        vrel - relative speed wrt the galactic rest frame in [km/s]
        sigv - velocity dispersion in [km/s]
        vesc - escape velocity in [km/s]
        output:
        x_T [Å]
        (dR/dx) in [1/Å/kg/Myr]
        """
        xvec = np.logspace(0, 4, 801)  # grid in track lengths for output [Å]
        # calculate recoil energy spectra for each element
        # these are not weighted by the mass fraction yet (done when converting to track length spectra)
        dRdE_MW = []
        for i in range(len(self.TargetList.nameT_list)):
            dRdE_MW.append(
                MWHaloSpec.dRdE(
                    mDM,
                    SIDD,
                    self.TargetList.mT_list[i],
                    self.TargetList.AT_list[i],
                    1.0,
                    rhoDM=rhoDM,
                    vrel=vrel,
                    sigv=sigv,
                    vesc=vesc,
                )
            )
        # calculate track length spectrum for each element (in units of 1/Å/kg/day)
        dRdx_MW = []
        for i in range(len(self.TargetList.nameT_list)):
            dRdx_MW.append(
                [
                    np.interp(
                        dRdE_MW[i][0],
                        self.Trange[i][:, 0],
                        self.Trange[i][:, 2],
                        left=0.0,
                        right=0.0,
                    ),
                    np.interp(
                        dRdE_MW[i][0],
                        self.Trange[i][:, 0],
                        self.Trange[i][:, 1],
                        left=0.0,
                        right=0.0,
                    )
                    * dRdE_MW[i][1],
                ]
            )
        # interpolate and sum
        dRdx_MW_sum = np.zeros(xvec.shape)
        for i in range(len(self.TargetList.nameT_list)):
            dRdx_MW_sum += np.interp(
                xvec, dRdx_MW[i][0], dRdx_MW[i][1], left=0.0, right=0.0
            )
        # convert units
        dRdx_MW_sum *= dayPERMyr
        return xvec, dRdx_MW_sum

    # -----------------------------------------------
    # neutrino background calculation
    def calc_dRdx_BkgNeu(self, Fnu):
        """
        returns a tuple
        [track lengths, differential Rate(per unit time) of tracks
        per unit track length and unit target mass]
        summed over the mineral from the neutrino background Fnu
        inputs:
        Fnu - neutrino flux in [1/cm^2/s/MeV] these can be found in
            BkgNeuSpec
        output:
        x_T [Å]
        (dR/dx) in [1/Å/kg/Myr]
        """
        xvec = np.logspace(0, 4, 801)  # grid in track lengths for output [Å]
        # calculate recoil energy spectra for each element
        # these are not weighted by the mass fraction yet (done when converting to track length spectra)
        dRdE_BkgNeu = []
        for i in range(len(self.TargetList.nameT_list)):
            dRdE_BkgNeu.append(
                BkgNeuSpec.dRdE(
                    self.TargetList.mT_list[i],
                    self.TargetList.AT_list[i],
                    self.TargetList.ZT_list[i],
                    1.0,
                    Fnu,
                )
            )
        # calculate track length spectrum for each element (in units of 1/Å/kg/day)
        dRdx_BkgNeu = []
        for i in range(len(self.TargetList.nameT_list)):
            dRdx_BkgNeu.append(
                [
                    np.interp(
                        dRdE_BkgNeu[i][0],
                        self.Trange[i][:, 0],
                        self.Trange[i][:, 2],
                        left=0.0,
                        right=0.0,
                    ),
                    np.interp(
                        dRdE_BkgNeu[i][0],
                        self.Trange[i][:, 0],
                        self.Trange[i][:, 1],
                        left=0.0,
                        right=0.0,
                    )
                    * dRdE_BkgNeu[i][1],
                ]
            )
        # interpolate and sum
        dRdx_BkgNeu_sum = np.zeros(xvec.shape)
        for i in range(len(self.TargetList.nameT_list)):
            dRdx_BkgNeu_sum += np.interp(
                xvec, dRdx_BkgNeu[i][0], dRdx_BkgNeu[i][1], left=0.0, right=0.0
            )
        # convert units
        dRdx_BkgNeu_sum *= dayPERMyr
        return xvec, dRdx_BkgNeu_sum

    def calc_dRdx_BkgNeu_solar(self):
        """
        returns a tuple
        [track lengths, differential Rate(per unit time) of tracks
        per unit track length and unit target mass]
        summed over the mineral from the solar neutrino background
        output:
        x_T [Å]
        (dR/dx) in [1/Å/kg/Myr]
        """
        return self.calc_dRdx_BkgNeu(BkgNeuSpec.Fnu_vec_solar)

    def calc_dRdx_BkgNeu_DSNB(self):
        """
        returns a tuple
        [track lengths, differential Rate(per unit time) of tracks
        per unit track length and unit target mass]
        summed over the mineral from the diffuse supernova
        neutrino background
        output:
        x_T [Å]
        (dR/dx) in [1/Å/kg/Myr]
        """
        return self.calc_dRdx_BkgNeu(BkgNeuSpec.Fnu_vec_DSNB)

    def calc_dRdx_BkgNeu_GSNB(self):
        """
        returns a tuple
        [track lengths, differential Rate(per unit time) of tracks
        per unit track length and unit target mass]
        summed over the mineral from the galactic supernova
        neutrino background
        output:
        x_T [Å]
        (dR/dx) in [1/Å/kg/Myr]
        """
        return self.calc_dRdx_BkgNeu(BkgNeuSpec.Fnu_vec_GSNB)

    def calc_dRdx_BkgNeu_atm(self):
        """
        returns a tuple
        [track lengths, differential Rate(per unit time) of tracks
        per unit track length and unit target mass]
        summed over the mineral from the atmospheric neutrino
        background
        output:
        x_T [Å]
        (dR/dx) in [1/Å/kg/Myr]
        """
        return self.calc_dRdx_BkgNeu(BkgNeuSpec.Fnu_vec_Atm)

    # -----------------------------------------------
    # neutron background calculation
    def calc_dRdx_Bkgn(self, C238):
        """
        returns a tuple
        [track lengths, differential Rate(per unit time) of tracks
        per unit track length and unit target mass]
        summed over the mineral from radiogenic neutrons
        input:
        C238 - uranium-238 concentration per weigth
        output:
        x_T [Å]
        (dR/dx) in [1/Å/kg/Myr]
        """
        xvec = np.logspace(0, 4, 801)  # grid in track lengths for output [Å]
        mmol238 = 238.051
        # load recoil spectrum file from neutrons
        in_dRdE = np.loadtxt(
            PATH_abs
            + "/Background_inputs/NeutronRecoilSpectra/"
            + self.name_mineral
            + "_ninduced_wan.dat"
        )
        # sort data
        lines = open(
            PATH_abs
            + "/Background_inputs/NeutronRecoilSpectra/"
            + self.name_mineral
            + "_ninduced_wan.dat",
            "r",
        ).readlines()
        header = (lines[0].split())[4:]
        del header[1], lines
        # sum over isotopes
        dRdE_n = []
        for name in self.TargetList.nameT_list:
            spec = np.copy(in_dRdE)[:, 0:2]
            spec[:, 1] = 0.0
            for isotope in header:
                if name in isotope:
                    i = header.index(isotope)
                    spec[:, 1] += in_dRdE[:, i + 1]
            dRdE_n.append(spec)
        # calculate track length spectrum for each element (in units of 1/Å/kg/day)
        dRdx_n = []
        for i in range(len(self.TargetList.nameT_list)):
            dRdx_n.append(
                [
                    np.interp(
                        dRdE_n[i][:, 0],
                        self.Trange[i][:, 0],
                        self.Trange[i][:, 2],
                        left=0.0,
                        right=0.0,
                    ),
                    np.interp(
                        dRdE_n[i][:, 0],
                        self.Trange[i][:, 0],
                        self.Trange[i][:, 1],
                        left=0.0,
                        right=0.0,
                    )
                    * dRdE_n[i][:, 1]
                    / self.TargetList.massFrac_list[i],
                ]
            )
        # interpolate and sum
        dRdx_n_sum = np.zeros(xvec.shape)
        for i in range(len(self.TargetList.nameT_list)):
            dRdx_n_sum += np.interp(
                xvec, dRdx_n[i][0], dRdx_n[i][1], left=0.0, right=0.0
            )
        # normalize to uranium concentration per weight C238
        conv_fac = self.TargetList.mmol / mmol238 * C238 / 1e-10
        dRdx_n_sum *= conv_fac
        return xvec, dRdx_n_sum

    def smear_and_bin_1a(self, C238, sigx, cutoff=0.5, xmin=-1, xmax=1e4, nbins=-1, logbins=False):
        """
        returns a tuple [bin edges, entries/bin]
        of the smeared and binned spectrum for single-alpha
        decay background
        inputs:
        C238 - uranium-238 concentration per weigth
        sigx - track length resolution
        cutoff - hard cutoff if the true track length
            in units of sigx below which tracks will
            no longer be considered in the smearing
            to avoid sensitivity relying purely
            on up-scattering of tracks
        xmin - lower edge of smallest track length
            bin returned. For the default value
            xmin=-1, the code automatically sets
            xmin = sigx/s
        xmax - upper edge of the largest track length bin
        nbins - number of bins returned. For the 
            default value nbins=-1 the code
            automatically sets 
            nbins = (xmax . xmin)/sigx
        logbins - boolean switch. For the default
            value logbins=False, the bins have linear
            width. If logbins=True, the function uses
            log spaced bins
        output:
        x_i - track length bin edges in Å
        n_i - entries/bin [1/kg]
        """
        mmol238 = 238.051 * 1e-3  # molar mass of uranium 238 in kg/mol
        NA = 6.022e23  # Avogadro number [1/mol]
        Thalf238 = 4.468e9  # half-life if uranium-238 in [yr]
        Thalf234 = 2.455e5  # half-life if uranium-234 in [yr]
        if xmin == -1:
            xmin = sigx / 2.
        if nbins == -1 and not logbins:
            nbins = int(np.floor((xmax - xmin) / sigx))
            x_i = np.linspace(xmin, xmin + sigx * nbins, nbins + 1)
        elif not logbins:
            x_i = np.linspace(xmin, xmax, nbins+1)
        else:
            x_i = np.geomspace(xmin, xmax, nbins+1)
        # make smeared spectrum
        if self.TargetList.Th_length >= cutoff * sigx:
            n_i = 0.5 * (
                erf((self.TargetList.Th_length - x_i[:-1]) / (np.sqrt(2.) * sigx))
                - erf((self.TargetList.Th_length - x_i[1:]) / (np.sqrt(2.) * sigx))
            )
        else:
            ni = np.zeros(len(x_i)-1)
        normfac = C238 * NA / mmol238 * Thalf234 / Thalf238
        return x_i, normfac * n_i


# -----------------------------------------------
# finite resolution smearing and binning
def smear_and_bin(spectrum, sigx, cutoff=0.5, xmin=-1, xmax=1e4, nbins=-1, logbins=False):
    """
    returns a tuple [bin edges, entries/bin]
    of the smeared and binned spectrum for the
    input spectrum defined by the spectrum spec
    at the track lengths xvec
    inputs:
       spectrum - tuple of (xvec, spec), where
          xvec - vector of track lengths at which
          spec - the input spectrum is calculated
       sigx - track length resolution
       cutoff - hard cutoff if the true track length
          in units of sigx below which tracks will
          no longer be considered in the smearing
          to avoid sensitivity relying purely
          on up-scattering of tracks
       xmin - lower edge of smallest track length
          bin returned. For the default value
          xmin=-1, the code automatically sets
          xmin = sigx/s
       xmax - upper edge of the largest track length bin
       nbins - number of bins returned. For the 
          default value nbins=-1 the code
          automatically sets 
          nbins = (xmax . xmin)/sigx
       logbins - boolean switch. For the default
          value logbins=False, the bins have linear
          width. If logbins=True, the function uses
          log spaced bins
    output:
       x_i - track length bin edges
       n_i - entries/bin
    note that the units of the output are determined by
    the input. For example, enter xvec, sigx, smin, xmax
    in Å, and spec in 1/Å/kg/Myr, then, the output has units
    of Å for the bin edges, and 1/kg/Myr for the entries/bin
    """
    if xmin == -1:
        xmin = sigx / 2.
    if nbins == -1 and not logbins:
        nbins = int(np.floor((xmax - xmin) / sigx))
        x_i = np.linspace(xmin, xmin + sigx * nbins, nbins + 1)
    elif not logbins:
        x_i = np.linspace(xmin, xmax, nbins+1)
    else:
        x_i = np.geomspace(xmin, xmax, nbins+1)
    # improve resolution of the spectrum where neccessary
    inds = np.where((spectrum[0][1:] - spectrum[0][:-1]) > sigx * 0.1)[0]
    xvec = np.concatenate(
        [
            spectrum[0][: inds[0]],
            np.linspace(
                spectrum[0][inds[0]],
                np.max(spectrum[0]),
                int((np.max(spectrum[0]) - spectrum[0][inds[0]]) / (0.1 * sigx)),
            ),
        ]
    )
    spec = np.interp(xvec, spectrum[0], spectrum[1], left=0.0, right=0.0)
    # remove tracks below the threshold
    spec[np.where(xvec[:-1] < cutoff * sigx)[0]] = 0.
    # compute the smeared and binned spectrum
    n_i = np.zeros(nbins)
    for i in range(nbins):
        Wvec = 0.5 * (
            erf((xvec - x_i[i]) / (np.sqrt(2.) * sigx))
            - erf((xvec - x_i[i + 1]) / (np.sqrt(2.) * sigx))
        )
        n_i[i] = np.trapz(Wvec * spec, x=xvec)
    return x_i, n_i
