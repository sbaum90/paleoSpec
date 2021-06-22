# paleoSpec
Calculates signal and background spectra for paleo detectors.

Jupyter Notebooks demonstrating how to use the package can be found in the "ExamplesAndTests" folder.

If you use **paleoSpec** for your work, please cite [arXiv:1806.05991](https://arxiv.org/abs/1806.05991), [arXiv:1811.06844](https://arxiv.org/abs/1811.06844), [arXiv:1906.05800](https://arxiv.org/abs/1906.05800), and [arXiv:2106.06559](https://arxiv.org/abs/2106.06559)

# Description
The high-level functions to compute the track length spectra are in "CalcSpectra.py". Currently, CalcSpectra allows to compute:

- the differential background rate (per unit target mass) dR/dx from the relevant neutrino flux components (calc_dRdx_BkgNeu_X), where X={solar, DSNB, GSNB atm}={all solar neutrinos, Diffuse Supernova Background, Galactic Supernova Background, Atmospheric neutrinos}

- the differential background rate (per unit target mass) dR/dx from radiogenic neutrons (calc_dRdx_Bkgn)
    
- the differential signal rate (per unit target mass) dR/dx from the local DM density (calc_dRdx_MW)

These functions return a tuple (x, dR/dx), where the first entry is a vector of track lengths, and the second entry is a vector containing the spectrum (at these track length). Note that depending on the source of the spectrum, these functions return either the differential rate of events, or the differential number of events! The inputs to the respective functions are described in the script.

To account for finite resolution effects, each spectrum produced by the above functions can be fed into the "smear_and_bin" function. This function accounts for finite resolution effects assuming a gaussian distribution of the track length error. It returns a tuple of (bin edges, entries/bin). Note that the input format of "smear_and_bin" is set up to take the output of the different functions to compute the spectra above. For example, to compute the smeared and binned spectrum from the local DM signal, call smear_and_bin(calc_dRdx_MW(args), sigx) where sigx is the track length resolution. Bins are chosen automatically, see the args of "smear_and_bin" and the comments for the input kwargs that allow the user to control the binning.

Note that there is a special function to get the smeared-and-binned spectrum from the single-alpha background (smear_and_bin_1a).

-------------------------------
For each spectrum, there is a script which contains the recoil energy spectrum calculation:

BkgNeuSpec.py has the relevant functions for calculating the recoil energy spectrum induced by neutrinos

MWHaloSpec.py has the relevant functions for calculating the recoil energy spectrum from the local DM density

helper_functions.py contains functions for the mean inverse speed of a Maxwell-Boltzmann distribution and the Helm form factor
# Minerals
Currently supported Minerals:
- Anhydrite
- Baddeleyite
- Bieberite
- Bischofite
- Borax
- Cattiite
- Diamond
- Epsomite
- Evenkite
- Gypsum
- Halite
- Mirabilite
- Natron
- Nchwaningite
- Nickelbischofite
- Olivine
- Phlogopite
- Sinjarite
- Trona
- Veatchite
