'''
Script to dealing with output data from Cloudy
Based in pyCloudy (Morisset, C., 2013, pyCloudy, Astrophysics Source Code Library)
Author: Luis A. Gutiérrez Soto
29/12/2022
'''
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
import os
import pyCloudy as pc
from astropy.io import fits
import numpy as np
import seaborn as sn
import argparse
import sys
from astropy.modeling import models
import astropy.units as u
from specutils import Spectrum1D
from specutils.analysis import line_flux
from specutils.fitting import fit_generic_continuum
from specutils import SpectralRegion
from specutils.analysis import equivalent_width
from specutils.analysis import centroid
from specutils.analysis import moment
from specutils.manipulation import noise_region_uncertainty
from specutils.fitting import find_lines_threshold
from specutils.manipulation import extract_region
from specutils.fitting import estimate_line_parameters
from specutils.fitting import fit_lines
from astropy.nddata import StdDevUncertainty
import warnings
sn.set_context("poster")


parser = argparse.ArgumentParser(
    description="""Reading the ouput cloudy models""")

parser.add_argument("source", type=str,
                    default="model_100000_36.58",
                    help="Name of input model ")


cmd_args = parser.parse_args()
file_ = cmd_args.source 


model_name = file_
# Reading the Cloudy outputs in the Mod CloudyModel object
Mod = pc.CloudyModel(model_name)

print(dir(Mod)) # This is the online answering way
print(Mod.print_stats())
print(Mod.get_ab_ion_vol_ne('O',2))

#Getting wl and flux
wl = Mod.get_cont_x(unit='Ang')
fluxx = Mod.get_cont_y(cont = 'total', unit = 'esAc')
#Ordered lambda and flux 
wll, fluxx = zip(*sorted(zip(wl, fluxx)))

data = Table([wll, fluxx], names=('Wl', 'Flux'), meta={'name': 'first table'})
mask = (data["Wl"] > 3000) & (data["Flux"] < 9000)
data_mask = data[mask]

# OUR PN
hdu = fits.open("../../../Spectra-lamostdr7/spec-56581-VB031N50V1_sp08-218.fits")
hdudata = hdu[0].data
wl = hdudata[2]
Flux = hdudata[0]


# Model
#Re-format this dataset into astropy quantities
lamb_m = data_mask["Wl"] * u.AA 
Flux_m = data_mask["Flux"] *  u.Unit('erg cm-2 s-1 AA-1')
spec_m = Spectrum1D(spectral_axis=lamb_m, flux=Flux_m)


#𝐅𝐢𝐧𝐝𝐢𝐧𝐠 𝐭𝐡𝐞 𝐥𝐢𝐧𝐞
hbeta_region = SpectralRegion((4863-50)*u.AA, (4863+50)*u.AA)
hbeta_spec_m = extract_region(spec_m, hbeta_region)
hbeta_lines_m = find_lines_threshold(hbeta_spec_m, noise_factor = 3)
print("Line center model:", hbeta_lines_m)

# Fitting Gaussian and measurent line Hbeta
# Estimate parameters
hbeta_region_m = SpectralRegion(hbeta_lines_m['line_center'] - 50 * u.AA, hbeta_lines_m['line_center'] + 50 * u.AA)
sub_spectrum_m = extract_region(spec_m, hbeta_region_m)
line_para_m = estimate_line_parameters(sub_spectrum_m, models.Gaussian1D())

print("Line parameter model:", line_para_m)

# Fit the spectrum and calculate the fitted flux values (``y_fit``)
g_init_m = models.Gaussian1D(amplitude=line_para_m.amplitude.value * u.Unit('erg cm-2 s-1 AA-1'), mean=line_para_m.mean.value * u.AA , stddev=line_para_m.stddev.value * u.AA )
g_fit_m = fit_lines(spec_m, g_init_m, window=(hbeta_lines_m['line_center'] - 50 * u.AA, hbeta_lines_m['line_center'] + 50 * u.AA))
y_fit_m = g_fit_m(lamb_m)
#plt.plot(data_mask["Wl"], data_mask["Flux"], linewidth=0.4)
#plt.plot(data_mask["Wl"], y_fit_m)
plt.xlim((4863-50), (4863+50))
#plt.xlim(3500, 9000)
#plt.title('Single fit peak window')
plt.grid(True)
#plt.show()

# Integrating
min_lamb_m = hbeta_lines_m['line_center'] - 3*line_para_m.stddev.value * u.AA
max_lamb_m = hbeta_lines_m['line_center'] + 3*line_para_m.stddev.value * u.AA
sub_region0_m = SpectralRegion(min_lamb_m,  max_lamb_m)
gauss_m = Spectrum1D(spectral_axis=lamb_m, flux=y_fit_m) 
sub_gauss_m = extract_region(gauss_m, sub_region0_m)
flux_line_Hbeta_m = np.trapz(sub_gauss_m.flux, sub_gauss_m.spectral_axis)
print("Emission line Hbeta model:", flux_line_Hbeta_m)

#Normalizing
flux_m = data_mask["Flux"] / flux_line_Hbeta_m

# Our PN
# Defining units astropy
#Re-format this dataset into astropy quantities
rel_flux = u.def_unit('Relative~flux')
print(rel_flux.decompose())
inve_var = hdudata[1]
sigma = 1 / np.sqrt(inve_var)
lamb_o = wl * u.AA 
flux_o = Flux * rel_flux
Sigma = StdDevUncertainty(sigma * rel_flux)
#Spectrun-1d
spec_o = Spectrum1D(spectral_axis=lamb_o, flux=flux_o, uncertainty=Sigma)

#𝐅𝐢𝐧𝐝𝐢𝐧𝐠 𝐭𝐡𝐞 𝐥𝐢𝐧𝐞
hbeta_spec_o = extract_region(spec_o, hbeta_region)
hbeta_lines_o = find_lines_threshold(hbeta_spec_o, noise_factor = 3)
print("Line center lamost:", hbeta_lines_o)

# Fitting Gaussian and measurent line Hbeta
# Estimate parameters
def closest(lst, K):
    '''find the closest number'''
    lst = np.array(lst)
    idx = (np.abs(lst - K)).argmin()
    return lst[idx]
wl_hbeta_o = closest(hbeta_lines_o['line_center'], 4861.333)
mask = hbeta_lines_o['line_center'] == wl_hbeta_o
hbeta_lines_o_mask = hbeta_lines_o[mask]
print("Line center lamost Hbeta:", hbeta_lines_o_mask)


hbeta_region_o = SpectralRegion(hbeta_lines_o_mask['line_center'] - 50 * u.AA, hbeta_lines_o_mask['line_center'] + 50 * u.AA)
sub_spectrum_o = extract_region(spec_o, hbeta_region_o)
line_para_o = estimate_line_parameters(sub_spectrum_o, models.Gaussian1D())
print("Line parameter obser:", line_para_o)

# Fit the spectrum and calculate the fitted flux values (``y_fit``)
g_init_o = models.Gaussian1D(amplitude=line_para_o.amplitude.value * rel_flux , mean=line_para_o.mean.value * u.AA , stddev=line_para_o.stddev.value * u.AA )
g_fit_o = fit_lines(spec_o, g_init_o, window=(hbeta_lines_o_mask['line_center'] - 50 * u.AA, hbeta_lines_o_mask['line_center'] + 50 * u.AA))
y_fit_o = g_fit_o(lamb_o)
plt.plot(wl, Flux, linewidth=0.4)
plt.plot(wl, y_fit_o)
plt.xlim((4863-30), (4863+30))
#plt.xlim(3500, 9000)
#plt.title('Single fit peak window')
plt.grid(True)
#plt.show()

# Integrating
min_lamb_o = hbeta_lines_o_mask['line_center'] - 3*line_para_o.stddev.value * u.AA
max_lamb_o = hbeta_lines_o_mask['line_center'] + 3*line_para_o.stddev.value * u.AA
sub_region0_o = SpectralRegion(min_lamb_o,  max_lamb_o)
gauss_o = Spectrum1D(spectral_axis=lamb_o, flux=y_fit_o) 
sub_gauss_o = extract_region(gauss_o, sub_region0_o)
flux_line_Hbeta_o = np.trapz(sub_gauss_o.flux, sub_gauss_o.spectral_axis)
print("Emission line Hbeta obser:", flux_line_Hbeta_o)

#Normalizing
flux_o = Flux / flux_line_Hbeta_o

fig, ax = plt.subplots(figsize=(11, 5))
#ax.set_title(namefile)
ax.set(xlim=[3600,9100])
#plt.ylim(ymin=-200,ymax=1500)
ax.set(xlabel='Wavelength $(\AA)$')
ax.set(ylabel='Flux')
plt.plot(data_mask["Wl"], flux_m,  c = "darkolivegreen", linewidth=0.7, label = 'Model')
plt.plot(wl, flux_o, c = "blueviolet", linewidth=0.7, label = 'Our PN')

ax.legend(loc="upper right")
sn.despine()
plt.tight_layout()
plt.savefig(model_name + "-hbetaInte.jpg")
