'''
Script to estimate the chi-square
Author: Luis A. Gutiérrez Soto
26/12/2022
'''
import numpy as np
from astropy.table import Table, QTable
import matplotlib.pyplot as plt
import os
import pyCloudy as pc
from astropy.io import fits
import seaborn as sn
import argparse
import sys
from astropy import units as u
from astropy.visualization import quantity_support
from astropy.wcs import WCS
from astropy.modeling import models
from astropy.nddata import StdDevUncertainty
# specutils packages
from specutils import Spectrum1D
from specutils.analysis import line_flux
from specutils.fitting import fit_generic_continuum
from specutils import SpectralRegion
from specutils.analysis import equivalent_width
from specutils.analysis import centroid
from specutils.analysis import moment
from specutils.manipulation import noise_region_uncertainty
from specutils.fitting import estimate_line_parameters
from specutils.manipulation import extract_region
from specutils.fitting import fit_lines
from specutils.fitting import find_lines_threshold
from specutils.analysis import gaussian_sigma_width, gaussian_fwhm, fwhm, fwzi
from specutils.analysis import centroid
from specutils.fitting.continuum import fit_continuum
from specutils.fitting import fit_generic_continuum
import warnings
with warnings.catch_warnings():  # Ignore warnings
    warnings.simplefilter('ignore')
quantity_support()
sn.set_context("poster")

parser = argparse.ArgumentParser(
    description="""Reading the ouput cloudy models""")

parser.add_argument("source", type=str,
                    default="model_100000_36.58",
                    help="Name of input model ")

cmd_args = parser.parse_args()
file_ = cmd_args.source 

dir_ = '../Model-PNNew/models/N2242_Npow_He091O382Ar624_itere_R166R1705_d3000/'
model_name = file_

# Reading the Cloudy outputs in the Mod CloudyModel object
Mod = pc.CloudyModel(model_name)
#Mod = pc.CloudyModel(os.path.join(dir_, model_name))
#Mod = pc.CloudyModel(model_name)

# print(dir(Mod)) # This is the online answering way
# print(Mod.print_stats())
# print(Mod.get_ab_ion_vol_ne('O',2))

#Getting wl and flux
wl = Mod.get_cont_x(unit='Ang')
flux = Mod.get_cont_y(cont = 'total', unit = 'esAc')
#Ordered lambda and flux 
wll, flux = zip(*sorted(zip(wl, flux)))

data = Table([wll, flux], names=('Wl', 'Flux'), meta={'name': 'first table'})
mask = (data["Wl"] > 3000) & (data["Wl"] < 9000)
data_mask = data[mask]

# OUR PN
hdu = fits.open("../Spectra-lamostdr7/spec-56581-VB031N50V1_sp08-218.fits")
hdudata = hdu[0].data
wl = hdudata[2]
Flux = hdudata[0]
inve_var = hdudata[1]
sigma = 1 / np.sqrt(inve_var)

def closest(lst, K):
    '''find the closest number'''
    lst = np.array(lst)
    idx = (np.abs(lst - K)).argmin()
    return lst[idx]

#testing
# Defining units astropy
rel_flux = u.def_unit('Relative~flux')
print(rel_flux.decompose())
lamb = wl * u.AA 
flux = Flux * rel_flux
Sigma = StdDevUncertainty(sigma * rel_flux)
spec = Spectrum1D(spectral_axis=lamb, flux=flux, uncertainty=Sigma)

# dispersion per pixel 
D = 1.0002302850208247
#print(spec.spectral_axis)

# Emission lines
# lines = {"[NeIII]": 3868.760,
#          "[NeIII]H7": 3967.470,
#          "Hdelta": 4101.742,
#          "Hgamma": 4340.471,
#          "HeII": 4685.99,
#          "Hbeta": 4861.333,
#          "[OIII]4958": 4958.911,
#          "[OIII]5006": 5006.843,
#          "[FeIII]": 5412.12,
#          "Halpha": 6564.614,
#          "[ArV]": 7005.87,
#          "[ArIII]7135": 7135.80,
#          "[ArIII]7751": 7751.06,
#          "HeII8236": 8236.8
#          }

# lines = {"[NeIII]+H7": 3967.470,
#          "Hdelta": 4101.742,
#          "Hgamma": 4340.471,
#          "HeII": 4685.99,
#          "Hbeta": 4861.333,
#          "[OIII]4958": 4958.911,
#          "[OIII]5006": 5006.843,
#          "[FeIII]": 5412.12,
#          "Halpha": 6564.614,
#          }
lines = {"NE_3_396747A": 3967.470,
         "H__1_410173A": 4101.742,
         "H__1_434046A": 4340.471,
         "HE_2_468568A": 4685.99,
         "H__1_486132A": 4861.333,
         "O__3_495891A": 4958.911,
         "O__3_500684A": 5006.843,
         "HE_2_541145A": 5412.12,
         "H__1_656281A": 6564.614,
         }

table = Table.read("parameters-lamost-pn-model_130000_37.12_3.60-selec-lines.ecsv-", format="ascii.ecsv")
ratio_lines = table["Ratio Flux"]
err_ratio_lines = table["Sigma"]

#################################################################################
# Model  ########################################################################
#################################################################################
flux_lines_models = []
for vv, tt in lines.items():
    flux_lines_models.append(Mod.get_emis_vol(vv))

Hbeta_models = np.array(flux_lines_models[4])
ratio_lines_models = np.array(flux_lines_models / Hbeta_models)
print("Name model:", model_name.split("l/")[-1])
print("Lamost:", ratio_lines)
print("Err Lamost:", err_ratio_lines)
print("model:", ratio_lines_models)
# Estimating the chi-square
chi = (ratio_lines_models - ratio_lines)**2 / err_ratio_lines**2
chi_sum = chi.sum()
print("Chi:", chi)
# Estimating the degree freedom
n = 9
np = 3
vv = n - np
chi_sum_red = chi_sum / vv

modell, chii, chii_red = [], [], [] 
if chi_sum_red <= 4:
    modell.append(model_name.split("l/")[-1])
    chii.append(chi_sum)
    chii_red.append(chi_sum_red)    
print(chii, chii_red)    
tab = Table([modell, chii, chii_red],
           names=('Name model', 'chi', 'Chi red'),
           meta={'name': 'first table'})

try:
    tab.write("better-" + str(modell).split("['")[-1].split("']")[0] +".ecsv", format="ascii.ecsv", overwrite=True)
except TypeError:
    pass

