'''
Script to estimate the chi-square
Author: Luis A. Gutiérrez Soto
26/12/2022
'''
import numpy as np
from astropy.table import Table
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
import warnings
with warnings.catch_warnings():  # Ignore warnings
    warnings.simplefilter('ignore')
quantity_support()
sn.set_context("poster")

# parser = argparse.ArgumentParser(
#     description="""Reading the ouput cloudy models""")

# parser.add_argument("source", type=str,
#                     default="model_100000_36.58",
#                     help="Name of input model ")


# cmd_args = parser.parse_args()
# file_ = cmd_args.source 


# model_name = file_
# # Reading the Cloudy outputs in the Mod CloudyModel object
# Mod = pc.CloudyModel(model_name)

# print(dir(Mod)) # This is the online answering way
# print(Mod.print_stats())
# print(Mod.get_ab_ion_vol_ne('O',2))

# #Getting wl and flux
# wl = Mod.get_cont_x(unit='Ang')
# flux = Mod.get_cont_y(cont = 'total', unit = 'esAc')
# #Ordered lambda and flux 
# wll, flux = zip(*sorted(zip(wl, flux)))

# data = Table([wll, flux], names=('Wl', 'Flux'), meta={'name': 'first table'})
# mask = (data["Wl"] > 3000) & (data["Flux"] < 9000)
# data_mask = data[mask]

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

def find_line(wl_vacuum, spec):
    '''
    This function find the line
    in a wavelenght range
    '''
    line_region = SpectralRegion((wl_vacuum-50)*u.AA, (wl_vacuum+50)*u.AA)
    line_spec = extract_region(spec, line_region)
    with warnings.catch_warnings():  # Ignore warnings
        warnings.simplefilter('ignore')
        line_spec = find_lines_threshold(line_spec, noise_factor = 3) # find the lines
        print("Number of line found:", len(line_spec))
    if len(line_spec) > 1:
        wl_line = closest(line_spec['line_center'].value, wl_vacuum)
        mask = line_spec['line_center'].value == wl_line
        line_spec_mask = line_spec[mask]
    else:
        line_spec_mask = line_spec
    return line_spec_mask

def measurent_line(wl_vacuum, spec, lamb, wl, Flux, NameLine):
    line_spec_mask = find_line(wl_vacuum,  spec)
    #Extract again the region using the lice center found
    line_region_ = SpectralRegion(line_spec_mask['line_center'] - 5 * u.AA, line_spec_mask['line_center'] + 5 * u.AA)
    sub_spectrum_line = extract_region(spec, line_region_)
    line_para_line = estimate_line_parameters(sub_spectrum_line, models.Gaussian1D())
    print("Parameters of the 1D-Gaussin:", line_para_line)
    # Fit the spectrum and calculate the fitted flux values (``y_fit``)
    g_init_line = models.Gaussian1D(amplitude=line_para_line.amplitude.value * rel_flux,
                                    mean=line_para_line.mean.value * u.AA , stddev=line_para_line.stddev.value * u.AA )
    g_fit_line = fit_lines(spec, g_init_line, window=(line_spec_mask['line_center'] - 5 * u.AA, line_spec_mask['line_center'] + 5 * u.AA))
    y_fit_line = g_fit_line(lamb)
    #Integrating along the fit 1D-Gaussian
    gauss = Spectrum1D(spectral_axis=lamb, flux=y_fit_line) 
    sub_gauss = extract_region(gauss, line_region_)
    min_lamb = line_para_line.mean.value - 3*line_para_line.stddev.value
    max_lamb = line_para_line.mean.value + 3*line_para_line.stddev.value
    sub_region_int = SpectralRegion(min_lamb * u.AA,  max_lamb * u.AA)
    sub_gauss_int = extract_region(gauss, sub_region_int)
    flux_line = np.trapz(sub_gauss_int.flux, sub_gauss_int.spectral_axis) 
    #Ploting the lina and fit Gaussian
    fig, ax = plt.subplots(figsize=(12, 12))
    plt.plot(wl, Flux, linewidth=5, label = "Obs.")
    plt.plot(wl, y_fit_line, label = "Model")
    plt.tight_layout()
    #plt.ylim(-100, (sub_spectrum_line.max() + 500*rel_flux))
    plt.xlim((line_spec_mask['line_center'].value-30), (line_spec_mask['line_center'].value+30))
    bbox_props = dict(boxstyle="round", fc="w", ec="0.88", alpha=0.6, pad=0.1)
    plt.text(0.1, 0.9, NameLine,
         transform=ax.transAxes, c="black", weight='bold', fontsize=24.8, bbox=bbox_props)
    ax.legend()
    plt.savefig(NameLine + ".pdf")
    plt.close()
    return flux_line

def err_line():
    '''
    Estimating the error of the emission line flux,
    I am using the eq. of Tresse et al. 1999
    ''' 

    

#testing
# Defining units astropy
rel_flux = u.def_unit('Relative~flux')
print(rel_flux.decompose())
lamb = wl * u.AA 
flux = Flux * rel_flux
Sigma = StdDevUncertainty(sigma * rel_flux)
spec = Spectrum1D(spectral_axis=lamb, flux=flux, uncertainty=Sigma)
#print(spec.uncertainty)

lines = {"[NeIII]": 3868.760,
         "[NeIII]H7": 3967.470,
         "Hdelta": 4101.742,
         "Hgamma": 4340.471,
         "HeII": 4685.99,
         "Hbeta": 4861.333,
         "[OIII]4958": 4958.911,
         "[OIII]5006": 5006.843,
         "Halpha": 6564.614,
         "[FeIII]": 5412.12,
         "[ArV]": 7005.87,
         "[ArIII]7135": 7135.80,
         "[ArIII]7751": 7751.06,
         "HeII8236": 8236.8
         }

flux_line = []
for v, t in lines.items():
    flux_line.append(measurent_line(t, spec, lamb, wl, Flux, v))
print(flux_line[8])

# Halpha = 6564.614
# Halpha_ = measurent_line(Halpha, spec, lamb, wl, Flux, "Halpha")
# # He II
# HeII = 4685.99
# HeII_ = measurent_line(HeII, spec, lamb, wl, Flux, "HeII")

# # Model
# wl_hbeta = closest(data_mask["Wl"], 4861.333)
# MaskHbeta = data_mask["Wl"] == wl_hbeta
# HBeta = data_mask[MaskHbeta]
# flux_m = data_mask["Flux"] / HBeta["Flux"]


# # Our PN
# wl_hbeta_our = closest(wl, 4861.333)

# MaskHbeta_our = wl == wl_hbeta_our
# flux_HBeta_our = Flux[MaskHbeta_our]
# flux_our = Flux / flux_HBeta_our

# # Estimating chi-square
# #{\displaystyle \chi ^{2}=\sum _{i=1}^{n}{\frac {(O_{i}-E_{i})^{2}}{E_{i}}}}

# # mask for the wavelenght model
# mask_wlm = (data_mask["Wl"] >= 4000) & (data_mask["Wl"] <= 5700)
# wlm = data_mask["Wl"][mask_wlm]
# print("Number of wavelenths of model:", len(wlm))

# # mask for the wavelenght observed
# mask_wlo = (wl >= 4000) & (wl <= 5700)
# wlo = wl[mask_wlo]
# print("Number of wavelenths of observed:", len(wlo))

# fig, ax = plt.subplots(figsize=(11, 5))
# #ax.set_title(namefile)
# ax.set(xlim=[3600,9100])
# #plt.ylim(ymin=-200,ymax=1500)
# ax.set(xlabel='Wavelength $(\AA)$')
# ax.set(ylabel='Flux')
# plt.plot(wlm, flux_m[mask_wlm],  c = "darkolivegreen", linewidth=0.7, label = 'Model')
# plt.plot(wlo, flux_our[mask_wlo], c = "blueviolet", linewidth=0.7, label = 'Our PN')

# ax.legend(loc="upper right")
# sn.despine()
# plt.tight_layout()
# plt.savefig(model_name + "limited.jpg")
