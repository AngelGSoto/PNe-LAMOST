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

def measurent_line(wl_vacuum, spec, lamb, wl, Flux, units_flux, type_spec, NameLine, saveplot = "y"):
    '''
    Meassurent the flux line
    '''
    line_spec_mask = find_line(wl_vacuum,  spec)
    #Extract again the region using the lice center found
    line_region_ = SpectralRegion(line_spec_mask['line_center'] - 5 * u.AA, line_spec_mask['line_center'] + 5 * u.AA)
    sub_spectrum_line = extract_region(spec, line_region_)
    line_para_line = estimate_line_parameters(sub_spectrum_line, models.Gaussian1D())
    print("Parameters of the 1D-Gaussin:", line_para_line)
    # Fit the spectrum and calculate the fitted flux values (``y_fit``)
    g_init_line = models.Gaussian1D(amplitude=line_para_line.amplitude.value * units_flux,
                                    mean=line_para_line.mean.value * u.AA , stddev=line_para_line.stddev.value * u.AA )
    g_fit_line = fit_lines(spec, g_init_line, window=(line_spec_mask['line_center'] - 5 * u.AA, line_spec_mask['line_center'] + 5 * u.AA))
    y_fit_line = g_fit_line(lamb)
    #Integrating along the fit 1D-Gaussian
    gauss = Spectrum1D(spectral_axis=lamb, flux=y_fit_line) 
    sub_gauss = extract_region(gauss, line_region_)
    min_lamb = line_para_line.mean.value - 5*line_para_line.stddev.value
    max_lamb = line_para_line.mean.value + 5*line_para_line.stddev.value
    sub_region_int = SpectralRegion(min_lamb * u.AA,  max_lamb * u.AA)
    sub_gauss_int = extract_region(gauss, sub_region_int)
    flux_line = np.trapz(sub_gauss_int.flux, sub_gauss_int.spectral_axis) 
    #Ploting the lina and fit Gaussian
    if saveplot == "y":
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
        plt.savefig(type_spec+ "_" + NameLine + ".pdf")
        plt.close()
    else:
        pass
    return flux_line

def measurent_line_model(wl_vacuum, spec, lamb, wl, Flux, units_flux, type_spec, NameLine, saveplot = "y"):
    '''
    Meassurent the flux line
    '''
    line_spec_mask = find_line(wl_vacuum,  spec)
    #Extract again the region using the lice center found
    line_region_ = SpectralRegion(line_spec_mask['line_center'] - 30 * u.AA, line_spec_mask['line_center'] + 30 * u.AA)
    sub_spectrum_line = extract_region(spec, line_region_)
    line_para_line = estimate_line_parameters(sub_spectrum_line, models.Gaussian1D())
    print("Parameters of the 1D-Gaussin:", line_para_line)
    # Fit the spectrum and calculate the fitted flux values (``y_fit``)
    g_init_line = models.Gaussian1D(amplitude=line_para_line.amplitude.value * units_flux,
                                    mean=line_para_line.mean.value * u.AA , stddev=line_para_line.stddev.value * u.AA )
    g_fit_line = fit_lines(spec, g_init_line, window=(line_spec_mask['line_center'] - 30 * u.AA, line_spec_mask['line_center'] + 30 * u.AA))
    y_fit_line = g_fit_line(lamb)
    #Integrating along the fit 1D-Gaussian
    gauss = Spectrum1D(spectral_axis=lamb, flux=y_fit_line) 
    sub_gauss = extract_region(gauss, line_region_)
    min_lamb = line_para_line.mean.value - 8*line_para_line.stddev.value
    max_lamb = line_para_line.mean.value + 8*line_para_line.stddev.value
    sub_region_int = SpectralRegion(min_lamb * u.AA,  max_lamb * u.AA)
    sub_gauss_int = extract_region(gauss, sub_region_int)
    flux_line = np.trapz(sub_gauss_int.flux, sub_gauss_int.spectral_axis) 
    #Ploting the lina and fit Gaussian
    if saveplot == "y":
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
        plt.savefig(type_spec + "_" + NameLine + ".pdf")
        plt.close()
    else:
        pass
    return flux_line

def ew(wl_vacuum, spec):
    '''
    Estimating the equivalent width of a line
    '''
    line_spec_mask = find_line(wl_vacuum,  spec)
    min_lamb = line_spec_mask['line_center'] - 5 * u.AA
    max_lamb = line_spec_mask['line_center'] + 5 * u.AA
    sub_region_int = SpectralRegion(min_lamb,  max_lamb)
    sub_spect_int = extract_region(spec, sub_region_int)
    #Estimating equivalent width
    with warnings.catch_warnings():  # Ignore warnings
        warnings.simplefilter('ignore')
        cont_norm_spec = spec / fit_generic_continuum(spec)(spec.spectral_axis)
    ew = equivalent_width(cont_norm_spec, regions=sub_region_int)
    return ew
        
def npx_errcont(wl_vacuum, spec):
    '''
    Find mean standar deviation both side of the line, and number of pixel cover for the 
    line
    '''
    line_spec_mask = find_line(wl_vacuum,  spec)
    min_lamb = line_spec_mask['line_center'] - 5 * u.AA
    max_lamb = line_spec_mask['line_center'] + 5 * u.AA
    sub_region_int = SpectralRegion(min_lamb, max_lamb)
    sub_spect_int = extract_region(spec, sub_region_int)
    line_para_line = estimate_line_parameters(sub_spect_int, models.Gaussian1D())
    min_lamb_ = line_para_line.mean.value - 3*line_para_line.stddev.value
    max_lamb_ = line_para_line.mean.value + 3*line_para_line.stddev.value
    sub_region_line_ = SpectralRegion(min_lamb_ * u.AA,  max_lamb_ * u.AA)
    sub_line_ = extract_region(spec, sub_region_line_)
    n_pixel = len(sub_line_.spectral_axis)
    # Determinante the median desviation standar in both side of the line
    min_lamb_cont = min_lamb_ - 20
    max_lamb_cont = min_lamb_ + 20
    sub_region_cont_left = SpectralRegion(min_lamb_cont * u.AA, min_lamb_ * u.AA) # extract spec on left side of the line
    sub_region_cont_right = SpectralRegion(max_lamb_ * u.AA, max_lamb_cont * u.AA) # extract spec on right side of the line
    sub_cont_left =  extract_region(spec, sub_region_cont_left)
    sub_cont_right =  extract_region(spec, sub_region_cont_right)
    err = []
    for errleft in sub_cont_left.uncertainty.array:
        err.append(errleft)
    for errright in sub_cont_right.uncertainty.array:
        err.append(errright)
    err = np.array(err).mean()
    # equivalent widht
    return n_pixel, err

def err_line(wl_vacuum, spec, D=1.0002302850208247):
    '''
    Estimate the uncertainty of the lines
    using the eq. of Tresse et al. 1999
    '''
    ew1 = ew(wl_vacuum, spec) #equivalent width
    npixel, err_cont = npx_errcont(wl_vacuum, spec) #number of pixel and error continuum
    err_line = err_cont * D * np.sqrt((2 * npixel) + (np.abs(ew1.value) / D))
    return err_line
    
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

lines = {"[NeIII]H7": 3967.470,
         "Hdelta": 4101.742,
         "Hgamma": 4340.471,
         "HeII": 4685.99,
         "Hbeta": 4861.333,
         "[OIII]4958": 4958.911,
         "[OIII]5006": 5006.843,
         "[FeIII]": 5412.12,
         "Halpha": 6564.614,
         }

nlines_ = []
lines_ = []
flux_lines = []
err_lines = []
EW = []
for v, t in lines.items():
    flux_lines.append(measurent_line(t, spec, lamb, wl, Flux, rel_flux, "Obs", v, saveplot = "n").value)
    err_lines.append(err_line(t, spec, D = D))
    EW.append(ew(t, spec))
    lines_.append(t)
    nlines_.append(v)

Hbeta = flux_lines[5]
ratio_lines = flux_lines / Hbeta

# And the error
err_Hbeta = err_lines[5]
err_ratio_lines = np.sqrt((ratio_lines**2) * ((np.array(err_lines) / np.array(flux_lines))**2 + (err_Hbeta / Hbeta)**2))

#creating table and save it
table = QTable([nlines_, lines_,  ratio_lines, err_ratio_lines, EW],
           names=('Line', 'Lambda', 'Flux', 'Sigma', "EW"),
           meta={'name': 'first table'})
#save the table
table.write("parameters-lamost-pn-selec-lines.ecsv", format="ascii.ecsv", overwrite=True)

#################################################################################
# Model  ########################################################################
#################################################################################
lamb_model = data_mask["Wl"] * u.AA 
flux_model = data_mask["Flux"] *  u.Unit('erg cm-2 s-1 AA-1') 
spec_model = Spectrum1D(spectral_axis=lamb_model, flux=flux_model)

#Estimating the flux line for the models
flux_lines_models = []
for vv, tt in lines.items():
    flux_lines_models.append(measurent_line_model(tt, spec_model, lamb_model, data_mask["Wl"], data_mask["Flux"],  u.Unit('erg cm-2 s-1 AA-1'), "Model", vv, saveplot = "n").value)


Hbeta_models = flux_lines_models[5]
ratio_lines_models = flux_lines_models / Hbeta_models

# Estimating the chi-square
chi = (ratio_lines_models - ratio_lines)**2 / err_ratio_lines**2
chi_sum = chi.sum()

# Estimating the degree freedom
n = 9
np = 3
vv = n - np
chi_sum_red = chi_sum / vv

if chi_sum_red <= 2:
    tab = QTable([model_name, hi_sum_red],
           names=('Name model', 'Chi red'),
           meta={'name': 'first table'})
    tab.write("better-models-chisquera.ecsv", format="ascii.ecsv", overwrite=True)

