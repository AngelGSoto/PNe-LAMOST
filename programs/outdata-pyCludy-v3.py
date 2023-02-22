'''
Script to dealing with output data from Cloudy
Based in pyCloudy (Morisset, C., 2013, pyCloudy, 
Astrophysics Source Code Library)
Author: Luis A. Gutiérrez Soto
10/12/2022
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
flux = Mod.get_cont_y(cont = 'total', unit = 'esAc')
#Ordered lambda and flux 
wll, flux = zip(*sorted(zip(wl, flux)))

data = Table([wll, flux], names=('Wl', 'Flux'), meta={'name': 'first table'})
mask = (data["Wl"] > 3000) & (data["Wl"] < 9000)
data_mask = data[mask]

# OUR PN
hdu = fits.open("../../../Spectra-lamostdr7/spec-56581-VB031N50V1_sp08-218.fits")
hdudata = hdu[0].data
wl = hdudata[2]
Flux = hdudata[0]


def closest(lst, K):
    '''find the closest number'''
    lst = np.array(lst)
    idx = (np.abs(lst - K)).argmin()
    return lst[idx]

# Model
MaskHbeta = (data_mask["Wl"] >= (4863 - 50)) & (data_mask["Wl"] <= (4863 + 50))
HBeta = data_mask[MaskHbeta]
max_HBeta = HBeta["Flux"].max()
flux_m = data_mask["Flux"] / max_HBeta


# Our PN
MaskHbeta_our = (wl >= (4863 - 50)) & (wl <= (4863 + 50))
HBeta_our = Flux[MaskHbeta_our]
max_HBeta_our = HBeta_our.max()
flux_our = Flux / max_HBeta_our

fig, ax = plt.subplots(figsize=(11, 5))
#ax.set_title(namefile)
ax.set(xlim=[3600,9100])
#plt.ylim(ymin=-200,ymax=1500)
ax.set(xlabel='Wavelength $(\AA)$')
ax.set(ylabel='Normalised flux')
plt.plot(data_mask["Wl"], flux_m,  c = "darkolivegreen", linewidth=0.7, label = 'Model 3')
plt.plot(wl, flux_our, c = "blueviolet", linewidth=0.7, label = 'J020808.63+491401.0')

bbox_props = dict(boxstyle="round", fc="w", ec="0.88", alpha=0.6, pad=0.1)

info_mod1 = [r"$\mathrm{T_{eff}=16\times10^{4}K}$", r"$\mathrm{L=1.6\times10^{3}L_{\odot}}$", r"$\chi^2=6.78$", r"$\chi^2_r=1.13$"]
info_mod2 = [r"$\mathrm{T_{eff}=14\times10^{4}K}$", r"$\mathrm{L=1.9\times10^{3}L_{\odot}}$", r"$\chi^2=9.55$", r"$\chi^2_r=1.59$"]
info_mod3 = [r"$\mathrm{T_{eff}=13\times10^{4}K}$", r"$\mathrm{L=2.2\times10^{3}L_{\odot}}$", r"$\chi^2=9.67$", r"$\chi^2_r=1.61$"]
loc_text = [(0.75, 0.9), (0.75, 0.9-0.13), (0.75,  0.9-2*0.13), (0.75, 0.9-3*0.13)]
for loc_, taxt in zip(loc_text, info_mod3):
    ax.text(loc_[0], loc_[1], taxt, fontsize=19,
            bbox=dict(facecolor='gray', alpha=0.0),
            transform=ax.transAxes)
ax.legend(loc="upper left", fontsize="x-small")
sn.despine()
plt.tight_layout()
plt.savefig(model_name + ".pdf")
