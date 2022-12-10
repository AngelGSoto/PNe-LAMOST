'''
Script to dealing with output data from Cloudy
Based in pyCloudy (Morisset, C., 2013, pyCloudy, Astrophysics Source Code Library)
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
mask = (data["Wl"] > 3000) & (data["Flux"] < 9000)
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
wl_hbeta = closest(data_mask["Wl"], 4861.333)
MaskHbeta = data_mask["Wl"] == wl_hbeta
HBeta = data_mask[MaskHbeta]
flux_m = data_mask["Flux"] / HBeta["Flux"]


# Our PN
wl_hbeta_our = closest(wl, 4861.333)

MaskHbeta_our = wl == wl_hbeta_our
flux_HBeta_our = Flux[MaskHbeta_our]
flux_our = Flux / flux_HBeta_our


fig, ax = plt.subplots(figsize=(11, 5))
#ax.set_title(namefile)
ax.set(xlim=[3600,9100])
#plt.ylim(ymin=-200,ymax=1500)
ax.set(xlabel='Wavelength $(\AA)$')
ax.set(ylabel='Flux')
plt.plot(data_mask["Wl"], flux_m,  c = "darkolivegreen", linewidth=0.7, label = 'Model')
plt.plot(wl, flux_our, c = "blueviolet", linewidth=0.7, label = 'Our PN')

ax.legend(loc="upper right")
sn.despine()
plt.tight_layout()
plt.savefig(model_name + ".jpg")


