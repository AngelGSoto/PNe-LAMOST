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
file_ = cmd_args.source + ".dat"

# Reading the Cloudy outputs in the Mod CloudyModel object
spec = Table.read(file_, format="ascii")


#Getting wl and flux
wl = spec["col1"]
flux = spec["col2"]
#Ordered lambda and flux 
wll, flux = zip(*sorted(zip(wl, flux)))

data = Table([wll, flux], names=('Wl', 'Flux'), meta={'name': 'first table'})
mask = (data["Wl"] > 3000) & (data["Flux"] < 9000)
data_mask = data[mask]

# OUR PN
hdu = fits.open("../../../../Spectra-lamostdr7/spec-56581-VB031N50V1_sp08-218.fits")
hdudata = hdu[0].data
wl = hdudata[2]
Flux = hdudata[0]


def closest(lst, K):
    '''find the closest number'''
    lst = np.array(lst)
    idx = (np.abs(lst - K)).argmin()
    return lst[idx]

# Model
MaskHbeta = (data_mask["Wl"] >= 4863 - 50) & (data_mask["Wl"] <= 4863 + 50)
HBeta = data_mask[MaskHbeta]
max_HBeta = HBeta["Flux"].max()
flux_m = data_mask["Flux"] / max_HBeta

# Our PN
MaskHbeta_our = (wl >= 4863 - 50) & (wl <= 4863 + 50)
HBeta_our = Flux[MaskHbeta_our]
max_HBeta_our = HBeta_our.max()
flux_our = Flux / max_HBeta_our

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
print("Spectra of:", file_)
asciifile = file_.replace(".dat", ".jpg") 
plt.savefig(asciifile)
