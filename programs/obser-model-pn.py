import numpy as np
import matplotlib.pyplot as plt
from  astropy.table import Table, vstack, hstack
import pandas as pd
import argparse
from astropy.io import fits
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=FutureWarning)
    import h5py


parser = argparse.ArgumentParser(
    description="""Estimate magnitudes of a spectrum for any photometryc system""")

parser.add_argument("source", type=str,
                    default="DdDm-1",
                    help="Name of source, taken the prefix ")

parser.add_argument("--debug", action="store_true",
                    help="Print out verbose debugging info about each line in region file")

# read de file
# Model
#dat = Table.read("DdDm1_L4_T200_output_SED.dat.tere_E0.2", format="ascii")

cmd_args = parser.parse_args()
asciifile = cmd_args.source + ".dat.tere_E0.2"

dat = Table.read(asciifile, format="ascii")

wl = dat["col1"]
flux = dat["col2"]

# Observed
hdu = fits.open("../../Spectra-lamostdr7/spec-56581-VB031N50V1_sp08-218.fits")
hdudata = hdu[0].data
wl_o = hdudata[2]
Flux_o = hdudata[0]

# Matching the two spectra
mask_m = (wl > 5990) & (wl < 6010)
wl_m = wl[mask_m]
flux_m = flux[mask_m].mean() 

mask_o = (wl_o > 5990) & (wl_o < 6010)
wl_o_m = wl_o[mask_o].mean()
flux_o_m = Flux_o[mask_o].mean()

# Finding the factor
factor = flux_m / flux_o_m

###################################################################
###################################################################
###################################################################
fig, ax = plt.subplots(figsize=(11, 5))
#ax.set_title(namefile)
ax.set(xlim=[3600,9100])
#plt.ylim(ymin=0.0,ymax=500)
ax.set(xlabel='Wavelength $(\AA)$')
ax.set(ylabel='Flux')
ax.plot(wl, flux, c = "darkolivegreen", linewidth=0.7, zorder=5, label="Model")
ax.plot(wl_o, Flux_o*factor, c = "blueviolet", linewidth=0.7, zorder=5, label="Obs")

ax.legend()
plt.tight_layout()

if cmd_args.debug:
    print("Comparing with:", asciifile.split(".dat")[0])


# Creates The JSON files with the magnitudes
pdffile = asciifile.replace(".dat.tere_E0.2", 
                  "-E02-comparing-spectra.jpg")
plt.savefig(pdffile)
