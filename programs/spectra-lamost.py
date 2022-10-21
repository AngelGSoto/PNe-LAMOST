from astropy.io import ascii
from astropy.table import Table
import numpy as np
import argparse
from astropy.io import fits
import matplotlib.pyplot as plt
import seaborn as sn
sn.set_context("poster")
import glob

parser = argparse.ArgumentParser(
    description="""Making the LOMOST spectra""")

parser.add_argument("source", type=str,
                    default="spec-56568-VB081N11V1_sp05-101",
                    help="Name of source-spectrum, taken the prefix ")


cmd_args = parser.parse_args()
file_ = cmd_args.source + ".fits"

hdu = fits.open(file_)
hdudata = hdu[0].data
wl = hdudata[2]
Flux = hdudata[0]

n = np.linspace(1, len(wl), num=len(wl), dtype = int)
fig, ax = plt.subplots(figsize=(11, 5))
#ax.set_title(namefile)
ax.set(xlim=[3600,9100])
#plt.ylim(ymin=0.0,ymax=500)
ax.set(xlabel='Wavelength $(\AA)$')
ax.set(ylabel='Flux')
ax.plot(wl, Flux, linewidth=0.9, zorder=5)
ax.axvline(6560.28, color='k', linewidth=0.5, linestyle='--', label=r"H$\alpha$")
ax.axvline(5000.7, color='k', linewidth=0.5, linestyle='-', label="[O III]")
ax.axvline(4861.33, color='g', linewidth=0.5, linestyle='-', label=r"H$\beta$")
ax.axvline(4686, color='r', linewidth=0.5, linestyle='-', label="He II")
ax.axvline(7751, color='y', linewidth=0.4, linestyle='-.', label="[Ar III] 7751")
ax.axvline(7135, color='k', linewidth=0.4, linestyle='-.', label="[Ar III] 7135")
ax.axvline(9069, color='b', linewidth=0.4, linestyle='-.', label="[S III] 9069")
ax.axvline(4072, color='b', linewidth=0.4, linestyle='-.', label="[S II] 4072")
ax.axvline(6716.42, color='b', linewidth=0.4, linestyle='-.', label="[S II] 6716")
ax.axvline(4363.21, color='r', linewidth=0.4, linestyle='-.', label="[O III] 4363")
ax.axvline(4958.92, color='r', linewidth=0.4, linestyle='-.', label="[O III] 4958")
ax.legend()
namefile = file_.replace(".fits", ".pdf")
plt.savefig(namefile)
ra = hdu[0].header["RA"]
dec = hdu[0].header["DEC"]
print("Nome file:", namefile, "=>", "RA and DEC:", ra, dec)
