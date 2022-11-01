from astropy.io import ascii
from astropy.table import Table
import numpy as np
import argparse
from astropy.io import fits
import matplotlib.pyplot as plt
import seaborn as sn
import glob
sn.set_context("poster")


parser = argparse.ArgumentParser(
    description="""Making the LAMOST spectra""")

parser.add_argument("source", type=str,
                    default="spec-56568-VB081N11V1_sp05-101",
                    help="Name of source-spectrum, taken the prefix")

cmd_args = parser.parse_args()
file_ = cmd_args.source + ".fits"

hdu = fits.open(file_)
hdudata = hdu[0].data
wl = hdudata[2]
Flux = hdudata[0]

# Emission lines
wv_lin = [3869.76, 3967.46, 4101.74, 4340.471, 4542, 4685.71, 4711.26, 4740.120, 4861.33, 4958.92, 5006.8, 5411, 6562.82, 6879,  7005.87, 7113, 7135, 7751, 8236.79, 9001.27]
#em_lin = [r"[Ne III]" + "\n" + "3869", "[Ne III] + H7" + "\n" + "3968", r"H$\delta$", r"H$\gamma$", "He II", "He II", "[Ar IV]" + "\n" + "4711", "[Ar IV]" + "\n" + "4740", r"H$\beta$", "[O III]" + "\n" + "4958", "[O III]" + "\n" + "5007", "He II" + "\n" + "5411", "[N II]", r"H$\alpha$", "[N II]", "[S II]" + "\n" + "6731", "[Ar V]" + "\n" + "7005", "?" + "\n" + "7113", "[Ar III]" + "\n" + "7135", "[Ar III]" + "\n" + "7751", "He II" + "\n" + "8236", "[S III]" + "\n" + "9069"]

em_lin = [r"[Ne III] 3869", "[Ne III] + H7 3968", r"H$\delta$", r"H$\gamma$", "He II", "He II", "[Ar IV] 4711", "[Ar IV] 4740", r"H$\beta$", "[O III] 4958", "[O III] 5007", "He II 5411",  r"H$\alpha$", "? 6879", "[Ar V] 7005", "? 7113", "[Ar III] 7135", "[Ar III] 7751", "He II 8236", "? 9001"]

#m = (5000 < wl) &  (6000 > wl)
#new_flux = Flux[m].max()
#print(new_flux)

max_flux = []
for i in wv_lin:
    j = i - 10
    k = i + 10
    mask = (j < wl) & (wl < k)
    wl_= wl[mask]
    flux = Flux[mask]
    try:
        max_flux.append(np.max(flux))
    except ValueError:
        max_flux.append(10)
               
# PLOT
n = np.linspace(1, len(wl), num=len(wl), dtype = int)
fig, ax = plt.subplots(figsize=(11, 5))
#ax.set_title(namefile)
ax.set(xlim=[3600,9100])
#plt.ylim(ymin=0.0,ymax=500)
ax.set(xlabel='Wavelength $(\AA)$')
ax.set(ylabel='Flux')
ax.plot(wl, Flux, c = "blueviolet", linewidth=0.7, zorder=5)
for wll in wv_lin:
    ax.axvline(wll, color='k', linewidth=0.4, alpha=0.5, linestyle='--')
# ax.axvline(6560.28, color='k', linewidth=0.5, linestyle='--', label=r"H$\alpha$")
# ax.axvline(5000.7, color='k', linewidth=0.5, linestyle='--', label="[O III]")
# ax.axvline(4861.33, color='k', linewidth=0.5, linestyle='--', label=r"H$\beta$")
# ax.axvline(4686, color='k', linewidth=0.5, linestyle='--', label="He II")
# ax.axvline(7751, color='k', linewidth=0.4, linestyle='--', label="[Ar III] 7751")
# ax.axvline(7135, color='k', linewidth=0.4, linestyle='--', label="[Ar III] 7135")
# ax.axvline(9069, color='k', linewidth=0.4, linestyle='--', label="[S III] 9069")
# ax.axvline(4072, color='k', linewidth=0.4, linestyle='--', label="[S II] 4072")
# ax.axvline(6716.42, color='k', linewidth=0.4, linestyle='--', label="[S II] 6716")
# ax.axvline(4363.21, color='k', linewidth=0.4, linestyle='--', label="[O III] 4363")
# ax.axvline(4958.92, color='k', linewidth=0.4, linestyle='--', label="[O III] 4958")

bbox_props = dict(boxstyle="round", fc="w", ec="0.88", alpha=0.6, pad=0.1)
for label_, x, y in zip(em_lin, wv_lin, max_flux):
    ax.annotate(label_, (x, y), alpha=1, size=4,
                   xytext=(3.0, 5.6), textcoords='offset points', ha='right', va='bottom', rotation=90, bbox=bbox_props, zorder=100)
    
#ax.legend()
plt.tight_layout()
namefile = file_.replace(".fits", ".pdf")
plt.savefig(namefile)
ra = hdu[0].header["RA"]
dec = hdu[0].header["DEC"]
print("Nome file:", namefile, "=>", "RA and DEC:", ra, dec)

############################################################
#----------------------------------------------------------#
############################################################
fig, ax = plt.subplots(figsize=(11, 5))
#ax.set_title(namefile)
ax.set(xlim=[3600,9100])
plt.ylim(ymin=-200,ymax=1500)
ax.set(xlabel='Wavelength $(\AA)$')
ax.set(ylabel='Flux')
ax.plot(wl, Flux, c = "blueviolet", linewidth=0.7, zorder=5)
for wll in wv_lin:
    ax.axvline(wll, color='k', linewidth=0.4, alpha=0.5, linestyle='--')

bbox_props = dict(boxstyle="round", fc="w", ec="0.78", alpha=0.6, pad=0.1)
for label_, x, y in zip(em_lin, wv_lin, max_flux):
    ax.annotate(label_, (x, y), alpha=1, size=4,
                   xytext=(3.0, 5.6), textcoords='offset points', ha='right', va='bottom', rotation=90, bbox=bbox_props, zorder=100)
    
#ax.legend()
plt.tight_layout()
namefile0 = file_.replace(".fits", "-zoom.pdf")
plt.savefig(namefile0)


############################################################
#----------------------------------------------------------#
############################################################
fig1, ax1 = plt.subplots(figsize=(12, 10))
#ax.set_title(namefile)
ax1.set(xlim=[6560.28-100, 6560.28+100])
#plt.ylim(ymin=0.0,ymax=500)
ax1.set(xlabel='Wavelength $(\AA)$')
ax1.set(ylabel='Flux')
ax1.plot(wl, Flux, c = "blueviolet", linewidth=0.7, zorder=5)

    
#ax.legend()
plt.tight_layout()
namefile1 = file_.replace(".fits", "-Halpha.pdf")
plt.savefig(namefile1)
