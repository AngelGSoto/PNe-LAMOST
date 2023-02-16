from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
#from mpl_toolkits.axes_grid1.inset_locator import mark_inset, inset_axes
from astropy.io import ascii
from astropy.table import Table
import numpy as np
import argparse
from astropy.io import fits
import matplotlib.pyplot as plt
import seaborn as sn
import glob
sn.set_context("poster")


# Data of the spectra
hdulist0 = fits.open("known-pn/NGC2242/Kwitter_NGC2242A_KP091296_id704.fits")
hdu0 = hdulist0[0]
nx0, wav00, i00, dwav0 = [hdu0.header[k] for k in ("NAXIS1", "CRVAL1", "CRPIX1", "CDELT1")]
wl0 = wav00 + (np.arange(nx0) - (i00 - 1))*dwav0

Flux0 = hdulist0[0].data
Flux0 /=10e-14
Flux0 +=0.8

hdulist1 = fits.open("known-pn/NGC4361/StenholmAcker_pn_g294_1+43_6_id910.fits")
hdu1 = hdulist1[0]
nx1, wav01, i01, dwav1 = [hdu1.header[k] for k in ("NAXIS1", "CRVAL1", "CRPIX1", "CDELT1")]
wl1 = wav01 + (np.arange(nx1) - (i01 - 1))*dwav1

Flux1 = hdulist1[0].data
Flux1 /=10e2
Flux1 +=1.6

hdulist2 = fits.open("known-pn/PNPRTM1/SAAO2016_PrTm1_SA061116_id797.fits")
hdu2 = hdulist2[0]
nx2, wav02, i02, dwav2 = [hdu2.header[k] for k in ("NAXIS1", "CRVAL1", "CRPIX1", "CDELT1")]
wl2 = wav02 + (np.arange(nx2) - (i02 - 1))*dwav2

Flux2 = hdulist2[0].data
Flux2 /=10e-14 
Flux2 += 2.3
# New object
hdu = fits.open("Spectra-lamostdr7/spec-56581-VB031N50V1_sp08-218.fits")
hdudata = hdu[0].data
wl = hdudata[2]
Flux = hdudata[0]
Flux /=10e3

#spliting the spectra in the blue and red arms
m_blue = (wl >= min(wl)) & (wl <= 5800)
m_red = (wl >= 6300) & (wl <= max(wl))
wl_blue = wl[m_blue]
wl_red = wl[m_red]
Flux_blue = Flux[m_blue]
Flux_red = Flux[m_red]

# Emission lines
wv_lin = [3869.76, 3967.46, 4101.74, 4340.471, 4542, 4685.71, 4711.26, 4740.120, 4861.33, 4958.92, 5006.8, 5411, 6562.82, 6879,  7005.87, 7113, 7135, 7751, 8236.79, 9001.27]
#em_lin = [r"[Ne III]" + "\n" + "3869", "[Ne III] + H7" + "\n" + "3968", r"H$\delta$", r"H$\gamma$", "He II", "He II", "[Ar IV]" + "\n" + "4711", "[Ar IV]" + "\n" + "4740", r"H$\beta$", "[O III]" + "\n" + "4958", "[O III]" + "\n" + "5007", "He II" + "\n" + "5411", "[N II]", r"H$\alpha$", "[N II]", "[S II]" + "\n" + "6731", "[Ar V]" + "\n" + "7005", "?" + "\n" + "7113", "[Ar III]" + "\n" + "7135", "[Ar III]" + "\n" + "7751", "He II" + "\n" + "8236", "[S III]" + "\n" + "9069"]

em_lin = [r"[Ne III] 3869", "[Ne III] + H7 3968", r"H$\delta$", r"H$\gamma$", "He II", "He II", "[Ar IV] 4711", "[Ar IV] 4740", r"H$\beta$", "[O III] 4958", "[O III] 5007", "He II 5411",  r"H$\alpha$", "? 6879", "[Ar V] 7005", "? 7113", "[Ar III] 7135", "[Ar III] 7751", "He II 8236", "? 9001"]

#m = (5000 < wl) &  (6000 > wl)
#new_flux = Flux[m].max()
#print(new_flux)

max_flux = [2.74, 2.74, 2.74, 2.74, 2.74, 2.74, 2.74, 2.74, 2.74, 2.74, 2.74, 2.74, 2.74, 2.74, 2.74, 2.74, 2.74, 2.74, 2.74, 2.74]
               
# PLOT
n = np.linspace(1, len(wl), num=len(wl), dtype = int)
fig, ax = plt.subplots(figsize=(12, 12))
#ax.set_title(namefile)
ax.set(xlim=[3600,9100])
ax.set(ylim=[-0.07,3.0])
#plt.ylim(ymin=0.0,ymax=500)
ax.set(xlabel='Wavelength $(\AA)$')
ax.set(ylabel='Normalised flux')
ax.plot(wl2, Flux2 , c = "red", linewidth=2., zorder=5)
ax.plot(wl1, Flux1 , c = "#FFA500", linewidth=2., zorder=5)
ax.plot(wl0, Flux0 , c = "green", linewidth=2., zorder=5)
ax.plot(wl, Flux , c = "white", linewidth=2., zorder=5)
ax.plot(wl_blue, Flux_blue, c = "blueviolet", linewidth=2, zorder=5)
ax.plot(wl_red, Flux_red, c = "blueviolet", linewidth=2, zorder=5)
for wll in wv_lin:
    ax.axvline(wll, color='k', linewidth=0.4, alpha=0.5, linestyle='--')

#ax.get_yaxis().set_visible(True)
plt.yticks([])
bbox_props = dict(boxstyle="round", fc="w", ec="0.88", alpha=0.6, pad=0.1)
for label_, x, y in zip(em_lin, wv_lin, max_flux):
    ax.annotate(label_, (x, y), alpha=1, size=6,
                   xytext=(3.0, 5.6), textcoords='offset points', ha='right', va='bottom', rotation=90, bbox=bbox_props, zorder=200)

plt.text(0.78, 0.1, 'J020808.63+491401.0',
         transform=ax.transAxes, c="blueviolet", weight='bold', fontsize=12.8, bbox=bbox_props)

plt.text(0.85, 0.35, 'NGC 2242',
         transform=ax.transAxes, c="g", weight='bold', fontsize=12.8, bbox=bbox_props)

plt.text(0.85, 0.6, 'NGC 4361',
         transform=ax.transAxes, c="orange", weight='bold', fontsize=12.8, bbox=bbox_props)

plt.text(0.85, 0.85, 'PN PRTM 1',
         transform=ax.transAxes, c="red", weight='bold', fontsize=12.8, bbox=bbox_props)
#ax.legend()
plt.tight_layout()
plt.savefig("Text/Figs/spectra-compare.pdf")
