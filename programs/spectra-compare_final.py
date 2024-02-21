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
#sn.set_context("poster")

# Data of the spectra
hdulist0 = fits.open("known-pn/NGC2242/Kwitter_NGC2242A_KP091296_id704.fits")
hdu0 = hdulist0[0]
nx0, wav00, i00, dwav0 = [hdu0.header[k] for k in ("NAXIS1", "CRVAL1", "CRPIX1", "CDELT1")]
wl0 = wav00 + (np.arange(nx0) - (i00 - 1))*dwav0

Flux0 = hdulist0[0].data
Flux0 /=10e-14
Flux0 +=2.3

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
Flux2 +=  0.8
# New object
hdu = fits.open("Spectra-lamostdr7/spec-56581-VB031N50V1_sp08-218.fits")
hdudata = hdu[0].data
wl = hdudata[2]
Flux = hdudata[0]
Flux /=10e3

# Emission lines
emission_lines = {
        "[Ne III] 3869": 3868.75,
        "[Ne III] + H7 3968": 3967.46,
        "Hδ": 4101.74,
        "Hγ": 4340.471,
        "He II 4542": 4541.59,
        "He II 4686": 4685.68,
        "[Ar IV] 4711": 4711.37,
        "[Ar IV] 4740": 4740.17,
        "Hβ": 4861.33,
        "[O III] 4958": 4958.911,
        "[O III] 5007": 5006.843,
        "He II 5412": 5411.52,
        "Hα": 6562.82,
        "? 6879": 6879,
        "[Ar V] 7006": 7005.67,
        "? 7113": 7113,
        "[Ar III] 7536": 7135.80,
        "He II 8237": 8236.79,
        "? 9001":  9000.99,  
    }

displace_lines = {
    "He II 4686": 0.15, 
    "[Ar IV] 4711": 0.26,
    "? 7113": 0.27 # Apply a displacement factor of 1.5 to the Hα line  # Apply a displacement factor of 3.0 to the OIII line
                 }

max_flux = []
for wavelength in emission_lines.values():
    lambda_ob = wavelength 
    j = lambda_ob - 10
    k = lambda_ob + 10
    mask = (j < wl0) & (wl0 < k)
    flux_values = Flux0[mask]
    try:
        max_flux.append(np.max(flux_values))
    except ValueError:
        max_flux.append(10)
               
# PLOT
n = np.linspace(1, len(wl), num=len(wl), dtype = int)
fig, ax = plt.subplots(figsize=(12, 12))
#ax.set_title(namefile)
ax.set(xlim=[3600,9100])
ax.set(ylim=[-0.07,3.5])
#plt.ylim(ymin=0.0,ymax=500)
plt.tick_params(axis='x', labelsize=18) 
plt.tick_params(axis='y', labelsize=18)
ax.set_xlabel('Wavelength $(\AA)$', fontsize=18)
ax.set_ylabel('Normalised flux', fontsize=18)
# ax.set(xlabel='Wavelength $(\AA)$')
# ax.set(ylabel='Normalised flux')
ax.plot(wl2, Flux2 , c = "red", linewidth=2., zorder=5)
ax.plot(wl1, Flux1 , c = "#FFA500", linewidth=2., zorder=5)
ax.plot(wl0, Flux0 , c = "green", linewidth=2., zorder=5)
ax.plot(wl, Flux , c = "blueviolet", linewidth=2., zorder=5)


for label, wavelength in emission_lines.items():
    wll_ob = wavelength 
    max_flux_val = max_flux[list(emission_lines.keys()).index(label)]
    bbox_props = dict(boxstyle="round", fc="w", ec="0.88", alpha=0.6, pad=0.1)

    # Check if label should be displaced
    if displace_lines and label in displace_lines:
        displacement_factor = displace_lines[label]
        # Adjust y-position based on displacement factor
        max_flux_val += displacement_factor
    else:
        displacement_factor = 2  # Default displacement factor

    ax.axvline(wll_ob, color='k', linewidth=0.8, alpha=0.5, linestyle='--')
    ax.annotate(label, (wll_ob, max_flux_val), alpha=1, size=8.2,
                xytext=(4.8, 5.6), textcoords='offset points', ha='right', va='bottom',
                rotation=90, bbox=bbox_props, zorder=200)

plt.yticks([])
plt.text(0.78, 0.1, 'J020808.63+491401.0',
         transform=ax.transAxes, c="blueviolet", weight='bold', fontsize=12.8, bbox=bbox_props)

plt.text(0.85, 0.85, 'NGC 2242',
         transform=ax.transAxes, c="g", weight='bold', fontsize=12.8, bbox=bbox_props)

plt.text(0.85, 0.6, 'NGC 4361',
         transform=ax.transAxes, c="orange", weight='bold', fontsize=12.8, bbox=bbox_props)

plt.text(0.85, 0.35, 'PN PRTM 1',
         transform=ax.transAxes, c="red", weight='bold', fontsize=12.8, bbox=bbox_props)
#ax.legend()
plt.tight_layout()
plt.savefig("Text/Figs/gr6.pdf")
