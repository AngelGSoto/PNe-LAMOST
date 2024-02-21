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
from matplotlib.patches import Ellipse
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
#sn.set_context("poster")

parser = argparse.ArgumentParser(
    description="""Making the LAMOST spectra""")

parser.add_argument("source", type=str,
                    default="spec-56568-VB081N11V1_sp05-101",
                    help="Name of source-spectrum, taken the prefix")

cmd_args = parser.parse_args()
file_ = cmd_args.source + ".fits"

#LAMOST DR8
hdu = fits.open(file_)
hdudata = hdu[0].data
wl = hdudata[2]
Flux = hdudata[0]

# DR8 LAMOST
#wl = hdudata["WAVELENGTH"][0]
#Flux =hdudata["FLUX"][0]
Flux /=1e3


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
        "C IV 5802": 5801.51,
        "Hα": 6562.82,
        "? 6879": 6879,
        "[Ar V] 7006": 7005.67,
        "? 7113": 7113,
        "[Ar III] 7536": 7135.80,
        "He II 8237": 8236.79,
        "? 9001":  9000.99,
       
    }

displace_lines = {
    "He II 4686": 0.1, 
    "[Ar IV] 4711": 1.2,
    "? 7113": 1.0 # Apply a displacement factor of 1.5 to the Hα line  # Apply a displacement factor of 3.0 to the OIII line
}


max_flux = []
for wavelength in emission_lines.values():
    lambda_ob = wavelength 
    j = lambda_ob - 10
    k = lambda_ob + 10
    mask = (j < wl) & (wl < k)
    flux_values = Flux[mask]
    try:
        max_flux.append(np.max(flux_values))
    except ValueError:
        max_flux.append(10)
        
# PLOT
n = np.linspace(1, len(wl), num=len(wl), dtype = int)
fig, ax = plt.subplots(figsize=(11, 5))
#ax.set_title(namefile)
ax.set(xlim=[3600,9100])
ax.set(ylim=[-0.5,5])
plt.tick_params(axis='x', labelsize=16) 
plt.tick_params(axis='y', labelsize=16)
#plt.ylim(ymin=0.0,ymax=500)
# ax.set(xlabel='Wavelength $(\AA)$')
# ax.set(ylabel='Relative flux')
ax.set_xlabel('Wavelength $(\AA)$', fontsize=16)
ax.set_ylabel('Relative flux', fontsize=16)
ax.plot(wl, Flux , c = "blueviolet", linewidth=0.7, zorder=5)

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



#sn.despine()
###########
#zoom plot#
###########
axins = zoomed_inset_axes(ax, 2.2, loc=2, bbox_to_anchor=(0.39, 0.94),
                   bbox_transform=ax.figure.transFigure) # zoom = 6
axins.plot(wl, Flux, c = "blueviolet", linewidth=0.7, zorder=5)
axins.set_xlim(4575, 4780) # Limit the region for zoom
axins.set_ylim(-0.05, 1)


for (label_, x), y in zip(emission_lines.items(), max_flux):
    wll_ob1 = x 
    axins.annotate(label_, (wll_ob1, y), alpha=1, size=8.2,
                   xytext=(4.8, 5.6), textcoords='offset points', ha='right', va='bottom', rotation=90, bbox=bbox_props, zorder=200)

plt.xticks(visible=False)  # Not present ticks
plt.yticks(visible=False)
#
## draw a bbox of the region of the inset axes in the parent axes and
## connecting lines between the bbox and the inset axes area
mark_inset(ax, axins, loc1=2, loc2=4, fc="none", lw=1.2,  ec="0.6", zorder=1)

###########
#Zoom other region
axins1 = zoomed_inset_axes(ax, 2.2, loc=2, bbox_to_anchor=(0.73, 0.94),
                   bbox_transform=ax.figure.transFigure) # zoom = 6
axins1.plot(wl, Flux, c = "blueviolet", linewidth=0.7, zorder=5)
axins1.set_xlim(6830, 7190) # Limit the region for zoom
axins1.set_ylim(-0.05, 1)

for (label_, x), y in zip(emission_lines.items(), max_flux):
    wll_ob2 = x 
    axins1.annotate(label_, (wll_ob2, y), alpha=1, size=8.2,
                    xytext=(4.8, 5.6), textcoords='offset points', ha='right', va='bottom', rotation=90, bbox=bbox_props, zorder=200)


plt.xticks(visible=False)  # Not present ticks
plt.yticks(visible=False)
#
## draw a bbox of the region of the inset axes in the parent axes and
## connecting lines between the bbox and the inset axes area
mark_inset(ax, axins1, loc1=2, loc2=4, fc="none",lw=1.2,  ec="0.6", zorder=1)

# axins = inset_axes(ax, 1, 1, loc=2, bbox_to_anchor=(0.2, 0.55),
#                   bbox_transform=ax.figure.transFigure)
# axins.plot(x, y)
# x1, x2 = .4, .6
# y1, y2 = x1 ** 2, x2 ** 2

# axins.set_xlim(x1, x2)
# axins.set_ylim(y1, y2)

# mark_inset(ax, axins, loc1=1, loc2=3, fc="none", ec="0.5")


#ax.legend()
plt.tight_layout()
namefile = file_.replace(".fits", ".pdf")
plt.savefig("gr4.pdf")
ra = hdu[0].header["RA"]
dec = hdu[0].header["DEC"]
print("Nome file:", namefile, "=>", "RA and DEC:", ra, dec)

