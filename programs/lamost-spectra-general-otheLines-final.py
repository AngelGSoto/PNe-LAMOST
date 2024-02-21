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

parser.add_argument("--name", type=str,
                    default="PSP",
                    help="Name of the objet")

parser.add_argument("--ymin", required=False, type=float, default=None,
                    help="""Value y-axis min""")

parser.add_argument("--ymax", required=False, type=float, default=None,
                    help="""Value y-axis max""")

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
Flux /=10e3

#spliting the spectra in the blue and red arms
m_blue = (wl >= min(wl)) & (wl <= 5800)
m_red = (wl >= 5800) & (wl <= max(wl))
wl_blue = wl[m_blue]
wl_red = wl[m_red]
Flux_blue = Flux[m_blue]
Flux_red = Flux[m_red]

# Emission lines
emission_lines = {
    "[Ne III] 3869": 3868.75,
    "[Ne III] + H7 3968": 3967.46,
    "Hδ": 4101.74,
    "Hγ": 4340.471,
    "[O III] 4363": 4363.21,
    "He I 4472": 4471.50 ,
    "O II 4639": 4638.86,
    "He II 4686": 4685.68,
    "[Ar IV] 4711": 4711.26,
    "[Ar IV] 4740": 4740.12,
    "Hβ": 4861.33,
    "He I 4922": 4921.93,
    "[O III] 4958": 4958.911,
    "[O III] 5007": 5006.843,
    "He II 5412": 5411.52,
    "He I 5876": 5875.624,
    "He II 6311": 6310.85,
    "[N II] 6548": 6548.05,
    "Hα": 6562.82,
    "[N II] 6583": 6583.46,
    "He II 6683": 6683.20,
    "[S II] 6731": 6730.82,
    "He I 7065": 7065.25,
    "[Ar III] 7135": 7135,
    "He I 7281": 7281.35,
    "[O II] 7319": 7318.92,
    "[Ar III] 7751": 7751,
    "[Cl IV] 8046 ": 8045.63 ,
    "He II 8236": 8236.79,
    "[S III] 9069": 9069
}

emission_lines1 = {
    "[Ne III] 3869": 3868.75,
    "[Ne III] + H7 3968": 3967.46,
    "Hδ": 4101.74,
    "Hγ": 4340.471,
    "[O III] 4363": 4363.21,
    "He I 4472": 4471.50 ,
    "O II 4639": 4638.86,
    "He II 4686": 4685.68,
    "[Ar IV] 4711": 4711.26,
    "[Ar IV] 4740": 4740.12,
    "Hβ": 4861.33,
    "He I 4922": 4921.93,
    "[O III] 4958": 4958.911,
    "[O III] 5007": 5006.843,
    "He II 5412": 5411.52,
    "C IV 5801": 5801.33,
    "He I 5876": 5875.624,
    "He II 6311": 6310.85,
    "[N II] 6548": 6548.05,
    "Hα": 6562.82,
    "[N II] 6583": 6583.46,
    "He II 6683": 6683.20,
    "[S II] 6731": 6730.82,
    "He I 7065": 7065.25,
    "[Ar III] 7135": 7135,
    "He I 7281": 7281.35,
    "[O II] 7319": 7318.92,
    "[Ar III] 7751": 7751,
    "[Cl IV] 8046 ": 8045.63 ,
    "He II 8236": 8236.79,
    "[S III] 9069": 9069
}

displace_lines = {
    "He II 4686": 1.5, 
    "[Ar IV] 4711": 0.8,
    "Hα": 0.8,
    "He II 6683": 0.5,
    "[N II] 6583": 0.8,
    "[O II] 7319": 0.8,
}

# Lines to displace to the middle of the line
lines_to_adjust = ["[O III] 4958", "[O III] 5007"]
##########################################################################################
displace_lines1 = {
    "[Ar IV] 4711": 4,
    "Hα": 2,
    "He II 6683": 4,
    "[N II] 6583": 4,
    "[O II] 7319": 0.8,
    "[O II] 7319": 3
}

# Lines to displace to the middle of the line
lines_to_adjust1 = ["[Ne III] 3869", "Hβ", "[O III] 4958", "He II 4686"]


#############################################33
max_flux = []
for wavelength in emission_lines1.values():
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
ax.set(ylim=[-0.05,0.7])
plt.tick_params(axis='x', labelsize=16) 
plt.tick_params(axis='y', labelsize=16)
#plt.ylim(ymin=0.0,ymax=500)
# ax.set(xlabel='Wavelength $(\AA)$')
# ax.set(ylabel='Relative flux')
ax.set_xlabel('Wavelength $(\AA)$', fontsize=16)
ax.set_ylabel('Relative flux', fontsize=16)
ax.plot(wl, Flux , c = "blueviolet", linewidth=1.3, zorder=5)

for label, wavelength in emission_lines1.items():
    wll_ob = wavelength 
    max_flux_val = max_flux[list(emission_lines1.keys()).index(label)]
    bbox_props = dict(boxstyle="round", fc="w", ec="0.88", alpha=0.6, pad=0.1)

    # Check if label should be displaced
    if displace_lines1 and label in displace_lines1:
        displacement_factor = displace_lines1[label]
        # Adjust y-position based on displacement factor
        max_flux_val += displacement_factor
    elif label in lines_to_adjust1:
        flux_values = Flux[(wl >= wll_ob - 10) & (wl <= wll_ob + 10)]
        y_position = (np.max(flux_values) + np.min(flux_values)) / 1.9
        ax.annotate(label, (wll_ob, y_position), alpha=1, size=10,
                    xytext=(4.8, 5.6), textcoords='offset points', ha='right', va='bottom',
                    rotation=90, bbox=bbox_props, zorder=200)
        continue

    ax.axvline(wll_ob, color='k', linewidth=0.8, alpha=0.5, linestyle='--')
    ax.annotate(label, (wll_ob, max_flux_val), alpha=1, size=10,
                xytext=(4.8, 5.6), textcoords='offset points', ha='right', va='bottom',
                rotation=90, bbox=bbox_props, zorder=200)



# Add text to the plot
plt.text(5500, cmd_args.ymax/1.2, cmd_args.name, fontsize=19, color='black',
        bbox=dict(facecolor='lightgray', alpha=0.5, edgecolor='gray', boxstyle='round'),
        horizontalalignment='center', verticalalignment='center')

# set Y-axis range (if applicable)
if cmd_args.ymin is not None and cmd_args.ymax is not None:
    plt.ylim(cmd_args.ymin, cmd_args.ymax)
elif cmd_args.ymin is not None:
    plt.ylim(ymin=cmd_args.ymin)
elif cmd_args.ymax is not None:
    plt.ylim(ymax=cmd_args.ymax)

#plt.xticks(visible=False)  # Not present ticks
    
#ax.legend()
plt.tight_layout()
namefile = file_.replace(".fits", ".pdf")
plt.savefig(namefile)
ra = hdu[0].header["RA"]
dec = hdu[0].header["DEC"]
print("Nome file:", namefile, "=>", "RA and DEC:", ra, dec)
