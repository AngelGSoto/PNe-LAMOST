import __future__
import numpy as np
import json
import matplotlib.pyplot as plt
from  astropy.table import Table, vstack, hstack
import pandas as pd
from matplotlib.colors import PowerNorm
import seaborn as sns
from scipy.optimize import fsolve
from pathlib import Path
ROOT_PATH = Path("data")


#Find the point inteception between two lines     
def findIntersection(m, y, m1, y1, x0):
    x = np.linspace(-10.0, 15.5, 200)
    return fsolve(lambda x : (m*x + y) - (m1*x + y1), x0)

# Planetary nebula
tab_pn = Table.read("PN-PS-Gaia.ecsv", format="ascii.ecsv")
tab_star = Table.read("Star-PS-Gaia.ecsv", format="ascii.ecsv")
tab_star2 = Table.read("Star2-PS-Gaia.ecsv", format="ascii.ecsv")

catalogs = {'SySt': {'file': 'SySt-PS-GaiaEDR3.ecsv',
                     'G': 'phot_g_mean_mag',
                     'BP-RP': 'bp_rp',
                     'r': 'rmag'},
            'CV': {'file': 'CV-PS-GaiaEDR3.ecsv',
                    'G': 'phot_g_mean_mag',
                   'BP-RP': 'bp_rp',
                   'r': 'rmag'},
            'SNR': {'file': 'SNR-PS-GaiaEDR3.ecsv',
                    'G': 'phot_g_mean_mag',
                    'BP-RP': 'bp_rp',
                    'r': 'rmag'},
            'YSO': {'file': 'YSO-PS-GaiaEDR3.ecsv',
                    'G': 'phot_g_mean_mag',
                    'BP-RP': 'bp_rp',
                    'r': 'rmag'},
            'AeBe': {'file': 'AeBe-PS-GaiaEDR3.ecsv',
                     'G': 'phot_g_mean_mag',
                     'BP-RP': 'bp_rp',
                     'r': 'rmag'},
            }

G_r_, bp_rp_ = [], []
for cat, metadata in catalogs.items():
    table = Table.read(metadata["file"], format="ascii.ecsv")
    #table_final = vstack([table])
    # color gaia and ps
    G_r = table[metadata["G"]] - table[metadata["r"]]
    bp_rp = table[metadata["BP-RP"]]
    for i, ii in zip(G_r, bp_rp):
        if (np.isfinite(i) and np.isfinite(ii)):
            G_r_.append(i)
            bp_rp_.append(ii)

# PN and stars
G_r_pn = tab_pn['Gmag_1'] - tab_pn['rmag']
bp_rp_pn = tab_pn['BP-RP_1']

G_r_star = tab_star['Gmag'] - tab_star['rmag']
bp_rp_star = tab_star['BP-RP']

G_r_star2 = tab_star2['Gmag_1'] - tab_star2['rmag_x']
bp_rp_star2 = tab_star2['BPmag_1'] -  tab_star2['RPmag_1']
 
#PLOTS
colors = [ "light brown", "light orange", "cerulean", "dark pink", "purple", 
           "forest green", ]
#colors = ["light orange", "purple", "purple", "purple", "purple", "purple", "forest green",]
colors = sns.xkcd_palette(colors)
zorder = [1, 9, 8, 5, 3, 1, 10]
size = [30, 10, 10, 10, 10, 10, 50]
alpha = [0.4, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6]

#Plotting
lgd_kws = {'frameon': True, 'fancybox': True, 'shadow': True}
sns.set_style('ticks')
fig, ax = plt.subplots(figsize=(12, 10))
plt.tick_params(axis='x', labelsize=30) 
plt.tick_params(axis='y', labelsize=30)
#fig = plt.figure(figsize=(10, 8))
#ax = fig.add_subplot(111)
ax.scatter(bp_rp_pn, G_r_pn, s=80, c = sns.xkcd_palette(["forest green"]), edgecolors= "g", zorder = 1, lw=0.5, alpha = 0.7, label = "PN")
#ax.scatter(bp_rp_star, G_r_star, s=50, c = sns.xkcd_palette(["forest green"]), edgecolors= "k", zorder = 11, alpha = 0.6, label = "Stars")

# sns.kdeplot(
#     bp_rp_star, G_r_star,
#     ax=ax,
#     bw_method='scott',
#     levels=[0.05, 0.08, 0.1, 0.2, 0.4, 0.5, 0.7, 0.8, 1],
#     norm=PowerNorm(0.5),
#     cmap="Blues",
#     zorder = 3
# )

sns.kdeplot(
    bp_rp_star2, G_r_star2,
    ax=ax,
    bw_method= 'scott',
    levels=[0.05, 0.08, 0.1, 0.2, 0.4, 0.5, 0.7, 0.8, 1],
    norm=PowerNorm(0.5),
    cmap="Blues",
    zorder = 4, 
)

sns.kdeplot(
    bp_rp_, G_r_,
    ax=ax,
    bw_method='scott',
    levels=[0.05, 0.08, 0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 1],
    norm=PowerNorm(0.5),
    cmap="Reds", zorder = 2,
     )
#ax.scatter(bp_rp_, G_r_, s=50, c = sns.xkcd_palette(["purple"]), edgecolors= "k", zorder = 10, alpha = 0.6, label = "Other EM")
plt.xlabel(r'$G_{BP} - G_{RP}$', fontsize=30)
plt.ylabel(r'$G - r$', fontsize=30)
ax.set_xlim(-1.2, 5.7)
ax.set_ylim(-3.3, 9.2)

# Region where are located the PNe
result = findIntersection(1/1.5, 0.15, -7, 23.5, 0.0)

x_new = np.linspace(-15.5, result,  200)
y = (1/1.5)*x_new + 0.15
yy = -7*x_new + 23.5
#Mask
#mask = y >= result_y - 0.5
ax.plot(x_new, y, color='k', linestyle='-.')
ax.plot(x_new, yy , color='k', linestyle='-.')

#reshape
x_new = x_new.ravel()
y = y.ravel()
yy = yy.ravel()
plt.fill_between(x_new, y, yy, color="k", alpha=0.1)

# Some texts
ax.text(0.05, 0.9, "PN zone", fontsize=30,
        bbox=dict(facecolor='gray', alpha=0.0),
        transform=ax.transAxes)

bbox_props = dict(boxstyle="round", fc="w", ec="0.78", alpha=0.5, pad=0.1)
plt.text(0.3, 0.3, 'Stars',
         transform=ax.transAxes, c="b", weight='bold', fontsize=18.8, bbox=bbox_props)

plt.text(0.1, 0.12, 'Others EMO',
         transform=ax.transAxes, c="orange", weight='bold', fontsize=14.8, bbox=bbox_props)

ax.legend(prop={'family': 'monospace', 'size': 30}, **lgd_kws)
#plt.gca().invert_yaxis()
fig.savefig("Figs/color-diagram-ps-gaiaEDR3.pdf")
plt.clf()
