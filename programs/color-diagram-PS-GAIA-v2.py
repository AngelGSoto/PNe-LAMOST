import numpy as np
import json
import matplotlib.pyplot as plt
from  astropy.table import Table, vstack, hstack
import pandas as pd
from matplotlib.colors import PowerNorm
import seaborn as sns
from pathlib import Path
ROOT_PATH = Path("data")

# Planetary nebula
tab_pn = Table.read("PN-PS-Gaia.ecsv", format="ascii.ecsv")
tab_star = Table.read("Star-PS-Gaia.ecsv", format="ascii.ecsv")

catalogs = {'SySt': {'file': 'SySt-PS-Gaia.ecsv',
                     'G': 'Gmag',
                     'BP-RP': 'BP-RP',
                     'r': 'rmag'},
            'CV': {'file': 'CV-PS-Gaia.ecsv',
                    'G': 'Gmag',
                   'BP-RP': 'BP-RP',
                   'r': 'rmag'},
            'SNR': {'file': 'SNR-PS-Gaia.ecsv',
                    'G': 'Gmag',
                    'BP-RP': 'BP-RP',
                    'r': 'rmag'},
            'YSO': {'file': 'YSO-PS-Gaia.ecsv',
                    'G': 'Gmag',
                    'BP-RP': 'BP-RP',
                    'r': 'rmag'},
            'AeBe': {'file': 'AeBe-PS-Gaia.ecsv',
                     'G': 'Gmag',
                     'BP-RP': 'BP-RP',
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
            print(i, ii)
            G_r_.append(i)
            bp_rp_.append(ii)

# PN and stars
G_r_pn = tab_pn['Gmag_1'] - tab_pn['rmag']
bp_rp_pn = tab_pn['BP-RP_1']

G_r_star = tab_star['Gmag'] - tab_star['rmag']
bp_rp_star = tab_star['BP-RP']
 
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
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)
ax.scatter(bp_rp_pn, G_r_pn, s=30, c = sns.xkcd_palette(["forest green"]), edgecolors= "w", zorder = 1, lw=0.5, alpha = 0.6, label = "PN")
#ax.scatter(bp_rp_star, G_r_star, s=50, c = sns.xkcd_palette(["forest green"]), edgecolors= "k", zorder = 11, alpha = 0.6, label = "Stars")

sns.kdeplot(
    bp_rp_star, G_r_star,
    ax=ax,
    bw_method='scott',
    levels=[0.1, 0.2, 0.4, 0.5, 0.7, 0.8, 1],
    norm=PowerNorm(0.5),
    cmap="Blues",
    zorder = 3
)

sns.kdeplot(
    bp_rp_, G_r_,
    ax=ax,
    bw_method='scott',
    levels=[0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 1],
    norm=PowerNorm(0.5),
    cmap="Reds", zorder = 2,
     )
#ax.scatter(bp_rp_, G_r_, s=50, c = sns.xkcd_palette(["purple"]), edgecolors= "k", zorder = 10, alpha = 0.6, label = "Other EM")

plt.xlabel(r'$G_{BP} - G_{RP}$')
plt.ylabel(r'$G - r$')
#ax.set_xlim(-30.0, 390.0)
#ax.set_ylim(-90.0, 90.0)
ax.legend(prop={'family': 'monospace', 'size': 'x-large'}, **lgd_kws)
#plt.gca().invert_yaxis()
fig.savefig("Figs/color-diagram-ps-gaia-v2.pdf")
plt.clf()
    
