import numpy as np
import json
import matplotlib.pyplot as plt
from  astropy.table import Table
import pandas as pd
import seaborn as sns
from pathlib import Path
ROOT_PATH = Path("data")

catalogs = {'PN': {'file': 'PN-PS-Gaia.ecsv',
                   'G': 'Gmag_1',
                   'BP-RP': 'BP-RP_1',
                   'r': 'rmag'},
            'SySt': {'file': 'SySt-PS-Gaia.ecsv',
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

colors = [ "light brown", "light orange", "cerulean", "dark pink", "purple", 
           "forest green", ]
colors = sns.xkcd_palette(colors)
zorder = [10, 9, 8, 5, 3, 1]
#Plotting
lgd_kws = {'frameon': True, 'fancybox': True, 'shadow': True}
sns.set_style('ticks')
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)

for (cat, metadata), color, zorder_ in zip(catalogs.items(), colors, zorder):
    table = Table.read(metadata["file"], format="ascii.ecsv")
    # color gaia and ps
    G_r = table[metadata["G"]] - table[metadata["r"]]
    bp_rp = table[metadata["BP-RP"]]
    
    ax.scatter(bp_rp, G_r, s=100, c = color, edgecolors= "k", zorder = zorder_, label = cat)
    
    #ax.scatter(bp_rpmag_wd, g_i_wd, color="gray", cmap="seismic", alpha = 0.5, s=5, label = "WD")

plt.xlabel(r'$G_{BP} - G_{RP}$')
plt.ylabel(r'$G - r$')
#ax.set_xlim(-30.0, 390.0)
#ax.set_ylim(-90.0, 90.0)
ax.legend(prop={'family': 'monospace', 'size': 'x-large'}, **lgd_kws)
#plt.gca().invert_yaxis()
fig.savefig("Figs/color-diagram-ps-gaia.pdf")
plt.clf()
    




    
