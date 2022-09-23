import numpy as np
import json
import matplotlib.pyplot as plt
from  astropy.table import Table
import pandas as pd
import seaborn as sns
from pathlib import Path
ROOT_PATH = Path("data")

#PNs
df_pn_all = pd.read_csv(ROOT_PATH / "Luis_hash-pn-gaia.csv")
mask = df_pn_all["PNstat"] == "T"
df_pn = df_pn_all[mask]
#CVs
df_cv = pd.read_csv(ROOT_PATH / "cvs-gaiadr3.csv")
#SNs
df_sn = pd.read_csv(ROOT_PATH / "TAP_4_VII_284_snrs-gaia.csv")

# PNe
Gmag_pn = df_pn["phot_g_mean_mag"]
bp_rpmag_pn = df_pn["bp_rp"]
colors = ["cerulean",]
colors = sns.xkcd_palette(colors)
  
#CV
Gmag_cv = df_cv["phot_g_mean_mag"]
bp_rpmag_cv = df_cv["bp_rp"]
colors1 = ["pale yellow"]
colors1 = sns.xkcd_palette(colors1)

#SN
Gmag_sn = df_sn["phot_g_mean_mag"]
bp_rpmag_sn = df_sn["bp_rp"]
colors2 = ["purple"]
colors2 = sns.xkcd_palette(colors2)

lgd_kws = {'frameon': True, 'fancybox': True, 'shadow': True}
sns.set_style('ticks')
fig = plt.figure(figsize=(6, 7))
ax = fig.add_subplot(111)
  
ax.scatter(bp_rpmag_pn, Gmag_pn, s=50, c = colors,edgecolor=['none'], label = "PN")
ax.scatter(bp_rpmag_cv, Gmag_cv, c = colors1, edgecolor=['black'], alpha = 0.4, s = 50, label = "CV")
plt.scatter(bp_rpmag_sn, Gmag_sn, c = colors2, alpha=0.8, label = "SN")
plt.xlabel(r'$G_{BP} - G_{RP}$')
plt.ylabel(r'$G$')
#ax.set_xlim(-30.0, 390.0)
#ax.set_ylim(-90.0, 90.0)
ax.legend(prop={'family': 'monospace', 'size': 'x-small'}, **lgd_kws)
plt.gca().invert_yaxis()
fig.savefig("Figs/color-diagram-gaia.pdf")    
