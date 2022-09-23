import numpy as np
import json
import matplotlib.pyplot as plt
from  astropy.table import Table
import pandas as pd
import seaborn as sns
from pathlib import Path
ROOT_PATH = Path("data")

#PNs
df_pn_all = pd.read_csv(ROOT_PATH / "Luis_hash-pn-gaia-panstarra.csv")
mask = df_pn_all["PNstat_1"] == "T"
df_pn = df_pn_all[mask]
#CVs
df_cv = pd.read_csv(ROOT_PATH / "cvs-gaiadr3-panstarrs.csv")
#SN
df_sn = pd.read_csv(ROOT_PATH / "snr-gaia-ps.csv")
#WD
# merged_table_list = []
# for chunk in pd.read_csv(ROOT_PATH / "WD-gaia-ps.csv", chunksize=50000):
#     merged_table_list.append(chunk)

# df_wd_all_c = pd.concat(merged_table_list)
# col = ["Pwd", "Gmag", "BP-RP", "gmag", "rmag", "imag"]
# df_wd_all = df_wd_all_c[col]
# mask_wd = df_wd_all["Pwd"] >= 0.95
# df_wd = df_wd_all[mask_wd]

# PNe
Gmag_pn = df_pn["phot_g_mean_mag"]
bp_rpmag_pn = df_pn["bp_rp"]
# PS
g_ps_pn = df_pn["gmag"]
r_ps_pn = df_pn["rmag"]
i_ps_pn = df_pn["imag"]
G_r_ps_pn = Gmag_pn - r_ps_pn
g_i_pn = g_ps_pn - i_ps_pn
g_r_pn = g_ps_pn - r_ps_pn
r_i_pn = r_ps_pn - i_ps_pn

colors = ["cerulean",]
colors = sns.xkcd_palette(colors)


#CV
Gmag_cv = df_cv["phot_g_mean_mag"]
bp_rpmag_cv = df_cv["bp_rp"]
# PS
g_ps_cv = df_cv["gmag"]
r_ps_cv = df_cv["rmag"]
i_ps_cv = df_cv["imag"]
G_r_ps_cv = Gmag_cv - r_ps_cv
g_i_cv = g_ps_cv - i_ps_cv
g_r_cv = g_ps_cv - r_ps_cv
r_i_cv = r_ps_cv - i_ps_cv

colors1 = ["pale yellow"]
colors1 = sns.xkcd_palette(colors1)

#SN
Gmag_sn = df_sn["phot_g_mean_mag"]
bp_rpmag_sn = df_sn["bp_rp"]
# PS
g_ps_sn = df_sn["gmag"]
r_ps_sn = df_sn["rmag"]
i_ps_sn = df_sn["imag"]
G_r_ps_sn = Gmag_sn - r_ps_sn
g_i_sn = g_ps_sn - i_ps_sn
g_r_sn = g_ps_sn - r_ps_sn
r_i_sn = r_ps_sn - i_ps_sn

colors2 = ["purple"]
colors2 = sns.xkcd_palette(colors2)

#WD
# Gmag_wd = df_wd["Gmag"]
# bp_rpmag_wd = df_wd["BP-RP"]
# # PS
# g_ps_wd = df_wd["gmag"]
# r_ps_wd = df_wd["rmag"]
# i_ps_wd = df_wd["imag"]
# G_r_ps_wd = Gmag_wd - r_ps_wd
# g_i_wd = g_ps_wd - i_ps_wd
# g_r_wd = g_ps_wd - r_ps_wd
# r_i_wd = r_ps_wd - i_ps_wd

colors3 = ["black"]
colors3 = sns.xkcd_palette(colors3)

#PLOTS
lgd_kws = {'frameon': True, 'fancybox': True, 'shadow': True}
sns.set_style('ticks')
fig = plt.figure(figsize=(7, 7))
ax = fig.add_subplot(111)
  
ax.scatter(bp_rpmag_pn, g_i_pn, s=50, c = colors,edgecolor=['none'], label = "PN")
ax.scatter(bp_rpmag_cv, g_i_cv, c = colors1, edgecolor=['black'], alpha = 0.4, s = 50, label = "CV")
ax.scatter(bp_rpmag_sn, g_i_sn, c = colors2, edgecolor=['black'], alpha = 0.4, s = 50, label = "CV")
#ax.scatter(bp_rpmag_wd, g_i_wd, color="gray", cmap="seismic", alpha = 0.5, s=5, label = "WD")

plt.xlabel(r'$G_{BP} - G_{RP}$')
plt.ylabel(r'$g - i$')
#ax.set_xlim(-30.0, 390.0)
#ax.set_ylim(-90.0, 90.0)
ax.legend(prop={'family': 'monospace', 'size': 'x-small'}, **lgd_kws)
#plt.gca().invert_yaxis()
fig.savefig("Figs/color-diagram-gai-PS-gi.pdf")
plt.clf()

ax1 = fig.add_subplot(111)
ax1.scatter(bp_rpmag_pn, r_i_pn, s=50, c = colors,edgecolor=['none'], label = "PN")
ax1.scatter(bp_rpmag_cv, r_i_cv, c = colors1, edgecolor=['black'], alpha = 0.4, s = 50, label = "CV")
ax1.scatter(bp_rpmag_sn, r_i_sn, c = colors2, edgecolor=['black'], alpha = 0.4, s = 50, label = "CV")
#ax1.scatter(bp_rpmag_wd, r_i_wd, color="gray", cmap="seismic", alpha = 0.5, s=5, label = "WD")
plt.xlabel(r'$G_{BP} - G_{RP}$')
plt.ylabel(r'$r - i$')
#ax.set_xlim(-30.0, 390.0)
#ax.set_ylim(-90.0, 90.0)
ax1.legend(prop={'family': 'monospace', 'size': 'x-small'}, **lgd_kws)
#plt.gca().invert_yaxis()
fig.savefig("Figs/color-diagram-gai-PS-ri.pdf")
plt.clf()

ax2 = fig.add_subplot(111)
ax2.scatter(bp_rpmag_pn, g_r_pn, s=50, c = colors,edgecolor=['none'], label = "PN")
ax2.scatter(bp_rpmag_cv, g_r_cv, c = colors1, edgecolor=['black'], alpha = 0.4, s = 50, label = "CV")
ax2.scatter(bp_rpmag_sn, g_r_sn, c = colors2, edgecolor=['black'], alpha = 0.4, s = 50, label = "CV")
#ax2.scatter(bp_rpmag_wd, g_r_wd, color="gray", cmap="seismic", alpha = 0.5, s=5, label = "WD")

plt.xlabel(r'$G_{BP} - G_{RP}$')
plt.ylabel(r'$g - r$')
#ax.set_xlim(-30.0, 390.0)
#ax.set_ylim(-90.0, 90.0)
ax2.legend(prop={'family': 'monospace', 'size': 'x-small'}, **lgd_kws)
#plt.gca().invert_yaxis()
fig.savefig("Figs/color-diagram-gai-PS-gr.pdf")
plt.clf()


ax3 = fig.add_subplot(111)
ax3.scatter(bp_rpmag_pn, G_r_ps_pn, s=50, c = colors,edgecolor=['none'], label = "PN")
ax3.scatter(bp_rpmag_cv, G_r_ps_cv, c = colors1, edgecolor=['black'], alpha = 0.4, s = 50, label = "CV")
ax3.scatter(bp_rpmag_sn, G_r_ps_sn, c = colors2, edgecolor=['black'], alpha = 0.4, s = 50, label = "CV")
#ax3.scatter(bp_rpmag_wd, G_r_ps_wd, color="gray", cmap="seismic", alpha = 0.5, s=5, label = "WD")

plt.xlabel(r'$G_{BP} - G_{RP}$')
plt.ylabel(r'$G - r$')
#ax.set_xlim(-30.0, 390.0)
#ax.set_ylim(-90.0, 90.0)
ax3.legend(prop={'family': 'monospace', 'size': 'x-small'}, **lgd_kws)
#plt.gca().invert_yaxis()
fig.savefig("Figs/color-diagram-gai-PS-Gr.pdf")
plt.clf()
