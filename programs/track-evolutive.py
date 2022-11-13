import __future__
import numpy as np
import json
import matplotlib.pyplot as plt
from  astropy.table import Table, vstack, hstack
import pandas as pd
from matplotlib.colors import PowerNorm
import seaborn as sns
from scipy.optimize import fsolve
#from scipy.interpolate import interp1d
from pathlib import Path
ROOT_PATH = Path("data")

table3 = Table.read("table3.dat", format="ascii")
table4 = Table.read("table4.dat", format="ascii")
table5 = Table.read("table5.dat", format="ascii")


m1 = (table3["M"] == 1) & (table3["Z"] == 0.016)
T1 = table3["log(Teff)"][m1]
L1 = table3["log(L)"][m1]
M1 =  table3["M"][m1]

m15 = (table3["M"] == 1.5) & (table3["Z"] == 0.016)
T15 = table3["log(Teff)"][m15]
L15 = table3["log(L)"][m15]
M15 = table3["M"][m15]

m2 = (table3["M"] == 2) & (table3["Z"] == 0.016)
T2 = table3["log(Teff)"][m2]
L2 = table3["log(L)"][m2]
M2 = table3["M"][m2]

m25 = (table3["M"] == 2.5) & (table3["Z"] == 0.016)
T25 = table3["log(Teff)"][m25]
L25 = table3["log(L)"][m25]
M25 = table3["M"][m25]

m35 = (table3["M"] == 3.5) & (table3["Z"] == 0.016)
T35 = table3["log(Teff)"][m35]
L35 = table3["log(L)"][m35]
M35 = table3["M"][m35]

m5 = (table3["M"] == 5) & (table3["Z"] == 0.016)
T5 = table3["log(Teff)"][m5]
L5 = table3["log(L)"][m5]
M5 = table3["M"][m5]

#list for label
Tmin = []
Lmax = []
Munique = [1.0, 1.5, 2.0, 2.5, 3.5, 5.0]

Tmin.append(T1.min())
Tmin.append(T15.min())
Tmin.append(T2.min())
Tmin.append(T25.min())
Tmin.append(T35.min())
Tmin.append(T5.min())

Lmax.append(L1.max())
Lmax.append(L15.max())
Lmax.append(L2.max())
Lmax.append(L25.max())
Lmax.append(L35.max())
Lmax.append(L5.max())


#Plotting
lgd_kws = {'frameon': True, 'fancybox': True, 'shadow': True}
sns.set_style('ticks')
fig, ax = plt.subplots(figsize=(10, 11))
plt.tick_params(axis='x', labelsize=30) 
plt.tick_params(axis='y', labelsize=30)
#fig = plt.figure(figsize=(10, 8))
#ax = fig.add_subplot(111)
ax.plot(T1, L1, 'k-', zorder = 1, lw=1, alpha = 0.9)
ax.plot(T15, L15, 'k-', zorder = 1, lw=1, alpha = 0.9)
ax.plot(T2, L2, 'k-', zorder = 1, lw=1, alpha = 0.9)
ax.plot(T25, L25, 'k-', zorder = 1, lw=1, alpha = 0.9)
ax.plot(T35, L35, 'k-', zorder = 1, lw=1, alpha = 0.9)
ax.plot(T5, L5, 'k-', zorder = 1, lw=1, alpha = 0.9)
#ax.scatter(bp_rp_star, G_r_star, s=50, c = sns.xkcd_palette(["forest green"]), edgecolors= "k", zorder = 11, alpha = 0.6, label = "Stars")

#ax.scatter(bp_rp_, G_r_, s=50, c = sns.xkcd_palette(["purple"]), edgecolors= "k", zorder = 10, alpha = 0.6, label = "Other EM")
plt.xlabel(r'$\mathrm{\log(T_{eff})(K)}$', fontsize=30)
plt.ylabel(r'$\mathrm{\log(L/L_{\odot})}$', fontsize=30)
ax.set_xlim(3.8, 5.7)
ax.set_ylim(0.8, 4.7)

# Putting the masses
bbox_props = dict(boxstyle="round", fc="w", ec="0.88", alpha=0.1, pad=0.1)
for label_, x, y in zip(Munique, Tmin, Lmax):
    ax.annotate(label_, (x, y), alpha=1, size=22,
                   xytext=(40.0, -13.0), textcoords='offset points', ha='right', va='bottom', rotation=0, bbox=bbox_props, zorder=200)
#Mask
#mask = y >= result_y - 0.5

# # Some texts
# ax.text(0.05, 0.9, "PN zone", fontsize=30,
#         bbox=dict(facecolor='gray', alpha=0.0),
#         transform=ax.transAxes)

# bbox_props = dict(boxstyle="round", fc="w", ec="0.78", alpha=0.5, pad=0.1)
# plt.text(0.3, 0.3, 'Stars',
#          transform=ax.transAxes, c="b", weight='bold', fontsize=12.8, bbox=bbox_props)

# plt.text(0.1, 0.2, 'Others EMO',
#          transform=ax.transAxes, c="orange", weight='bold', fontsize=10.8, bbox=bbox_props)

#ax.legend(prop={'family': 'monospace', 'size': 30}, **lgd_kws)
plt.gca().invert_xaxis()
plt.tight_layout()
fig.savefig("tarck-evolution.pdf")
plt.clf()
