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

#Definition
#Find the point inteception between two lines     
def findIntersection(m, y, m1, y1, x0):
    x = np.linspace(-10.0, 15.5, 200)
    return fsolve(lambda x : (m*x + y) - (m1*x + y1), x0)

# Readeing the file
#tab1 = Table.read("cans-hue-PS-GaiaEDR3.ecsv", format="ascii.ecsv")
tab = Table.read("cans-new-PS-GaiaEDR3.ecsv", format="ascii.ecsv")
#tab3 = Table.read("cans-sim-PS-GaiaEDR3.ecsv", format="ascii.ecsv")

#tab_ff = vstack([tab1, tab2, tab3])

G_r = tab['phot_g_mean_mag'] - tab['rmag_x']
bp_rp = tab['bp_rp']

#applying the color criteria
c_eq1 = (1/1.5)*bp_rp + 0.15
c_eq2 = -7*bp_rp + 23.5

mask = (G_r >= c_eq1) & (G_r <= c_eq2)
tab_f = tab[mask]
for i, ii in zip(tab_f["RAJ2000_1"], tab_f["DEJ2000_1"]):
    print(i, ii)

G_r_f = tab_f['phot_g_mean_mag'] - tab_f['rmag_x']
bp_rp_f = tab_f['bp_rp']

tab_f["G_r"] = G_r_f

#Columns
col = ["recno_1", "SpecL_1", "LAMOST_1", "RAJ2000_1", "DEJ2000_1", "RAJ2000_x", "DEJ2000_x",
       "gmag", "e_gmag", "rmag_x", "e_rmag", "imag", "e_imag", "zmag", "e_zmag", "ymag", "e_ymag",       "ra", "dec", "parallax", "parallax_error", "parallax_over_error", "phot_g_mean_mag",
       "phot_bp_mean_mag", "phot_rp_mean_mag", "G_r", "bp_rp", "angDist_2"]

#Saving file
for i in tab_f:
        print(i["LAMOST_1"], i["G_r"], i["bp_rp"])
        
mask_inf = tab_f["G_r"] != float('inf')
tab_f[mask_inf][col].write("pn-candidates-gaiaDR3.ecsv", format="ascii.ecsv", overwrite=True)

# For lamost
n = len(tab_f["RAJ2000_1"])
print("Numbers:", n)
print("Objects:", tab_f["LAMOST_1"])
print("G-r:", tab_f["G_r"][mask_inf])

sep = np.linspace(2.0, 2.0, num=n)
ra = tab_f["RAJ2000_1"]
dec = tab_f["DEJ2000_1"]
table = Table([ra, dec, sep], names=('ra', 'dec', 'radius'), meta={'name': 'first table'})
asciifile = "cans-new-PS-GaiaEDR3-coorLamost.dat"
table.write(asciifile, format="ascii.commented_header", delimiter=',', overwrite=True)

#Plotting
lgd_kws = {'frameon': True, 'fancybox': True, 'shadow': True}
sns.set_style('ticks')
fig, ax = plt.subplots(figsize=(12, 10))
plt.tick_params(axis='x', labelsize=30) 
plt.tick_params(axis='y', labelsize=30)
plt.xlabel(r'$G_{BP} - G_{RP}$', fontsize=30)
plt.ylabel(r'$G - r$', fontsize=30)
ax.set_xlim(-1.2, 5.7)
ax.set_ylim(-3.3, 9.2)
#fig = plt.figure(figsize=(10, 8))
#ax = fig.add_subplot(111)
ax.scatter(bp_rp_f, G_r_f, s=250, c = sns.xkcd_palette(["denim blue"]), edgecolors= "b", zorder = 1, lw=3, alpha = 0.7, label = "PN candidates")

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

bbox_props = dict(boxstyle="round", fc="w", ec="0.78", alpha=0.6, pad=0.1)
for label_, x, y in zip(tab_f["LAMOST_1"], bp_rp_f, G_r_f):
    ax.annotate(label_, (x, y), alpha=1, size=8,
                   xytext=(25.0, 10.6), textcoords='offset points', ha='right', va='bottom', bbox=bbox_props, zorder=100)

ax.legend(prop={'family': 'monospace', 'size': 30}, **lgd_kws)
#plt.gca().invert_yaxis()
fig.savefig("Figs/pn-candidates-gaiaDR3.pdf")
plt.clf()
