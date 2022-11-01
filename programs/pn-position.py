import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from  astropy.table import Table, vstack, hstack
import seaborn as sns
from astropy.coordinates import SkyCoord
from astropy import units as u

df = pd.read_csv('Luis_afxcge1f2c.csv')
dat = Table.read("Halo-PNe.dat", format="ascii")
dat_cand = Table.read("pn-candidates-gaiaDR3.ecsv", format="ascii.ecsv")

l = df["Glon"]
b = df[ "Glat"]

# Converting pandas astropy
dat_ = Table.from_pandas(df)
id1 = dat_["PNG"]

#Halo
id2 = dat["ID"]

mask = np.array([source in id2 for source in id1])

# Converting pandas astropy and masking
dat_ = Table.from_pandas(df)
dat_h = dat_[mask]
l_ = dat_h["Glon"]
b_ = dat_h[ "Glat"]

print(len(dat_h))

# PN candidates
id_cand = dat_cand["LAMOST_1"]


PN = ["LAMOST J020808.63+491401.0"]
mask1 = np.array([source in PN for source in id_cand])
dat_cand_ = dat_cand[mask1]
ra = dat_cand_["RAJ2000_1"]
dec = dat_cand_["DEJ2000_1"]

icrs = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
gal = icrs.galactic  

l_rad = gal.l.radian
l_rad[l_rad > np.pi] -= 2. * np.pi
b_rad = gal.b.radian

#latitude in degrees
b_deg = b_rad * (180/np.pi)#gal.b.degree#b_rad * (180/np.pi)
l_deg = l_rad * (180/np.pi)

#Plotting
lgd_kws = {'frameon': True, 'fancybox': True, 'shadow': True}
sns.set_style('ticks')
fig, ax = plt.subplots(figsize=(12, 6))
plt.tick_params(axis='x', labelsize=30) 
plt.tick_params(axis='y', labelsize=30)
plt.xlabel(r'$l(Gal)$', fontsize=30)
plt.ylabel(r'$b(Gal)$', fontsize=30)
# ax.set_xlim(-1.2, 5.7)
# ax.set_ylim(-3.3, 9.2)
#fig = plt.figure(figsize=(10, 8))
#ax = fig.add_subplot(111)
ax.scatter(l, b, marker = ".", c = sns.xkcd_palette(["forest green"]), edgecolors= "g", zorder = 1, lw=1, alpha = 0.7)
ax.scatter(l_, b_, s = 100, c = sns.xkcd_palette(["purple"]), edgecolors= "purple", zorder = 2, lw=2)
ax.scatter(l_deg, b_deg, s = 100, c = sns.xkcd_palette(["azure"]), edgecolors= "b", zorder = 3, lw=2)

bbox_props = dict(boxstyle="round", fc="w", ec="0.78", alpha=0.6, pad=0.1)
for label_, x, y in zip(dat_h["Name"], l_, b_):
    ax.annotate(label_, (x, y), alpha=1, size=8,
                   xytext=(25.0, 8.6), textcoords='offset points', ha='right', va='bottom', bbox=bbox_props, zorder=100)

#ax.legend(prop={'family': 'monospace', 'size': 30}, **lgd_kws)
#plt.gca().invert_yaxis()
plt.tight_layout()
fig.savefig("Figs/galctic-coord-pn.pdf")
plt.clf()

