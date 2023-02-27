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
true = dat_["PNstat"] == "T"
proba = dat_["PNstat"] == "P"
like = dat_["PNstat"] == "L"
dat_true = dat_[true]
dat_proba = vstack([dat_[proba], dat_[like]])

l_true = dat_true["Glon"]
b_true = dat_true[ "Glat"]

l_proba = dat_proba["Glon"]
b_proba = dat_proba[ "Glat"]
print(len(dat_true))
print(len(dat_proba))

ra_true = dat_true["DRAJ2000"]
dec_true = dat_true[ "DDECJ2000"]
icrs_true = SkyCoord(ra=ra_true*u.degree, dec=dec_true*u.degree, frame='icrs')
gal_true = icrs_true.galactic
l_rad_true = gal_true.l.radian
l_rad_true[l_rad_true > np.pi] -= 2. * np.pi
b_rad_true = gal_true.b.radian
b_deg_true = b_rad_true * (180/np.pi)#gal.b.degree#b_rad * (180/np.pi)
l_deg_true = l_rad_true * (180/np.pi)

ra_proba = dat_proba["DRAJ2000"]
dec_proba = dat_proba[ "DDECJ2000"]
icrs_proba = SkyCoord(ra=ra_proba*u.degree, dec=dec_proba*u.degree, frame='icrs')
gal_proba = icrs_proba.galactic
l_rad_proba = gal_proba.l.radian
l_rad_proba[l_rad_proba > np.pi] -= 2. * np.pi
b_rad_proba = gal_proba.b.radian
b_deg_proba = b_rad_proba * (180/np.pi)#gal.b.degree#b_rad * (180/np.pi)
l_deg_proba = l_rad_proba * (180/np.pi)

#Halo
id2 = dat["ID"]
mask = np.array([source in id2 for source in id1])

# Converting pandas astropy and masking
dat_h = dat_[mask]
l_ = dat_h["Glon"]
b_ = dat_h[ "Glat"]

#print(len(dat_h))

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
# disntace fform GAIA
D = 2.31300952
z = D * np.sin(b_rad)

print("IzI =", z)

#latitude in degrees
b_deg = b_rad * (180/np.pi)#gal.b.degree#b_rad * (180/np.pi)
l_deg = l_rad * (180/np.pi)
print(b_deg, l_deg)
#Plotting
lgd_kws = {'frameon': True, 'fancybox': True, 'shadow': True}
sns.set_style('ticks')
fig, ax = plt.subplots(figsize=(12, 6))
# fig = plt.figure(figsize=(14,7))
# ax = fig.add_subplot(1,1,1, projection='aitoff')
plt.tick_params(axis='x', labelsize=30) 
plt.tick_params(axis='y', labelsize=30)
plt.xlabel(r'$l$(Galactic longitude)[deg]', fontsize=28)
plt.ylabel(r'$b$(Galactic latitude)[deg]', fontsize=28)
# ax.set_xlim(-1.2, 5.7)
# ax.set_ylim(-3.3, 9.2)
#fig = plt.figure(figsize=(10, 8))
#ax = fig.add_subplot(111)
ax.scatter(l_deg_true, b_deg_true, marker = "o", s = 30, c = sns.xkcd_palette(["forest green"]), edgecolors= "g", zorder = 1, lw=1, alpha = 0.6, label="Confirmed HASH PNe")
ax.scatter(l_deg_proba, b_deg_proba, marker = "o", s = 30, c = sns.xkcd_palette(["white"]), edgecolors= "k", zorder = 2, lw=0.7, alpha = 0.5, label="No confirmed HASH PNe")
#ax.scatter(l_, b_, s = 100, c = sns.xkcd_palette(["purple"]), edgecolors= "purple", zorder = 2, lw=2)
ax.scatter(l_deg, b_deg, s = 150, c = sns.xkcd_palette(["azure"]), edgecolors= "b", zorder = 3, lw=2, label="New PN")

bbox_props = dict(boxstyle="round", fc="w", ec="0.78", alpha=0.3, pad=0.1)
for x, y in zip(l_deg, b_deg):
    ax.annotate("J020808.63+491401.0", (x, y), alpha=1, size=10, c="blue",
                   xytext=(45.0, -20.6), textcoords='offset points', ha='right', va='bottom', bbox=bbox_props, zorder=100)


ax.grid(True, linestyle='-.', linewidth=0.7)
ax.legend(prop={'family': 'monospace', 'size': 15}, **lgd_kws)
#plt.gca().invert_yaxis()
plt.tight_layout()
fig.savefig("Figs/galctic-coord-pn.pdf")
plt.clf()

