'''
Script to dealing with output data from Cloudy
Based in pyCloudy (Morisset, C., 2013, pyCloudy, Astrophysics Source Code Library)
Author: Luis A. Gutiérrez Soto
30/01/2023
'''
import numpy as np
from astropy.table import Table
from astropy.io import ascii
import matplotlib.pyplot as plt
import os
import pyCloudy as pc
from astropy.io import fits
import numpy as np
import seaborn as sn
import argparse
import sys
#sn.set_context("poster")

parser = argparse.ArgumentParser(
    description="""Reading the ouput cloudy models""")

parser.add_argument("source", type=str,
                    default="model_100000_36.58",
                    help="Name of input model ")

cmd_args = parser.parse_args()
file_ = cmd_args.source 

#model_name_ = file_ + ".in"
model_name = file_

# Open file
f = open(model_name + ".in", 'r')
header1 = f.readline()
header2 = f.readline()
header3 = f.readline()
header4 = f.readline()
header5 = f.readline()
Te, Lu, Denss = [], [], []
for line in f:
    line = line.strip()
    columns = line.split()
    T = (columns[0] == "blackbody") 
    L = columns[0] == "luminosity"
    Dens = columns[0] == "hden"
    Te.append(columns[T])
    Lu.append(columns[L])
    Denss.append(columns[Dens])
    
Tf = "{:.0f}".format(float(Te[0]))
Lf = "{:.0f}".format((10**float(Lu[1])) / 3.839e33)
densf = "{:.0f}".format(10**float(Denss[8]))

# Chi-square and chi-reduce
tab = Table.read("better-" + model_name + "_0.0.ecsv", format="ascii.ecsv")
chi = "{:.2f}".format(float(tab["chi"]))
chi_r = "{:.2f}".format(float(tab["Chi red"]))

# Reading the Cloudy outputs in the Mod CloudyModel object
Mod = pc.CloudyModel(model_name)

print(dir(Mod)) # This is the online answering way
print(Mod.print_stats())
print(Mod.get_ab_ion_vol_ne('O',2))
#print("Halpha/Hbeta", Mod.get_emis_vol('H__1_6562.81A')/Mod.get_emis_vol('H__1_4861.33A'))
#print(Mod.print_lines())
print(Mod.Teff)
print(Mod.nH_mean)
print(Mod.te.mean())

#Getting wl and flux
wl = Mod.get_cont_x(unit='Ang')
flux = Mod.get_cont_y(cont = 'total', unit = 'esAc')
#Ordered lambda and flux 
wll, flux = zip(*sorted(zip(wl, flux)))

data = Table([wll, flux], names=('Wl', 'Flux'), meta={'name': 'first table'})
mask = (data["Wl"] > 3000) & (data["Wl"] < 9000)
data_mask = data[mask]

# OUR PN
hdu = fits.open("../Spectra-lamostdr7/spec-56581-VB031N50V1_sp08-218.fits")
hdudata = hdu[0].data
wl = hdudata[2]
Flux = hdudata[0]

def closest(lst, K):
    '''find the closest number'''
    lst = np.array(lst)
    idx = (np.abs(lst - K)).argmin()
    return lst[idx]

# Model
MaskHbeta = (data_mask["Wl"] >= (4863 - 50)) & (data_mask["Wl"] <= (4863 + 50))
HBeta = data_mask[MaskHbeta]
max_HBeta = HBeta["Flux"].max()
flux_m = data_mask["Flux"] / max_HBeta

# Our PN
MaskHbeta_our = (wl >= (4863 - 50)) & (wl <= (4863 + 50))
HBeta_our = Flux[MaskHbeta_our]
max_HBeta_our = HBeta_our.max()
flux_our = Flux / max_HBeta_our

#spliting the spectra in the blue and red arms
m_blue = (wl >= min(wl)) & (wl <= 5800)
m_red = (wl >= 6300) & (wl <= max(wl))
wl_blue = wl[m_blue]
wl_red = wl[m_red]
Flux_blue = flux_our[m_blue]
Flux_red = flux_our[m_red]

fig, ax = plt.subplots(figsize=(11, 5))
#ax.set_title(namefile)
ax.set(xlim=[3600, 5900])
#plt.ylim(ymin=-200,ymax=1500)
plt.tick_params(axis='x', labelsize=16) 
plt.tick_params(axis='y', labelsize=16)
#plt.ylim(ymin=0.0,ymax=500)
# ax.set(xlabel='Wavelength $(\AA)$')
# ax.set(ylabel='Relative flux')
ax.set_xlabel('Wavelength $(\AA)$', fontsize=16)
ax.set_ylabel('Normalised flux', fontsize=16)
# ax.set(xlabel='Wavelength $(\AA)$')
# ax.set(ylabel='Normalised flux')
plt.plot(data_mask["Wl"], flux_m,  c = "darkolivegreen", linewidth=0.8, zorder=3, label = 'Model ')
#plt.plot(wl, flux_our, c = "white", linewidth=900, alpha =0.6, zorder=1)
ax.plot(wl_blue, Flux_blue, c = "blueviolet", linewidth=2.5, alpha =0.7, zorder=5, label = 'J020808.63+491401.0')
#ax.plot(wl_red, Flux_red, c = "blueviolet", linewidth=0.7, zorder=5)
bbox_props = dict(boxstyle="round", fc="w", ec="0.88", alpha=0.6, pad=0.1)

info_mod = [r"$\mathrm{T_{eff}=}$" + Tf + "K", r"L=" + Lf + r"$\mathrm{L_{\odot}}$", r"$n_{H}=$" + densf + r"$\mathrm{cm^{-3}}$", r"$\chi^2=$" + chi]
loc_text = [(0.70, 0.7), (0.70, 0.7-0.13), (0.70, 0.7-2*0.13), (0.70, 0.7-3*0.13)]
for loc_, taxt in zip(loc_text, info_mod):
    ax.text(loc_[0], loc_[1], taxt, fontsize=19,
            bbox=dict(facecolor='gray', alpha=0.0),
            transform=ax.transAxes)
ax.legend(loc="upper left", fontsize="x-large")
#sn.despine()
plt.tight_layout()
plt.savefig(model_name + "-blue.pdf")

fig, ax1 = plt.subplots(figsize=(11, 5))
#ax.set_title(namefile)
ax1.set(xlim=[5900, 9100])
#plt.ylim(ymin=-200,ymax=1500)
plt.tick_params(axis='x', labelsize=16) 
plt.tick_params(axis='y', labelsize=16)
#plt.ylim(ymin=0.0,ymax=500)
# ax.set(xlabel='Wavelength $(\AA)$')
# ax.set(ylabel='Relative flux')
ax1.set_xlabel('Wavelength $(\AA)$', fontsize=16)
ax1.set_ylabel('Normalised flux', fontsize=16)
plt.plot(data_mask["Wl"], flux_m,  c = "darkolivegreen", linewidth=0.8, zorder=3, label = 'Model ')
#plt.plot(wl, flux_our, c = "white", linewidth=0.7, zorder=1)
#ax1.plot(wl_blue, Flux_blue, c = "blueviolet", linewidth=0.7, zorder=5, label = 'J020808.63+491401.0')
ax1.plot(wl_red, Flux_red, c = "blueviolet", linewidth=2.5, alpha =0.7, zorder=5)
# bbox_props = dict(boxstyle="round", fc="w", ec="0.88", alpha=0.6, pad=0.1)
# info_mod = [r"$\mathrm{T_{eff}=}$" + Tf + "K", r"L=" + Lf + r"$\mathrm{L_{\odot}}$", r"$n_{H}=$" + densf + r"$\mathrm{cm^{-3}}$", r"$\chi^2=$" + chi]
# loc_text = [(0.70, 0.7), (0.70, 0.7-0.13), (0.70, 0.7-2*0.13), (0.70, 0.7-3*0.13)]
# for loc_, taxt in zip(loc_text, info_mod):
#     ax1.text(loc_[0], loc_[1], taxt, fontsize=19,
#             bbox=dict(facecolor='gray', alpha=0.0),
#             transform=ax1.transAxes)

#ax1.legend(loc="upper left", fontsize="x-small")
#sn.despine()
plt.tight_layout()
plt.savefig(model_name + "-red.pdf")
