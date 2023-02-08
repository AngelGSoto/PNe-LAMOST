import numpy as np
import numpy.ma as ma
import os
import shutil
import glob
from astropy.table import Table
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 25})
import seaborn as sns

file_list_escv = glob.glob("*.ecsv")

Name, Chi = [], []
for file_name in file_list_escv:
    tab = Table.read(file_name)
    if  (tab["Chi red"] >= 1):
        Namee = file_name.split("ter-")[-1].split(".ecs")[0]
        chi = tab["chi"]
        Name.append(Namee)
        Chi.append(chi)

pattern = "../Model-PNNew/models/N2242_He091O382Ar624_itere_R166R1705_final/*.in"
file_list_in = glob.glob(pattern)

Name_, Te, Lu, Denss = [], [], [], []
for filename in file_list_in:
    Namee = filename.split("final/")[-1].split(".in")[0]
    Name_.append(Namee)
    f = open(filename, 'r')
    header1 = f.readline()
    header2 = f.readline()
    header3 = f.readline()
    header4 = f.readline()
    header5 = f.readline()

    for line in f:
        line = line.strip()
        columns = line.split()
        T = (columns[0] == "blackbody") 
        L = (columns[0] == "luminosity")
        Dens = columns[0] == "hden"
        Te.append(columns[T])
        Lu.append(columns[L])
        Denss.append(columns[Dens])

Teff = np.array(Te)
L_ = np.array(Lu)
Dens_ = np.array(Denss)
mask_value_T = ((Teff != "luminosity") & (Teff != "abundances") & (Teff != "save") & (Teff != "continue") & (Teff != "hden") & (Teff != "radius") & (Teff != "iterate") & (Teff != "sphere"))
mask_value_L = ((L_ != "blackbody") & (L_ != "abundances") & (L_ != "save") & (L_ != "continue") & (L_ != "hden") & (L_ != "radius") & (L_ != "iterate") & (L_ != "sphere"))
mask_value_dens = ((Teff != "luminosity") & (L_ != "blackbody") & (L_ != "abundances") & (L_ != "save") & (L_ != "continue") & (L_ != "radius") & (L_ != "iterate") & (L_ != "sphere"))

Tf = Teff[mask_value_T]
Lf = L_[mask_value_L]
densf = Dens_[mask_value_dens]

Teffpn, Lpn, densn = [], [], []
for ii, jj, kk in zip(Tf, Lf, densf):
    Teffpn.append(np.log10(float(ii)))#/10**5)
    Lpn.append(np.log10(10**float(jj) / 3.839e33))
    densn.append(10**float(kk))

Nombre, X, TT, LL, denss = [], [], [], [], []
for a, c in zip(Name, Chi):
    for b, t, l, d in zip(Name_, Teffpn, Lpn, densn):
        if a==b:
            Nombre.append(a)
            X.append(c)
            TT.append(t)
            LL.append(l)
            denss.append(d)

#Getting the x-square min:
X_TT = np.array(list(zip(X, TT)), dtype=object)

m_Xmin = X_TT[:,0] == X_TT[:,0].min()
X_TT_min = X_TT[m_Xmin]
print(X_TT[m_Xmin])
   
#print(Nombre, X, TT, LL, denss)
print(len(LL), len(denss))
#print("Minimum chi-square and Teff:", X_TT[:, 1])
fig, ax = plt.subplots(figsize=(18, 8))
#ax.axhline(62000, color="k", lw=0.5)
#ax.axhline(75000, color="k", lw=0.5)
#ax.axvline(0.0, color="k", lw=0.5)
scat = ax.scatter(X, TT,
         c=LL, marker="o", s=150, label="The best CLOUDY models")

ax.scatter(X_TT_min[:,0], X_TT_min[:,1], s=280, lw=1.5, facecolors='none', edgecolors='r')

plt.text(0.7, 0.95, 'CLOUDY models', horizontalalignment='left',
         verticalalignment='top', fontsize = 25, transform=ax.transAxes)

ax.annotate(r"$\mathrm{T_{eff}= 1.4\times10^5}K,~ \mathrm{L=1.6\times10^3L_{\odot}}, ~ \chi^2_{\mathrm{min}}=18.10$", 
    xy=(X_TT_min[:,0], X_TT_min[:,1]), xycoords='data', color='red', fontsize=16.5, zorder=100,
      xytext=(0.2, 0.52), textcoords='axes fraction',
      arrowprops=dict(arrowstyle="->",
                      connectionstyle="arc3,rad=0.2",  color='red',
                      ))
    
#ax.legend()
ax.set(
xlabel="$\chi^2$",
    ylabel="$\log_{10}\, T_{\mathrm{eff}}$")
#xlim=[-10000, 10000],
#ylim=[3000, 3e5],
#              )
# ax4.set_xscale("symlog")#, linthreshx=100)
# ax4.set_yscale("log")
#cb = plt.colorbar(ax=ax)
cb = fig.colorbar(scat, extend='both', ax=ax)#
cb.set_label("$\log_{10}\, L/L_\odot$", fontsize=25)
cb.ax.tick_params(labelsize=25)
#plt.gca().invert_xaxis()
sns.despine()
plt.tight_layout()
fig.savefig("../Text/Figs/chi-temperature.pdf")
            
            
    
