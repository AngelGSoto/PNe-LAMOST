'''
Adapted from: https://github.com/will-henney/teresa-turtle/blob/master/cspn-notebook/Turtle%20CSPN%20Miller-Bertolami%20models.py
'''
import numpy as np
from astropy.table import Table
from pathlib import Path

datadir = Path("Model-Bertolami")

byte_by_byte_description = """
Byte-by-byte Description of file:
--------------------------------------------------------------------------------
   Bytes Format Units     Label       Explanations
--------------------------------------------------------------------------------
   1-  5  I5    ---       N           Track point number
   7- 15  F9.6  [Lsun]    logL        logarithm of the stellar luminosity
  17- 25  F9.6  [K]       logTeff     logarithm of the effective temperature
  27- 35  F9.6  [cm/s2]   logg        logarithm of the surface gravity
  40- 51  F12.4 yr        t           Age since the point at LogTeff=3.85
  53- 61  F9.6  ---       Menv        Fractional mass of the envelope
  63- 71  F9.6  Msun      Mstar       Total mass of the star
  73- 82  F10.6 [Msun/yr] log(-dM/dt)  Logarithm of the Mass Loss Rate,
                                       log(-dMstar/dt)
--------------------------------------------------------------------------------
"""

def read_tracks(datafile):
    """Read each Miller–Bertolami track into a separate astropy.table
    
    Input argument `datafile` is a CDS file containing all tracks 
    for a given metallicity, e.g., "0100_t03.dat"
    
    Returns list of tables. Each table has a metadata "comments" field 
    that contains additional info (mass and surface composition). 
    """
    with open(datafile) as f:
        # Each track is separated by two blank lines
        tracks = f.read().split("\n\n\n")[:-1]
        tables = []
        for track in tracks:
            lines = track.split("\n")
            metadata = lines[:6]
            data = lines[7:]
            datastring = "\n".join(
                [byte_by_byte_description] + data
            )
            table = Table.read(datastring, format="ascii.cds")
            table.meta["comments"] = metadata
            tables.append(table)
    return tables

tabs = read_tracks(datadir / "0100_t03.dat")

from matplotlib import pyplot as plt
import seaborn as sns
# %matplotlib inline
sns.set_context("talk")
sns.set_color_codes()


#Plot of effective temperature versus gravity, which we compare with the turtle observed values.

def extract_masses(data):
    _, Mi, Mf, _ = data.meta["comments"][2].split()
    return round(float(Mi), 2), round(float(Mf), 3)


fig, ax = plt.subplots(figsize=(8, 8))
ax.axhline(4.875, color="k", lw=0.5)
ax.axhline(4.792, color="k", lw=0.5)
ax.axvspan(4.6, 5.0, color="k", alpha=0.1)
ax.axvline(4.8, color="k", lw=0.5)
for data in tabs:
    try:
        Mi, Mf = extract_masses(data)
        label = f"({Mi}, {Mf})"
    except:
        continue
    ax.plot(
        "logg", "logTeff",
        data=data, label=label,
    )
ax.legend()
ax.set(
    xlabel="$\log_{10}\, g$",
    ylabel="$\log_{10}\, T_{\mathrm{eff}}$",
    xlim=[4.0, 6.0],
    ylim=[4.7, 5.0],
)
sns.despine()
fig.savefig("hr-g-T.pdf")

# From this we conclude that the $2\,M_\odot$ model is the best fit, but we cannot rule out $1.25\,M_\odot$ to $3\,M_\odot$.

# Next, look at the timescales.

fig, [axL, axMd, ax] = plt.subplots(3, 1, sharex=True, figsize=(12, 8))
ax.axhline(62000, color="k", lw=0.5)
ax.axhline(75000, color="k", lw=0.5)
ax.axhline(20000, color="r", ls="--", lw=0.8)
for axx in axL, axMd, ax:
    axx.axvline(0.0, color="k", lw=0.5)
    axx.axvspan(-100, 100, color="m", alpha=0.05)
for data in tabs:
    try:
        Mi, Mf = extract_masses(data)
        label = f"({Mi}, {Mf})"
    except:
        continue
    data["Teff"] = 10**data["logTeff"]
    data["L"] = 10**data["logL"]
    data["Mdot"] = 10**data["log(-dM/dt)"]
    ax.plot(
        "t", "Teff",
        data=data, label=label,
    )
    axMd.plot(
        "t", "Mdot",
        data=data, label=label,
    )
    axL.plot(
        "t", "L",
        data=data,
    )
ax.legend()
ax.set(
    xlabel="time, years",
    ylabel="$T_{\mathrm{eff}}$, K",
    xlim=[-100000, 100000],
    ylim=[3000, 3e5],
)
ax.set_xscale("symlog")#, linthreshx=100)
ax.set_yscale("log")
axL.set(
    ylim=[0.0, 1.5e4],
    ylabel="$L / L_\odot$",
)
axMd.set(
    ylim=[1e-8, 1e-5],
    ylabel="$d M / dt$, $M_\odot$/yr",
    yscale="log",
)
sns.despine()
fig.savefig("time.pdf")

# So, conclusion from this is that we need the mass to be 1.25 to 1.5.  If it is 1.25, then the MB models have it not going through a C-star phase, which is needed to explain the low C/O ratio in nebula.
#
# On the other hand, if it was a triple interaction that ejected the envelope, it might have happened before the AGB got to end of its evolution. 

def make_table_of_times(tabs, Teff):
    logTeff = np.log10(Teff)
    tTkey = f"t({Teff})"
    rslts = {
        "Mi": [],
        "Mf": [],
        tTkey: [],
        "t_cross": [],
        "t_tr": [],
    }
    for data in tabs:
        Mi, Mf = extract_masses(data)
        rslts["Mi"].append(Mi)
        rslts["Mf"].append(Mf)
        # First time to reach given Teff
        mask = data["logTeff"] >= logTeff
        tT = data[mask]["t"].min()
        rslts[tTkey].append(tT)
        # Time to cross to maximum Teff
        icross = data["logTeff"].argmax()
        rslts["t_cross"].append(data[icross]["t"])
        # Transition time before t = 0
        rslts["t_tr"].append(-data["t"].min())
    return Table(rslts)


times = make_table_of_times(tabs, 75000)

print("times=", times)

# ## HR diagram with the high ionization nebula
weidmantab = Table.read("../PN-high-ion-weidman.dat", format="ascii")

###
# Teff and L of our PN after modeling with Cloudy.
logTeffpn = [np.log10(13e4), np.log10(14e4), np.log10(16e4)]
logLpn = [np.log10(8.44579999999994e+36/3.839e33), np.log10(7.29410000000006e+36/3.839e33), np.log10(6.142400000000004e+36/3.839e33)]  #36.86297171307024, 3.7238299999999735e+37
models_cloudy = ["Model 3", "Model 2", "Model 1"]

fig, ax = plt.subplots(figsize=(8, 8))
#ax.axvspan(4.7, 5.0, 0.6, 0.9, color="k", alpha=0.1)
lw = 0.5
tkin = 3500.0
logTion = 4.301029995663981
logTtkin = []
logLtkin = []
for data in tabs:
    try:
        Mi, Mf = extract_masses(data)
        label = f"({Mi}, {Mf})"
        labeli = f"({Mi})"
    except:
        continue
    ax.plot(
        "logTeff", "logL",
        data=data, label=label,
        zorder=-100, c="k", lw=lw,
    )
    t0 = np.interp(logTion, data["logTeff"], data["t"])
    logT = np.interp(tkin + t0, data["t"], data["logTeff"])
    logL = np.interp(tkin + t0, data["t"], data["logL"])
    logTtkin.append(logT)
    logLtkin.append(logL)
    #ax.plot(logT, logL, "*", c="k")
    m = (data["t"] > t0 + tkin/1.5) & (data["t"] < t0 + tkin*1.5)
    # ax.plot(
    #     "logTeff", "logL",
    #     data=data[m], label="_nolabel_",
    #     zorder=-100, c="b", lw=7, alpha=0.4,
    # )
    
    lw += 0.2
    Tmin = data["logTeff"].min()
    Lmax = data["logL"].max()

    bbox_props = dict(boxstyle="round", fc="w", ec="0.88", alpha=0.1, pad=0.1)
    # for label_, x, y in zip(labeli, Tmin, Lmax):
    ax.annotate(r"$M = $" + labeli.split("(")[-1].split(")")[0], (Tmin, Lmax), alpha=1, size=13,
                xytext=(-15, -15), textcoords='offset points', ha='right', va='bottom', rotation=0, bbox=bbox_props, zorder=200)


def my_annotate(ax, s, xy_arr=[], *args, **kwargs):
    '''
    Script to make several arrows
    to match only one taxt. 
    Taken from https://stackoverflow.com/questions/14542232/annotate-several-points-with-one-text-in-matplotlib
    '''
    ans = []
    an = ax.annotate(s, xy_arr[0], *args, **kwargs)
    ans.append(an)
    d = {}
    try:
        d['xycoords'] = kwargs['xycoords']
    except KeyError:
        pass
    try:
        d['arrowprops'] = kwargs['arrowprops']
    except KeyError:
        pass
    for xy in xy_arr[1:]:
        an = ax.annotate(s, xy, alpha=0.0, xytext=(0,0), textcoords=an, **d)
        ans.append(an)
    return ans
    
# my_annotate(ax, r"$\mathrm{t_{evo} = 3500~yr \pm 50\%}$", xy_arr=[(logTtkin[0], logLtkin[0]), (logTtkin[1], logLtkin[1]), (logTtkin[2], logLtkin[2])], size=15, xycoords='data',
#             xytext=(-50, -100), textcoords='offset points',
#             arrowprops=dict(arrowstyle="->",
#                              connectionstyle="arc3,rad=-0.2"),
#                 )
ax.scatter(
    "logTeff", "logL", data=weidmantab, 
    color="#ff7f0e", marker="o", s=130, label="Triple",
    edgecolors="k")
bbox_props1 = dict(boxstyle="round", fc="w", ec="0.78", alpha=0.6, pad=0.2)
for label_, x, y in zip(weidmantab["Name_2"], weidmantab["logTeff"], weidmantab["logL"]):
    ax.annotate(label_, (x, y), alpha=1, size=10,
                   xytext=(45.0, -17), textcoords='offset points', ha='right', va='bottom', bbox=bbox_props1, zorder=-100)

#Our PN

ax.scatter(logTeffpn, logLpn, 
    color="#377eb8", marker="*", s=550, label="Triple",
    edgecolors="k")

for a, b, c in zip(models_cloudy, logTeffpn, logLpn):
    print(a, b, c)
    ax.annotate(a, (b, c), alpha=1, size=10,
                       xytext=(-10.0, 8), textcoords='offset points', ha='right', va='bottom', bbox=bbox_props1, zorder=-100)
#ax.legend()
ax.set(
    ylabel="$\log_{10}\, L/L_\odot$",
    xlabel="$\log_{10}\, T_{\mathrm{eff}}$",
    xlim=[3.58, 5.5],
    ylim=[2.0, None]
    ,
)

plt.gca().invert_xaxis()
sns.despine()
fig.savefig("hr-planetarieNebula.pdf")
None
# Now, we try the appendix B tracks (how are they different?)

tabs_tb2 = read_tracks(datadir / "0100_tb2.dat")

tabs_tb2[0].meta

fig1, ax1 = plt.subplots(figsize=(8, 8))
ax1.axhline(4.875, color="k", lw=0.5)
ax1.axhline(4.792, color="k", lw=0.5)
ax1.axvspan(4.6, 5.0, color="k", alpha=0.1)
ax1.axvline(4.8, color="k", lw=0.5)
for data in tabs_tb2:
    try:
        _, Mi, Mf, _ = data.meta["comments"][2].split()
        Mi = round(float(Mi), 2)
        Mf = round(float(Mf), 3)
        label = f"({Mi}, {Mf})"
    except:
        continue
    ax1.plot(
        "logg", "logTeff",
        data=data, label=label,
    )
ax1.legend()
ax1.set(
    xlabel="$\log_{10}\, g$",
    ylabel="$\log_{10}\, T_{\mathrm{eff}}$",
    xlim=[4.0, 6.0],
    ylim=[4.7, 5.0],
)
sns.despine()
fig1.savefig("hr-gT-0100tb2.pdf")
None

fig2, ax2 = plt.subplots(figsize=(8, 8))
ax2.axhline(62000, color="k", lw=0.5)
ax2.axhline(75000, color="k", lw=0.5)
ax2.axvline(0.0, color="k", lw=0.5)
for data in tabs_tb2:
    try:
        _, Mi, Mf, _ = data.meta["comments"][2].split()
        Mi = round(float(Mi), 2)
        Mf = round(float(Mf), 3)
        label = f"({Mi}, {Mf})"
    except:
        continue
    data["Teff"] = 10**data["logTeff"]
    ax2.plot(
        "t", "Teff",
        data=data, label=label,
    )
ax2.legend()
ax2.set(
    xlabel="time, years",
    ylabel="$T_{\mathrm{eff}}$, K",
    xlim=[-10000, 10000],
    ylim=[3000, 3e5],
)
ax2.set_xscale("symlog")#, linthreshx=100)
ax2.set_yscale("log")
sns.despine()
fig2.savefig("hr-tT-0100tb2.pdf")
None
# Now, the low metallicity tracks

tabs_Z0010 = read_tracks(datadir / "0010_t03.dat")

fig3, ax3 = plt.subplots(figsize=(8, 8))
ax3.axhline(4.875, color="k", lw=0.5)
ax3.axhline(4.792, color="k", lw=0.5)
ax3.axvspan(4.6, 5.0, color="k", alpha=0.1)
ax3.axvline(4.8, color="k", lw=0.5)
for data in tabs_Z0010:
    try:
        _, Mi, Mf, _ = data.meta["comments"][2].split()
        Mi = round(float(Mi), 2)
        Mf = round(float(Mf), 3)
        label = f"({Mi}, {Mf})"
    except:
        continue
    ax3.plot(
        "logg", "logTeff",
        data=data, label=label,
    )
ax3.legend()
ax3.set(
    xlabel="$\log_{10}\, g$",
    ylabel="$\log_{10}\, T_{\mathrm{eff}}$",
    xlim=[4.0, 6.0],
    ylim=[4.7, 5.0],
)
sns.despine()
fig3.savefig("hr-logglogT-0100tb2.pdf")
None
None

fig4, ax4 = plt.subplots(figsize=(8, 8))
ax4.axhline(62000, color="k", lw=0.5)
ax4.axhline(75000, color="k", lw=0.5)
ax4.axvline(0.0, color="k", lw=0.5)
for data in tabs_Z0010:
    try:
        _, Mi, Mf, _ = data.meta["comments"][2].split()
        Mi = round(float(Mi), 2)
        Mf = round(float(Mf), 3)
        label = f"({Mi}, {Mf})"
    except:
        continue
    data["Teff"] = 10**data["logTeff"]
    ax4.plot(
        "t", "Teff",
        data=data, label=label,
    )
ax4.legend()
ax4.set(
    xlabel="time, years",
    ylabel="$T_{\mathrm{eff}}$, K",
    xlim=[-10000, 10000],
    ylim=[3000, 3e5],
)
ax4.set_xscale("symlog")#, linthreshx=100)
ax4.set_yscale("log")
sns.despine()
fig4.savefig("hr-tTeff-0100tb2.pdf")


