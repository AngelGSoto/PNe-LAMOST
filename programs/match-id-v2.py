from astropy.table import Table, hstack
import numpy as np
import argparse
import os
import pandas as pd

from pathlib import Path
  
ROOT1 = Path("PSDR1_tap")
ROOT2 = Path("match_gaia_tap")

catalogs_ps = {'SySt': {'file': 'PSDR1_tap/syst-PS1dr1.csv',
                    'ID': 'Name'},
            'PN': {'file': 'PSDR1_tap/PN-PS1dr1.csv',
                     'ID': 'OName'},
            'CV': {'file': 'PSDR1_tap/CV-PS1dr1.csv',
                    'ID': 'MAX'},
            'SNR': {'file': 'PSDR1_tap/snr-PS1dr1.csv',
                    'ID': 'SNR'},
            'YSO': {'file': 'PSDR1_tap/yso-rigliaco-PS1dr1.csv',
                    'ID': 'Name'},
            'AeBe': {'file': 'PSDR1_tap/AeBe-PS1dr1.csv',
                    'ID': 'Object'},
            }

catalogs_gaia = {'SySt': {'file': 'PSDR1_tap/syst-PS1dr1.csv',
                    'ID': 'Name'},
            'PN': {'file': 'PSDR1_tap/PN-PS1dr1.csv',
                     'ID': 'OName'},
            'CV': {'file': 'PSDR1_tap/CV-PS1dr1.csv',
                    'ID': 'MAX'},
            'SNR': {'file': 'PSDR1_tap/snr-PS1dr1.csv',
                    'ID': 'SNR'},
            'YSO': {'file': 'PSDR1_tap/yso-rigliaco-PS1dr1.csv',
                    'ID': 'Name'},
            'AeBe': {'file': 'PSDR1_tap/AeBe-PS1dr1.csv',
                    'ID': 'Object'},
            }

df = pd.read_csv(ROOT1 / file1)

# Table with the IDs
df2 = pd.read_csv(ROOT2 / file2)

# Converting pandas in astropy table
tab1 = Table.from_pandas(df)
tab2 = Table.from_pandas(df2)

# Making mask and applying
id1 = tab1["Name"]
id2 = tab2["Name"]
mask1 = np.array([source in id2 for source in id1])
mask2 = np.array([source in id1 for source in id2])

# aplying the mask
tab_new1 = tab1[mask1]
tab_new2 = tab2[mask2]

# Selecting columns of Gaia table
col = ["Plx", "e_Plx", "Gmag", "e_Gmag", "BPmag", "e_BPmag", "RPmag", "e_RPmag", "BP-RP", "BP-G", "G-RP", "Teff"]

tab_new2_ = tab_new2[col]
tab_final = hstack([tab_new1, tab_new2_])

print("Number of objects is:", len(tab_final))
#df_new2 = df2[mask]

#Save the final file (ASCII)
tab_final.write("SySt-PS-Gaia.ecsv", format="ascii.ecsv")
