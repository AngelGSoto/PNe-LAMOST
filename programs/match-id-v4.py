'''
Luis A. Gutiérrez Soto
This script allows to find the coincidences between list using the ID
'''
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
            'CV': {'file': 'PSDR1_tap/CV-PS1dr1.csv',
                    'ID': 'recno'},
            'SNR': {'file': 'PSDR1_tap/snr-PS1dr1.csv',
                    'ID': 'SNR'},
            'YSO': {'file': 'PSDR1_tap/yso-rigliaco-PS1dr1.csv',
                    'ID': 'Name'},
            'AeBe': {'file': 'PSDR1_tap/AeBe-PS1dr1.csv',
                    'ID': 'Object'},
            'cans-hue': {'file': 'PSDR1_tap/cans-hue-PS1dr1.csv',
                    'ID': 'LAMOST'},
            'cans-new': {'file': 'PSDR1_tap/cans-new-PS1dr1.csv',
                    'ID': 'LAMOST'},
            'cans-sim': {'file': 'PSDR1_tap/cans-sim-PS1dr1.csv',
                    'ID': 'LAMOST'},
            }

catalogs_gaia = {'SySt': {'file': 'gaiaEDR3-tap/SySt-gaiadr3.csv',
                    'ID': 'Name'},
            'CV': {'file': 'gaiaEDR3-tap/CV-gaiadr3.csv',
                    'ID': 'recno'},
            'SNR': {'file': 'gaiaEDR3-tap/snrs-gaiadr3.csv',
                    'ID': 'SNR'},
            'YSO': {'file': 'gaiaEDR3-tap/YSO-rialco-gaiadr3.csv',
                    'ID': 'Name'},
            'AeBe': {'file': 'gaiaEDR3-tap/AeBe-gaiadr3.csv',
                    'ID': 'Object'},
            'cans-hue': {'file': 'gaiaEDR3-tap/cans-hue-gaiadr3.csv',
                     'ID': 'LAMOST'},
            'cans-new': {'file': 'gaiaEDR3-tap/cans-new-gaiadr3.csv',
                     'ID': 'LAMOST'},
            'cans-sim': {'file': 'gaiaEDR3-tap/cans-sim-gaiadr3.csv',
                     'ID': 'LAMOST'},
            }

for (cat_ps, metadata_ps), (cat_gaia, metadata_gaia) in zip(catalogs_ps.items(),
                                                            catalogs_gaia.items()):
    
    df_ps = pd.read_csv(metadata_ps["file"])
    df_gaia = pd.read_csv(metadata_gaia["file"])
    # # Converting pandas in astropy table
    tab_ps = Table.from_pandas(df_ps)
    tab_gaia = Table.from_pandas(df_gaia)

    # Making mask and applying
    id1 = tab_ps[metadata_ps["ID"]]
    id2 = tab_gaia[metadata_gaia["ID"]]
    
    mask1 = np.array([source in id2 for source in id1])
    mask2 = np.array([source in id1 for source in id2])

    # aplying the mask
    tab_new1 = tab_ps[mask1]
    tab_new2 = tab_gaia[mask2]
    
    tab_final = hstack([tab_new1, tab_new2])

    #Save the final file (ASCII.ECSV)
    fileascii = "{}-PS-GaiaEDR3.ecsv".format(cat_ps) 
    tab_final.write(fileascii, format="ascii.ecsv", overwrite=True)
