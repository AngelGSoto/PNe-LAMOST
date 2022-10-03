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
            'PN': {'file': 'PSDR1_tap/PN-PS1dr1.csv',
                     'ID': 'OName'},
            'CV': {'file': 'PSDR1_tap/CV-PS1dr1.csv',
                    'ID': 'recno'},
            'SNR': {'file': 'PSDR1_tap/snr-PS1dr1.csv',
                    'ID': 'SNR'},
            'YSO': {'file': 'PSDR1_tap/yso-rigliaco-PS1dr1.csv',
                    'ID': 'Name'},
            'AeBe': {'file': 'PSDR1_tap/AeBe-PS1dr1.csv',
                    'ID': 'Object'},
            }

catalogs_gaia = {'SySt': {'file': 'match_gaia_tap/SySt_gaia.csv',
                    'ID': 'Name'},
            'PN': {'file': 'match_gaia_tap/PN_gaia.csv',
                     'ID': 'OName'},
            'CV': {'file': 'match_gaia_tap/cv_gaia.csv',
                    'ID': 'recno'},
            'SNR': {'file': 'match_gaia_tap/snr_gaia.csv',
                    'ID': 'SNR'},
            'YSO': {'file': 'match_gaia_tap/yso_gaia-rialco.csv',
                    'ID': 'Name'},
            'AeBe': {'file': 'match_gaia_tap/AeBe_gaia.csv',
                    'ID': 'Object'},
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
    fileascii = "{}-PS-Gaia.ecsv".format(cat_ps) 
    tab_final.write(fileascii, format="ascii.ecsv", overwrite=True)

    

