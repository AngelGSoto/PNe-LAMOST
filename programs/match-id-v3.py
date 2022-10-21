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

catalogs_ps = {'SySt': {'file': 'gaia_distance_tab/SySt-dist-gaia.csv',
                    'ID': 'Name'},
            'PN': {'file': 'gaia_distance_tab/PN-santamaria-dist-gaia.csv',
                     'ID': 'PNG'},
            'CV': {'file': 'gaia_distance_tab/CV-dist-gaia.csv',
                    'ID': 'recno'},
            'SNR': {'file': 'gaia_distance_tab/snrs-dist-gaia.csv',
                    'ID': 'SNR'},
            'YSO': {'file': 'gaia_distance_tab/YSO-dist-gaia.csv',
                    'ID': 'Name'},
            'AeBe': {'file': 'gaia_distance_tab/AeBe-dist-gaia.csv',
                    'ID': 'Object'},
            'Star': {'file': 'gaia_distance_tab/Stars-MNRAS-dist-gaia.csv',
                    'ID': 'Name'},
            'Star2': {'file': 'gaia_distance_tab/Star-smart-dist-gaia.csv',
                    'ID': 'GaiaEDR3'},
            'cans-hue': {'file': 'gaia_distance_tab/cans-hue-dist-gaia.csv',
                    'ID': 'LAMOST'},
            'cans-new': {'file': 'gaia_distance_tab/cans-new-dist-gaia.csv',
                    'ID': 'LAMOST'},
            'cans-sim': {'file': 'gaia_distance_tab/cans-sim-dist-gaia.csv',
                    'ID': 'LAMOST'},
            }

catalogs_gaia = {'SySt': {'file': 'match_gaia_tap/SySt_gaia.csv',
                    'ID': 'Name'},
            'PN': {'file': 'match_gaia_tap/PN-santamaria_gaia.csv',
                     'ID': 'PNG'},
            'CV': {'file': 'match_gaia_tap/cv_gaia.csv',
                    'ID': 'recno'},
            'SNR': {'file': 'match_gaia_tap/snr_gaia.csv',
                    'ID': 'SNR'},
            'YSO': {'file': 'match_gaia_tap/yso_gaia-rialco.csv',
                    'ID': 'Name'},
            'AeBe': {'file': 'match_gaia_tap/AeBe_gaia.csv',
                    'ID': 'Object'},
            'Star': {'file': 'match_gaia_tap/stars_gaia.csv',
                     'ID': 'Name'},
            'Star2': {'file': 'match_gaia_tap/stars-Smart_gaia.csv',
                     'ID': 'GaiaEDR3'},
            'cans-hue': {'file': 'match_gaia_tap/cans-hue_gaia.csv',
                     'ID': 'LAMOST'},
            'cans-new': {'file': 'match_gaia_tap/cans-new_gaia.csv',
                     'ID': 'LAMOST'},
            'cans-sim': {'file': 'match_gaia_tap/cans-sim_gaia.csv',
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
    fileascii = "{}-Gaia-distance.ecsv".format(cat_ps) 
    tab_final.write(fileascii, format="ascii.ecsv", overwrite=True)
