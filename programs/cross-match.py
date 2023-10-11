import numpy as np
from astropy.table import Table, hstack
import astropy.coordinates as coord
import pandas as pd
from astropy import units as u
import argparse
import sys


parser = argparse.ArgumentParser(
    description="""Reading the ouput cloudy models""")

parser.add_argument("table", type=str,
                    default="DR7lamost-True-PN-duplicate-GAIA-PS1",
                    help="Name of input model ")

####################################################################################
df = pd.read_csv('../Luis_arzd099bej.csv')                                         #
#                    fill_values=('--', np.nan) ).filled(np.nan)                   #
pntab = Table.from_pandas(df)                                                      #
                                                                                   #
                                                                                   #
pntab['coord'] = coord.SkyCoord(ra=pntab['DRAJ2000'], dec=pntab['DDECJ2000'],      #
                                 unit=('deg', 'deg'))                              #
####################################################################################

cmd_args = parser.parse_args()
file_ = cmd_args.table + ".csv"

# Reading the
truedf = pd.read_csv(file_)
truetab = Table.from_pandas(truedf)
truetab['coord'] = coord.SkyCoord(ra=truetab['RA'], dec=truetab['DEC'],
                                     unit=('deg', 'deg'))

# Perform the cross-match
idx, d2d, _ = pntab['coord'].match_to_catalog_sky(truetab['coord'])

# Define the search radius in arcseconds
search_radius = 1.5 * u.arcsecond

# Select matches within the search radius
matches_within_radius = d2d < search_radius

# Get the matched objects from catalog1 and catalog2
matched_pntab = pntab[matches_within_radius]

matched_pntab.sort('DRAJ2000')
truetab.sort('RA')

# Selecting columns in PN
col = ["Name", "MajDiam", "MinDiam"]
matched_pntab_seleCol = matched_pntab[col]

# Final table
truetab_final = hstack([matched_pntab_seleCol, truetab])

# Name of the column you want to remove
column_to_remove = 'coord'

# Remove the column if it exists, without checking
truetab_final.remove_column(column_to_remove)

# Save
truetab_final.write(file_.replace(".csv", "-final.ecsv"), format="ascii.ecsv", overwrite=True)
