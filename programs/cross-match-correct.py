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
tab_base = Table.from_pandas(df)                                                      #
                                                                                   #
                                                                                   #
tab_base['coord'] = coord.SkyCoord(ra=tab_base['DRAJ2000'], dec=tab_base['DDECJ2000'],      #
                                 unit=('deg', 'deg'))                              #
####################################################################################

cmd_args = parser.parse_args()
file_ = cmd_args.table + ".csv"

# Reading the
truedf = pd.read_csv(file_)
tab = Table.from_pandas(truedf)
tab['coord'] = coord.SkyCoord(ra=tab['RA'], dec=tab['DEC'],
                                      unit=('deg', 'deg'))
MAXSEP = 2.0
sourceseps = []
for pn in tab_base:
    seps = pn['coord'].separation(tab['coord']).arcsec
    iclosest = seps.argmin()
    sourceseps.append(seps[iclosest])
    

#mask
mask = np.array(sourceseps) <= MAXSEP
tab_match = tab_base[mask]

# Create an empty list to store duplicate indices
duplicate_indices = []

# Removing duplicate
# Compare each pair of objects to find duplicates
for i in range(len(tab_match['coord'])):
    for j in range(i + 1, len(tab_match['coord'])):
        if tab_match['coord'][i].separation(tab_match['coord'][j]).arcsecond < 1.5:  # Adjust the threshold as needed
            duplicate_indices.append((i, j))

# Print the duplicate object pairs
for i, j in duplicate_indices:
    print(f'Duplicate Objects: Object {tab_match["Name"][i]} and Object {tab_match["Name"][j]}')
    

# Create a new table without duplicates
tab_match_ = tab_match.copy()
tab_match_.remove_rows([index[1] for index in duplicate_indices])

#Shorted by RA
tab_match_.sort('DRAJ2000')
tab.sort('RA')

# Selecting columns in PN
col = ["Name", "PNstat", "MajDiam", "MinDiam"]
matched_pntab_seleCol = tab_match_[col]
#print(matched_pntab_seleCol)

# Final table
truetab_final = hstack([matched_pntab_seleCol, tab])
print(truetab_final)

# Name of the column you want to remove
column_to_remove = 'coord'

# Remove the column if it exists, without checking
truetab_final.remove_column(column_to_remove)

# Save
truetab_final.write(file_.replace(".csv", "-final.ecsv"), format="ascii.ecsv", overwrite=True)
