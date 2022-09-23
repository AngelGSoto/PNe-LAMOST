from astropy.table import Table, hstack
import numpy as np
import argparse
import os
import pandas as pd

from pathlib import Path
  
ROOT1 = Path("PSDR1_tap")
ROOT2 = Path("match_gaia_tap")

parser = argparse.ArgumentParser(
    description="""Firts table from the S-PLUS catalogs """)

parser.add_argument("table1", type=str,
                    default=" teste-program",
                    help="Name of catalog, taken the prefix ")
parser.add_argument("table2", type=str,
                    default=" teste-program",
                    help="Name of catalog, taken the prefix ")

cmd_args = parser.parse_args()
file1 = cmd_args.table1 + ".csv"

cmd_args = parser.parse_args()
file2 = cmd_args.table2 + ".csv"

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
#df_new2 = df2[mask]


# Save the final file (ASCII)

#df.to_csv(file_df, index=False)

