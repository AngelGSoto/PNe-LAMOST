import numpy as np
import json
import matplotlib.pyplot as plt
from  astropy.table import Table
import pandas as pd
import seaborn as sns
import argparse
from pathlib import Path
ROOT_PATH = Path("../gaia_match")

parser = argparse.ArgumentParser(
    description="""Make list format GAIA""")

parser.add_argument("source", type=str,
                    default="teste-program",
                    help="Name of catalog, taken the prefix ")


cmd_args = parser.parse_args()
file_ = cmd_args.source + ".csv"

df = pd.read_csv(file_)

col = ["RAJ2000", "DEJ2000"]
df_col = df[col]

tab = Table.from_pandas(df_col)

#Save the file
asciifile = file_.replace(".csv", "-list-match.txt")
tab.write(ROOT_PATH / asciifile, format='ascii.no_header') 
