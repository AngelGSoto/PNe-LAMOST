import numpy as np
import json
import matplotlib.pyplot as plt
from  astropy.table import Table
import pandas as pd
import seaborn as sns
import argparse
from pathlib import Path
ROOT_PATH = Path("../PS")

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

#Rename columns to put in PS format
df_ = df_col.rename(columns={"RAJ2000": "ra", "DEJ2000": "dec"})

asciifile = file_.replace(".csv", "-ps-list-match.csv")
df_.to_csv(ROOT_PATH / asciifile, index=False)
