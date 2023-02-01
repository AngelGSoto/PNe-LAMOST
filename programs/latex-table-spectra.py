'''
Create file.tex with several figures (table).
'''
from __future__ import print_function
import numpy as np
from astropy.io import fits
import os
import glob
import json
import matplotlib.pyplot as plt
import pandas as pd
#import StringIO
from astropy.table import Table
import seaborn as sns
import sys
from scipy.optimize import fsolve
import colours


#Read de files
pattern = "*.pdf"
file_list = glob.glob(pattern)
best_model = ["model_140000_37.15_3.70.pdf", "model_150000_36.98_3.60.pdf", "model_140000_37.25_3.78.pdf"]

# Discard the best models
main_files = list(set(file_list) - set(best_model))

# Number of objects
nrows = len(main_files)

fig_template = r'\includegraphics[width=0.24\linewidth, clip]{{{:s}}}'

NCOLS = 4
thefiles = []
thistable = []

for filename in main_files:
    thefiles.append(fig_template.format("Figs/" + filename))
    if len(thefiles) == NCOLS:
        thistable.append(thefiles)
        thefiles = []

if thefiles:
    # catch any partial row at the end
    thistable.append(thefiles)
             
def output_row(row):
    return " & ".join(row) + r" \\"

def output_figtable(table):
    result = r"\begin{table*}" + "\n"
    result += r"\centering" + "\n"
    result += r"  \caption{The best fitted models, on which the $\chi^{2}_r$ have value between 1 and 2." + " " + "\label{tab:best-model12}}" + "\\" + "\n"
    result += r"  \begin{tabular}{" + "l "*NCOLS + "}" + "\n"
    for row in table:
        result += "    " + output_row(row) + "\n"
    result += r"  \end{tabular}" + "\n"
    result += r"\end{table*}" + "\n"
    return result

with open("../Text/models-spectra.tex", 'w') as f:
    f.write(output_figtable(thistable))
