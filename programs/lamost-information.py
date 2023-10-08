'''
Getting things from the lamost file spectra.
'''
from astropy.io import ascii
from astropy.table import Table
import numpy as np
import argparse
from astropy.io import fits
import matplotlib.pyplot as plt
import seaborn as sn
import glob

patterns = "*.fits"
files_ = glob.glob(patterns)

filename = []
RA = []
DEC = []
# Print the list of matching files
for files in files_:
    hdu = fits.open(files)
    # Get the header of the primary data unit (usually the first HDU)
    header = hdu[0].header
    filename.append(files)
    RA.append(header["RA"])
    DEC.append(header["DEC"])

t = Table([filename, RA, DEC], names=('FileName', 'RA', 'DEC'), meta={'name': 'first table'})
# Save the table to a CSV file
t.write('DR7lamost-Likely-PN.csv', format='csv', overwrite=True)
