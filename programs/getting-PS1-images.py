import numpy as np
import time
from astropy.table import Table, hstack
import requests
import ps1bulk
import argparse
import os
import sys

parser = argparse.ArgumentParser(
    description="""Reading the ouput cloudy models""")

parser.add_argument("table", type=str,
                    default="DR7lamost-True-PN-duplicate-GAIA-PS1-final",
                    help="Name of input model ")

parser.add_argument("--Object", type=str,
                    default=None,
                    help="Id object of the source under interest ")

parser.add_argument("--size", type=str,
                    default="200",
                    help="Base name of FITS image that contains the source")

cmd_args = parser.parse_args()
file_ = cmd_args.table + ".ecsv"

tab = Table.read(file_, format="ascii.ecsv")

idfile = tab["FileName"]
    
if cmd_args.Object is not None:
    Object_ = cmd_args.Object
    mask = np.array([source in Object_ for source in idfile])
    tab = tab[mask]
else:
    tab = tab

t0 = time.time()
 
# create a test set of image positions
tdec =  tab["DEC"]
tra = tab["RA"]

#radi = str(5 * (tab["MajDiam"] / 0.25))

# get the PS1 info for those positions
table = ps1bulk.getimages(tra, tdec, size = cmd_args.size, filters="gri")


print("{:.1f} s: got list of {} images for {} positions".format(time.time()-t0,len(table),len(tra)))
 
# if you are extracting images that are close together on the sky,
# sorting by skycell and filter will improve the performance because it takes
# advantage of file system caching on the server
table.sort(['projcell','subcell','filter'])
 
# extract cutout for each position/filter combination
for row in table:
    ra = row['ra']
    dec = row['dec']
    projcell = row['projcell']
    subcell = row['subcell']
    filter = row['filter']
 
    # create a name for the image -- could also include the projection cell or other info
    fname = "t{:08.4f}{:+07.4f}.{}.fits".format(ra,dec,filter)
 
    url = row["url"]
    print("%11.6f %10.6f skycell.%4.4d.%3.3d %s" % (ra, dec, projcell, subcell, fname))
    r = requests.get(url)
    open(fname,"wb").write(r.content)
    
print("{:.1f} s: retrieved {} FITS files for {} positions".format(time.time()-t0,len(table),len(tra)))
