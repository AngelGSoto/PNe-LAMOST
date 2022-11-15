'''
Autor: L. A. Gutiérrez Soto
13/11/2022
'''

from __future__ import print_function
import json
from astropy.io import fits
from astropy import coordinates as coord
from astropy import units as u
import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import aplpy
import pyregion
import sys
from astropy.coordinates import SkyCoord
# from misc_utils import run_info, update_json_file
# from fits_utils import get_image_hdu

parser = argparse.ArgumentParser(
    description="""Plot side-by-side pure image of source and image with overlays""")

parser.add_argument("image", type=str,
                    default="1000001-JPLUS-00873-v2_rSDSS_swp-crop",
                    help="Name of original FITS image (section in database)")

parser.add_argument("--vmin", type=float, default=None,
                    help="""Set minimum brightness directly - overrides minfactor""")
parser.add_argument("--vmax", type=float, default=None,
                    help="""Set maximum brightness directly - overrides maxfactor""")

parser.add_argument("--debug", action="store_true",
                    help="Print out verbose debugging info")

MAXFACTOR_DEFAULT = 3.0
MINFACTOR_DEFAULT = 3.0
cmd_args = parser.parse_args()

# Maximum
vmax = cmd_args.vmax

# Minimum
vmin = cmd_args.vmin

image_name = cmd_args.image + ".fits"

hdul = fits.open(image_name)
instrument = hdul[0].header['BANDPASS']

                      
# Plot image of the FITS array of this object
# 
plt.clf()
f = plt.figure(figsize=(18,9))

ax1 = aplpy.FITSFigure(hdul, figure=f, subplot=(1, 1, 1))#, north=True)
#ax1.recenter(ra0, dec0, 4*R0/cmd_args.zoom)
ax1.show_grayscale(invert=True, vmin=vmin, vmax=vmax)
           

# With the mask regions, the file may not exist
# position = "position.reg"
# ra, dec = [], []

# f = open(position, 'r')
# header1 = f.readline()
# header2 = f.readline()
# header3 = f.readline()
# for line in f:
#     line = line.strip()
#     columns = line.split()
#     coor = line.split("(")[-1].split("\"")[0]
#     ra1, dec1 = coor.split(",")[0:2]
#     c = SkyCoord(ra1, dec1, unit=(u.hourangle, u.deg))
#     ra.append(c.ra.degree)
#     dec.append(c.dec.degree)

ax1.axis_labels.set_xtext('RA (J2000)')
#img.axis_labels.hide_x()
ax1.axis_labels.set_ytext('Dec (J2000)')
ax1.axis_labels.set_font(size=18, weight='medium', stretch='normal', family='sans-serif', style='normal', variant='normal')
ax1.tick_labels.set_font(size=18, weight='medium', stretch='normal', family='sans-serif', style='normal', variant='normal')
#ax1.recenter(ra, dec, width=0.05, height=0.05)    
ax1.add_scalebar(30.0/3600)
ax1.scalebar.set_label('30"')
ax1.scalebar.set(color='black', linewidth=4, alpha=0.9)
ax1.scalebar.set_font(size=45, weight='bold',
                      stretch='normal', family='sans-serif',
                      style='normal', variant='normal')

ax1.add_label(0.1, 0.9, "Red", color="black",
              horizontalalignment='left',
              weight='bold', size=30, relative=True, zorder=1000)
dx, dy = 0.001, -0.001
ax1.add_label(0.7+dx, 0.89+dy, "PRTM 1", color="black", alpha=0.9,
              horizontalalignment='left',
              bbox={"facecolor": "black", "edgecolor": "none",# "pad": 20,
                    "alpha": 0.1, "boxstyle": "round, pad=0.5"},
              weight='bold', size=18, relative=True, zorder=999)

# #ax1.list_layers()
# ax1.show_markers(ra, dec, layer='marker', edgecolor='green', facecolor='none', marker='o', s=10, alpha=0.9, linewidths=60)#, layer='marker_set_1', edgecolor='black', facecolor='none', s=30, alpha=0.5, linewidths=20)


# # ax1.axis_labels.hide_y()
# # ax1.tick_labels.hide_y()
# ax1.add_colorbar()
# ax1.colorbar.set_axis_label_text("counts")
# ax1.colorbar.set_location("top")
# #ax2.colorbar.set_box([0.95, 0.1, 0.015, 0.8])
ax1.set_theme('publication')
# #f.tight_layout()
# #f.savefig("-".join([image_name, "images.pdf"]))

plt.savefig(image_name.replace(".fits", ".pdf"))

# # record the --maxfactor and the --minfactor in the *-xycb.json file
# # also their respective help section
