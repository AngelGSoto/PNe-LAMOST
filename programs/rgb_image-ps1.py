'''
Making RGB images from PLUS .fits
Based on original: rgb_image.py
Autor: L. A. Gutiérrez Soto
02/09/20
'''

from __future__ import print_function
import aplpy
import numpy
import sys
from astropy import coordinates as coord
from astropy import units as u
from astropy.coordinates import SkyCoord
import argparse
import matplotlib.pyplot as plt
from astropy.io import fits
import matplotlib
matplotlib.use("Agg")

parser = argparse.ArgumentParser(
    description="""Plot side-by-side RGB images of sources""")

parser.add_argument("image_r", type=str,
                    default="1000001-JPLUS-01485-v2_iSDSS_swp-crop",
                    help="Name of original FITS image (section in database) in i")

parser.add_argument("image_g", type=str,
                    default="1000001-JPLUS-01485-v2_rSDSS_swp-crop",
                    help="Name of original FITS image (section in database) in r")

parser.add_argument("image_b", type=str,
                    default="1000001-JPLUS-01485-v2_gSDSS_swp-crop",
                    help="Name of original FITS image (section in database) in g")

# parser.add_argument("--name", type=str,
#                     default="PSP",
#                     help="Name of the objet")

parser.add_argument("--vmin_r", type=float, default=None,
                    help="""Set minimum brightness directly - overrides minfactor - r""")
parser.add_argument("--vmax_r", type=float, default=None,
                    help="""Set maximum brightness directly - overrides maxfactor - r""")

parser.add_argument("--vmin_g", type=float, default=None,
                    help="""Set minimum brightness directly - overrides minfactor - g""")
parser.add_argument("--vmax_g", type=float, default=None,
                    help="""Set maximum brightness directly - overrides maxfactor - g""")

parser.add_argument("--vmin_b", type=float, default=None,
                    help="""Set minimum brightness directly - overrides minfactor - b""")
parser.add_argument("--vmax_b", type=float, default=None,
                    help="""Set maximum brightness directly - overrides maxfactor - b""")

parser.add_argument("--zoom", type=float, default=None,
                    help="""\Zoom factor to adjust size of plot box - values > 1.0 mean to zoom in""")

parser.add_argument("--position", type=str,
                    default="HYDRA-0026-000010640-position",
                    help="Find the DS9 region")

parser.add_argument("--debug", action="store_true",
                    help="Print out verbose debugging info")


cmd_args = parser.parse_args()
image_r = cmd_args.image_r + ".fits"
image_g = cmd_args.image_g + ".fits"
image_b = cmd_args.image_b + ".fits"

hdul_r = fits.open(image_r)
instrument_r = hdul_r[0].header['HIERARCH FPA.FILTER'].split('.')[0]
hdul_g = fits.open(image_g)
instrument_g = hdul_g[0].header['HIERARCH FPA.FILTER'].split('.')[0]
hdul_b = fits.open(image_b)
instrument_b = hdul_b[0].header['HIERARCH FPA.FILTER'].split('.')[0]

#aplpy.make_rgb_cube(['1000001-JPLUS-01485-v2_iSDSS_swp-crop.fits', '1000001-JPLUS-01485-v2_rSDSS_swp-crop.fits',
                     #'1000001-JPLUS-01485-v2_gSDSS_swp-crop.fits'], 'JPLUS_cube.fits')

aplpy.make_rgb_cube([image_r, image_g, image_b], image_r.replace('.fits', '_cube.fits'))

aplpy.make_rgb_image(image_r.replace('.fits', '_cube.fits'),
                              image_r.replace('.fits', '_rgb.png'),
                      vmin_r=cmd_args.vmin_r, vmax_r=cmd_args.vmax_r, vmin_g=cmd_args.vmin_g,
                                                      vmax_g=cmd_args.vmax_g, vmin_b=cmd_args.vmin_b, vmax_b=cmd_args.vmax_b)

#aplpy.make_rgb_image('JPLUS_cube.fits','JPLUS_linear.png')
#hdul = fits.open('JPLUS_cube_2d.fits')
# aplpy.make_rgb_image('JPLUS_cube.fits','JPLUS_rgb.png',
#                       stretch_r='arcsinh', stretch_g='arcsinh',
#                       stretch_b='arcsinh')


# With the mask regions, the file may not exist
position = cmd_args.position + ".reg"
ra, dec = [], []

try:
    f = open(position, 'r')
    header1 = f.readline()
    header2 = f.readline()
    header3 = f.readline()
    for line in f:
        line = line.strip()
        columns = line.split()
        coor = line.split("(")[-1].split("\"")[0]
        ra1, dec1 = coor.split(",")[0:2]
        c = SkyCoord(ra1, dec1, unit=(u.hourangle, u.deg))
        ra.append(c.ra.hourangle)
        dec.append(c.dec.degree)

except FileNotFoundError:
    print("File", position,
                "not found - is not necesary now")
    
# Launch APLpy figure of 2D cube
img = aplpy.FITSFigure(image_r.replace('.fits', '_cube_2d.fits')) 
img.show_rgb(image_r.replace('.fits', '_rgb.png'))

# Maybe we would like the arcsinh stretched image more?
#img.show_rgb('ic348_color_arcsinh.png')

# Modify the tick labels for precision and format
# img.tick_labels.set_xformat('hh:mm:ss')
# img.tick_labels.set_yformat('dd:mm')
img.axis_labels.set_xtext('RA (J2000)')
#img.axis_labels.hide_x()
img.axis_labels.set_ytext('Dec (J2000)')
img.axis_labels.set_font(size=18, weight='medium', stretch='normal', family='sans-serif', style='normal', variant='normal')
img.axis_labels.hide()
img.axis_labels.hide_y()

img.tick_labels.set_font(size=18, weight='medium', stretch='normal', family='sans-serif', style='normal', variant='normal')
#img.axis_labels.set_yposition('right')
#img.tick_labels.set_yposition('right')
#img.tick_labels.hide()
#img.tick_labels.hide_x()  # Hide the x axis
#img.tick_labels.hide_y()  # Hide the y axis
# Let's add a scalebar to it
img.add_scalebar(30.0/3600.) #20
img.scalebar.set_label('30"')
img.scalebar.set(color='white', linewidth=4, alpha=0.9)
img.scalebar.set_font(size=45, weight='bold',
                      stretch='normal', family='sans-serif',
                      style='normal', variant='normal')

#Filter names
img.add_label(0.1, 0.9, str(instrument_b) + "+" + str(instrument_g) + "+" + str(instrument_r), color="white",
              horizontalalignment='left',
              weight='bold', size=20, relative=True, zorder=1000)
dx, dy = 0.001, -0.001
# img.add_label(0.1+dx, 0.9+dy, instrument_b, color="black", alpha=0.6,
#               horizontalalignment='left',
#               bbox={"facecolor": "black", "edgecolor": "none",# "pad": 20,
#                     "alpha": 0.5, "boxstyle": "round, pad=0.5"},
#               weight='bold', size=55, relative=True, zorder=999)

try:
    img.show_regions(position)
except FileNotFoundError: 
    print('Region not found')
# except FileNotFoundError:
#     print("File", position,
#                 "not found - is not necesary now")

# img.show_markers(ra, dec, layer='marker_set_1', edgecolor='red',
#                  facecolor='red', marker='o', s=10, alpha=1.)
#img.show_markers(ra, dec , layer="marker_set_1", edgecolor="red", facecolor="none", marker="o", s=10,  alpha=0.5)
try:
    img.recenter(ra, dec, radius = cmd_args.zoom/3600.) #degree, for this reazon a split into 3600. #zoom ax2.recenter(ra0, dec0, 4*R0/cmd_args.zoom)
except TypeError:
    print("No zoom here")
#img.show_markers(ra, dec, layer='marker', edgecolor='red', facecolor='none', marker='o', s=10, alpha=0.9, linewidths=100.)#, layer='marker_set_1', edgecolor='black', facecolor='none', s=30, alpha=0.5, linewidths=20)
# img.scalebar.set_font(size=23, weight='bold',
#                       stretch='normal', family='sans-serif',
#                       style='normal', variant='normal')

# We may want to lengthen the scalebar, move it to the top left,
# and apply a physical scale
#img.scalebar.set_corner('top left')
# img.scalebar.set_length(20/3600.)
# img.scalebar.set_label('20 arcsec')

if cmd_args.debug:
    print("Creating of PDF image of:", position.split('-p')[0])

img.set_theme('publication')
if image_r.endswith("_swp-crop.fits"):
    img.save(image_r.replace('_swp-crop.fits', '-'+str(instrument_r)+str(instrument_g)+str(instrument_b)+'-RGB.png'))
else:
    img.save(image_r.replace('.fits', '-'+str(instrument_r)+str(instrument_g)+str(instrument_b)+'-RGB.pdf'))
