import matplotlib.pyplot as plt
from astropy.visualization import SqrtStretch
from astropy.visualization import ZScaleInterval
import numpy as np
from astropy.visualization import make_lupton_rgb
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename

forCasting = np.float_()

# Read in the three images downloaded from here:
g_name = 'blue.fits'
r_name = 'red.fits'
i_name = 'infrared.fits'
g = np.array(fits.open(g_name)[0].data, dtype=float)
r = np.array(fits.open(r_name)[0].data, dtype=float)
i = np.array(fits.open(i_name)[0].data, dtype=float)

rgb_default = make_lupton_rgb(i, r, g, minimum=-1000, Q=100, stretch=1, filename="PNPRTM1.pdf")
plt.imshow(rgb_default, origin='lower')
