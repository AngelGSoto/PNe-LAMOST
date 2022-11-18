import matplotlib.pyplot as plt
from astropy.visualization import SqrtStretch
from astropy.visualization import ZScaleInterval
import numpy as np
from astropy.visualization import make_lupton_rgb
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename

forCasting = np.float_()

# Read in the three images downloaded from here:
# g_name = 'blue.fits'
# r_name = 'red.fits'
# i_name = 'infrared.fits'
# g = np.array(fits.open(g_name)[0].data, dtype=float)
# r = np.array(fits.open(r_name)[0].data, dtype=float)
# i = np.array(fits.open(i_name)[0].data, dtype=float)

# rgb_default = make_lupton_rgb(i, r, g, minimum=-1000, Q=100, stretch=1, filename="PNPRTM1.pdf")
# plt.imshow(rgb_default, origin='lower')


def RGB_lupton(x, rgb_q=15, rgb_stretch=0.5, rgb_min=0):
    if x.ndim==3:
        x = make_lupton_rgb(x[:,:,2], x[:,:,1], x[:,:,0],
                      Q=rgb_q, stretch=rgb_stretch, minimum=rgb_min)
    elif x.ndim==4:
        x = np.array([make_lupton_rgb(xi[:,:,2], xi[:,:,1], xi[:,:,0],
                      Q=rgb_q, stretch=rgb_stretch, minimum=rgb_min)
                      for xi in x])
    else:
        raise ValueError(f"Wrong number of dimensions! Gave {x.ndim}, need 3 or 4")
    return x


# α=146.080065° δ=-0.653004° ∠=0.101864
# 146.080065, -0.653004
fits_dir = '../data/fits_files'
fits_fns = ['dss_search_blue.fits',
            'dss_search_red.fits',
            'dss_search_IR.fits'
            ]

# fits_fns = ['cutout-HSC-G-9325-pdr2_wide-201204-201802.fits',
#             'cutout-HSC-R-9325-pdr2_wide-201204-201806.fits',
#             'cutout-HSC-I-9325-pdr2_wide-201204-153311.fits'
#             ]

data = []
for fits_fn in fits_fns:
  hdul = fits.open(f'{fits_fn}')
  data.append(hdul[1].data)

data = np.array(data)

print(data.shape)
#data = data.transpose(1, 2, 0)
#print(data.shape)

datalup = RGB_lupton(data)
