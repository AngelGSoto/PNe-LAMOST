import matplotlib
matplotlib.use('Agg')

import aplpy
from astropy.io import fits


def isolate_image_extension(fits_file, extension):
    '''
        Saves the data + header of the specified extension as
        new FITS file

        input
        ------
        fits_file: file path to FITS image
        extension: Number of HDU extension containing the image data
    '''

    header = fits.getheader(fits_file, extension)
    data = fits.getdata(fits_file, extension)

    fits.writeto('%s_image.fits' % fits_file.rstrip('.fits'), data, header)


# Create new FITS files containing only the image extension
for exposure in ['dss_search_IR.fits', 'dss_search_red.fits', 'dss_search_blue.fits']:
    isolate_image_extension(exposure, 1)

# Call mage_rgb_cube on the isolated image extension
aplpy.make_rgb_cube(['id1k8c010_drz_image.fits', 'id1k8i010_drz_image.fits',
                     'id1k8n010_drz_image.fits'], 'id_rgb.fits')
