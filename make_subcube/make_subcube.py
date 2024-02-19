import getpass
import itertools
import os
import socket
import warnings

import numpy as np

from astropy import units as u
from astropy.convolution import Gaussian1DKernel, Gaussian2DKernel, convolve
from astropy.io import fits
from astropy.wcs import WCS
from datetime import datetime
from tqdm import tqdm



def make_subcube(slice_params, path_to_file=None, hdu=None, dtype='float32',
                 save=False, overwrite=True, path_to_output_file=None,
                 get_hdu=False, get_data=True, get_header=True):
    """Extract subcube from a spectral cube.

    Parameters
    ----------
    slice_params : list
        List of slice parameters for cube, if only a subset of the data should be used.
    path_to_file : str
        Filepath to the FITS cube.
    hdu : astropy.io.fits.HDUList
        Header/Data unit of the FITS cube.
    dtype : str
        Data type to which the array should be transformed. Default is `float32`.
    save : bool
        Default is `False`. If set to `True`, the resulting FITS cube is saved under 'path_to_output_file'.
    overwrite : bool
        If set to `True`, overwrites any already existing files saved in `path_to_output_file`.
    path_to_output_file : type
        Filepath to which subcube gets saved.
    get_hdu : bool
        Default is `False`. If set to `True`, an astropy.io.fits.HDUList is returned. Overrides 'get_data' and 'get_header'.
    get_data : bool
        Default is `True`. Returns a numpy.ndarray of the FITS array.
    get_header : bool
        Default is `True`. Returns a astropy.io.fits.Header of the FITS array.

    """
    def save_fits(data, header, path_to_file, verbose=True):
        """Save data array and header as FITS file.

        Parameters
        ----------
        data : numpy.ndarray
            Data array.
        header : astropy.io.fits.Header
            Header of the FITS array.
        path_to_file : str
            Filepath to which FITS array should get saved.
        verbose : bool
            Default is `True`. Writes message to terminal about where the FITS file was saved.

        """
        if not os.path.exists(os.path.dirname(path_to_file)):
            os.makedirs(os.path.dirname(path_to_file))
        fits.writeto(path_to_file, data, header=header, overwrite=True)

    print('\nmaking subcube with the slice parameters {}...'.format(
        slice_params))

    # check_if_all_values_are_none(hdu, path_to_file, 'hdu', 'path_to_file')

    if path_to_file is not None:
        hdu = fits.open(path_to_file)[0]

    data = hdu.data
    header = hdu.header

    data = data[slice_params[0], slice_params[1], slice_params[2]]
    data = data.astype(dtype)
    wcs = WCS((header))
    wcs_cropped = wcs[slice_params[0], slice_params[1], slice_params[2]]
    header.update(wcs_cropped.to_header())

    header['NAXIS1'] = data.shape[2]
    header['NAXIS2'] = data.shape[1]
    header['NAXIS3'] = data.shape[0]
    header['COMMENT'] = "Cropped FITS file ({})".format(
            datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%S'))

 

    save_fits(data, header, path_to_output_file, verbose=True)

    # return return_hdu_options(
    #     fits.PrimaryHDU(data, header), get_hdu=get_hdu, get_data=get_data, get_header=get_header)



    # if verbose:
    #     save_file(os.path.basename(path_to_file), os.path.dirname(path_to_file))