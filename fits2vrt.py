"""
This file is part of fits2vrt
A FITS to GDAL Virtual Header tool
09/01/2016
Author : Chiara Marmo (chiara.marmo@u-psud.fr)
Copyright : CNRS, Universite Paris-Sud
"""

import os
import os.path
import astropy
from astropy.io import fits
from astropy import wcs
import numpy as np
from osgeo import gdal
from osgeo import osr

"""
The FITS Image metadata
"""
class fitskeys(object):

    # Image object initialization #
    def __init__(self,imname):
        self.__name = imname
        if os.path.isfile(imname):
            hdulist = fits.open(imname)
            self.__header = hdulist[0].header
            hdulist.close()
            self.__wcs = wcs.WCS(self.__header)

    # FITS metadata conversion to GDAL VRT Header
    def fits2vrt(self):
        format = "vrt"
        driver = gdal.GetDriverByName( format )
        vrtname, fitsext = os.path.splitext(self.__name)

        # Defining dataset dimensions
        dimx = self.__header['NAXIS1']
        dimy = self.__header['NAXIS2']
        if (self.__header['NAXIS'] > 2):
            dimz = self.__header['NAXIS3']
        else:
            dimz = 1

        # Defining dataset bit type
        fbittype = self.__header['BITPIX']
        if (fbittype == 8):
            gbittype = gdal.GDT_Byte
        elif (fbittype == 16):
            gbittype = gdal.GDT_Int16
        dst_ds = driver.Create( vrtname + '.vrt', dimx, dimy, dimz, gbittype )

        # Top Left pixel is bottom left in FITS
        #topleftx = 
        #toplefty = 
        #dst_ds.SetGeoTransform( [ topleftx, wepixres, rotation, toplefty, rotation, nspixres] )
        #srs = osr.SpatialReference()
        #srs.SetWellKnownGeogCS( '' )
        #dst_ds.SetProjection( srs.ExportToWkt() )
        # we're done, close properly the dataset
        dst_ds = None

