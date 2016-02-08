"""
This file is part of fits2vrt
A FITS to GDAL Virtual Header conversion tool
20/01/2016
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
            if (BZERO <= 0):
              gbittype = gdal.GDT_Int16
            elif (BZERO > 0):
              gbittype = gdal.GDT_UInt16
        elif (fbittype == 32):
            if (BZERO <= 0):
              gbittype = gdal.GDT_Int32
            elif (BZERO > 0):
              gbittype = gdal.GDT_UInt32
        elif (fbittype == -32):
            gbittype = gdal.GDT_Float32
        elif (fbittype == -64):
            gbittype = gdal.GDT_Float64
        else:
            print "Bit Type not supported"
            print fbittype

        dst_ds = driver.Create( vrtname + '.vrt', dimx, dimy, dimz, gbittype )

        # Setting all non mandatory FITS keywords as metadata
        metadata = {}
        header = self.__header
        nometadata = ['SIMPLE','EXTEND','BITPIX','BZERO','BSCALE','COMMENT','NAXIS']
        for i in range(1, header['NAXIS']+1):
            nometadata.append('NAXIS'+str(i))
        for key in header.keys():
            if key not in nometadata:
                metadata[key] = str(header[key])
        dst_ds.SetMetadata( metadata )

        # Defining projection type
        wcsproj = (self.__header['CTYPE1'])[-3:]
        if ( wcsproj == 'SFL' ):
            print wcsproj
            #dst_ds.SetSinusoidal (double dfCenterLong, double dfFalseEasting, double dfFalseNorthing)
        elif ( wcsproj == 'ZEA' ):
            print wcsproj
            #dst_ds.SetLAEA (double dfCenterLat, double dfCenterLong, double dfFalseEasting, double dfFalseNorthing)
        elif ( wcsproj == 'COO' ):
            print wcsproj
            #dst_ds.SetLCC (double dfStdP1, double dfStdP2, double dfCenterLat, double dfCenterLong, double dfFalseEasting, double dfFalseNorthing)
        elif ( wcsproj == 'CAR' ):
            print wcsproj
            #dst_ds.SetEquirectangular (double dfCenterLat, double dfCenterLong, double dfFalseEasting, double dfFalseNorthing)
        elif ( wcsproj == 'MER' ):
            print wcsproj
            #dst_ds.SetTM (double dfCenterLat, double dfCenterLong, double dfScale, double dfFalseEasting, double dfFalseNorthing)
        elif ( wcsproj == 'SIN' ):
            print wcsproj
            #dst_ds.SetOrthographic (double dfCenterLat, double dfCenterLong, double dfFalseEasting, double dfFalseNorthing)
        elif ( wcsproj == 'AZP' ):
            print wcsproj
            #dst_ds.zenithalperspective point perspective?
        elif ( wcsproj == 'STG' ):
            print wcsproj
            #dst_ds.SetStereographic (double dfCenterLat, double dfCenterLong, double dfScale, double dfFalseEasting, double dfFalseNorthing)
        else:
            print "Unknown projection"
            print wcsproj

        # Top Left pixel is bottom left in FITS
        #topleftx = 
        #toplefty = 
        #dst_ds.SetGeoTransform( [ topleftx, wepixres, rotation, toplefty, rotation, nspixres] )
        #srs = osr.SpatialReference()
        #srs.SetWellKnownGeogCS( '' )
        #dst_ds.SetProjection( srs.ExportToWkt() )

        # Close properly the dataset
        dst_ds = None

