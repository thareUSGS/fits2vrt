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
        #dst_ds.SetProjParm(false_northing,0.0)
        #dst_ds.SetProjParm(false_easting,0.0)
        wcsproj = (self.__header['CTYPE1'])[-3:]
        if ( wcsproj == 'SFL' ):
            #dfCenterLong = self.__header['CRVAL1']
            dst_ds.SetProjection("Sinusoidal")
            #dst_ds.SetProjParm(longitude_of_center,dfCenterLong)
        elif ( wcsproj == 'ZEA' ):
            dst_ds.SetProjection("Lambert_Azimuthal_Equal_Area")
        elif ( wcsproj == 'COO' ):
            dst.SetProjection("Lambert_Conformal_Conic_1SP")
        elif ( wcsproj == 'CAR' ):
            dst.SetProjection("Equirectangular")
        elif ( wcsproj == 'MER' ):
            dst.SetProjection("Transverse_Mercator")
        elif ( wcsproj == 'SIN' ):
            dst.SetProjection("Orthographic")
        elif ( wcsproj == 'AZP' ):
            dst.SetProjection("perspective_point_height")
        elif ( wcsproj == 'STG' ):
            dst.SetProjection("Stereographic")
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

