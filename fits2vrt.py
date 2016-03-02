"""
This file is part of fits2vrt
A FITS to GDAL Virtual Header conversion tool
02/03/2016
Author : Chiara Marmo (chiara.marmo@u-psud.fr)
Copyright : CNRS, Universite Paris-Sud
"""

import os
import os.path
import re
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

            # alternate linear WCS
            ctype1 = re.compile("CTYPE1")
            px = re.compile("PX-")
            for key in self.__header.keys():
                alt = re.search(ctype1,key)
                if re.search(ctype1,key):
                    ctype = str(self.__header[key])
                    if re.search(px,ctype):
                        altkey = key[alt.end(0):]
                        self.__altkey = altkey


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
        # Following http://www.gdal.org/gdal_8h.html#a22e22ce0a55036a96f652765793fb7a4 (GDAL)
        # and https://heasarc.gsfc.nasa.gov/docs/software/fitsio/c/c_user/node20.html (FITS)
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
        # COMMENT and HISTORY keywords are excluded too
        metadata = {}
        header = self.__header
        nometadata = ['SIMPLE','EXTEND','BITPIX','BZERO','BSCALE','COMMENT','NAXIS','HISTORY']
        for i in range(1, header['NAXIS']+1):
            nometadata.append('NAXIS'+str(i))
        for key in header.keys():
            if key not in nometadata:
                metadata[key] = str(header[key])
        dst_ds.SetMetadata( metadata )

        # Defining Geotransform: if linear WCS is defined 
        # GeoTransform[1] = CD1_1
        # GeoTransform[2] = CD1_2
        # GeoTransform[4] = CD2_1
        # GeoTransform[5] = CD2_2
        # GeoTransform[0] and GeoTransform[3] must be computed.
        if self.__altkey:
            altkey = self.__altkey
            geot1 = header['CD1_1'+altkey]
            geot2 = header['CD1_2'+altkey]
            geot4 = header['CD2_1'+altkey]
            geot5 = header['CD2_2'+altkey]

            # Top Left pixel is startx endy in FITS
            #leftx, lefty = self.__wcs.wcs_pix2world(1., dimy, 1)
            topleftx = header['CRVAL1'+altkey] + geot1 * (1-header['CRPIX1'+altkey])
            toplefty = header['CRVAL2'+altkey] + geot5 * (-dimy-header['CRPIX2'+altkey])
            dst_ds.SetGeoTransform( [ topleftx, geot1, geot2, toplefty, geot4, geot5] )

        # Defining projection type
        # Following http://www.gdal.org/ogr__srs__api_8h.html (GDAL)
        # and http://www.aanda.org/component/article?access=bibcode&bibcode=&bibcode=2002A%2526A...395.1077CFUL (FITS)
        fe = 0.0
        fn = 0.0
        srs=osr.SpatialReference()
        wcsproj = (self.__header['CTYPE1'])[-3:]
        # Sinusoidal / SFL projection
        if ( wcsproj == 'SFL' ):
            clong = self.__header['CRVAL1']
            srs.SetProjection('Sinusoidal')
            srs.SetProjParm('longitude_of_center',clong)
            srs.SetProjParm('false_easting',fe)
            srs.SetProjParm('false_northing',fn)
        elif ( wcsproj == 'ZEA' ):
            dst_ds.SetProjection("Lambert_Azimuthal_Equal_Area")
        elif ( wcsproj == 'COO' ):
            dst_ds.SetProjection("Lambert_Conformal_Conic_1SP")
        elif ( wcsproj == 'CAR' ):
            dst_ds.SetProjection("Equirectangular")
        elif ( wcsproj == 'MER' ):
            dst_ds.SetProjection("Transverse_Mercator")
        elif ( wcsproj == 'SIN' ):
            dst_ds.SetProjection("Orthographic")
        elif ( wcsproj == 'AZP' ):
            dst_ds.SetProjection("perspective_point_height")
        elif ( wcsproj == 'STG' ):
            dst_ds.SetProjection("Stereographic")
        else:
            print "Unknown projection"
            print wcsproj
        wkt = srs.ExportToWkt()
        dst_ds.SetProjection(wkt)

        # Close properly the dataset
        dst_ds = None

