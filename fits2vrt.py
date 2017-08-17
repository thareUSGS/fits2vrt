"""
This file is part of fits2vrt
A FITS to GDAL Virtual Header conversion tool
17/08/2017
Authors : Chiara Marmo (chiara.marmo@u-psud.fr)
          Trent Hare (thare@usgs.gov)
Copyright : CNRS, Universite Paris-Sud - USGS
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
            #not sure what this is (? from Trent)
            #irs.SetProjection(gdalproj)
            lenheader = (len(hdulist[0].header)+1)*80
            [blocks,remainder] = divmod(lenheader,2880)
            self.__offset = 2880 * (blocks+1)
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
        #src_ds = gdal.Open(self.__name)
        driver = gdal.GetDriverByName('vrt')
        vrtname, fitsext = os.path.splitext(self.__name)

        # Defining dataset dimensions
        dimx = self.__header['NAXIS1']
        dimy = self.__header['NAXIS2']
        if (self.__header['NAXIS'] > 2):
            dimz = self.__header['NAXIS3']
        else:
            dimz = 1

        # Defining Scale and Offset
        #bzero = self.__header['BZERO']
        #bscale = self.__header['BSCALE']
        bzero = 0.0
        bscale = 1.0

        # Defining dataset bit type
        # Following http://www.gdal.org/gdal_8h.html#a22e22ce0a55036a96f652765793fb7a4 (GDAL)
        # and https://heasarc.gsfc.nasa.gov/docs/software/fitsio/c/c_user/node20.html (FITS)
        fbittype = self.__header['BITPIX']
        if (fbittype == 8):
            gbittype = gdal.GDT_Byte
            nodata = 0
            #pxoffset needed for raw VRT type (points into FITS image)
            pxoffset = 1
        elif (fbittype == 16):
            pxoffset = 2
            if (bzero <= 0):
              gbittype = gdal.GDT_Int16
              nodata = -32768
            elif (bzero > 0):
              gbittype = gdal.GDT_UInt16
              nodata = 0
        elif (fbittype == 32):
            pxoffset = 4
            if (bzero <= 0):
              gbittype = gdal.GDT_Int32
            elif (bzero > 0):
              gbittype = gdal.GDT_UInt32
        elif (fbittype == -32):
            pxoffset = 4
            gbittype = gdal.GDT_Float32
            nodata = -3.40282e+38
        elif (fbittype == -64):
            pxoffset = 8
            gbittype = gdal.GDT_Float64
        else:
            print "Bit Type not supported"
            print fbittype
            return 'fail' # need better error text or method

        # Addressing FITS as raw raster: this will work without CFITSIO GDAL dependence.
        dst_ds = driver.Create( vrtname + '.vrt', dimx, dimy, 0 )

        # The next lines only work if gdal has cfitsio properly configured
        #src_ds = gdal.Open(self.__name)
        #dst_ds = driver.CreateCopy( vrtname + '.vrt', src_ds, 0)

        #lineoffset only needed for raw VRT type
        lnoffset = dimx * pxoffset
        src_filename_opt = 'SourceFileName=' + self.__name
        im_offset_opt = 'ImageOffset=' + str(self.__offset)
        px_offset_opt = 'PixelOffset=' + str(pxoffset)
        ln_offset_opt = 'LineOffset=' + str(lnoffset)
        options = [
           'subClass=VRTRawRasterBand',
           src_filename_opt,
           'relativeToVRT=1',
           im_offset_opt,
           px_offset_opt,
           ln_offset_opt,
           'ByteOrder=MSB' #FITS is always MSB
        ]

        result = dst_ds.AddBand( gbittype, options )
        if result != gdal.CE_None:
            print 'AddBand() returned error code'
            return 'fail'

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

        # Defining Target object
        try:
            target = header['OBJECT'] 
        except:
            print "OBJECT keyword is missing"
            target = 'Undefined'


        # Defining Geotransform: if linear WCS is defined 
        # GeoTransform[1] = CD1_1
        # GeoTransform[2] = CD1_2
        # GeoTransform[4] = CD2_1
        # GeoTransform[5] = CD2_2
        # GeoTransform[0] and GeoTransform[3] must be computed.
        try:
            altkey = self.__altkey
            geot1 = header['CD1_1'+altkey]
            geot2 = header['CD1_2'+altkey]
            geot4 = header['CD2_1'+altkey]
            geot5 = header['CD2_2'+altkey]

            # FITS rasters are still read upside-down by GIS software UpperLeftCorner is LowerLeftCorner for now
            topleftx = header['CRVAL1'+altkey] + geot1 * ( - header['CRPIX1'+altkey]) + geot2 * ( - header['CRPIX2'+altkey])
            toplefty = header['CRVAL2'+altkey] + geot5 * (1 - header['CRPIX2'+altkey]) + geot4 * (1 - header['CRPIX1'+altkey])
            dst_ds.SetGeoTransform( [ topleftx, geot1, geot2, toplefty, geot4, geot5] )
        except:
            print "WARNING! No linear keyword available, geotransformation matrix will not be calculated."

        # Defining projection type
        # Following http://www.gdal.org/ogr__srs__api_8h.html (GDAL)
        # and http://www.aanda.org/component/article?access=bibcode&bibcode=&bibcode=2002A%2526A...395.1077CFUL (FITS)

        srs = osr.SpatialReference()
        try:
            # Get radius values (maybe add an external method (e.g. string or URI)
            # new FITS keywords A_RADIUS, C_RADIUS 
            # note B_RADIUS (to define a triaxial) not generally used for mapping applications
            semiMajor = header['A_RADIUS']
            semiMinor = header['C_RADIUS']
            if ((semiMajor - semiMinor) > 0.00001):
                invFlattening= semiMajor / ( semiMajor - semiMinor)
            else:
                invFlattening= 0.0

            gcsName   = 'GCS_' + target
            datumName = 'D_' + target
            srsGeoCS  = 'GEOGCS["'+gcsName+'",DATUM["'+datumName+'",SPHEROID["'+ \
                        target+'",' + str(semiMajor) + ',' + str(invFlattening) + \
                        ']], PRIMEM["Reference_Meridian",0],' + \
                        'UNIT["degree",0.0174532925199433]]'
                   
            #print srsGeoCS
            srs.ImportFromWkt(srsGeoCS)
        except:
            print "WARNING! No Radii keyword available, metadata will not contain DATUM information."
            

        wcsproj = (self.__header['CTYPE1'])[-3:]
        # Sinusoidal / SFL projection
        if ( wcsproj == 'SFL' ):
            gdalproj = 'Sinusoidal'
            srs.SetProjection(gdalproj)
            clong = self.__header['CRVAL1']
            srs.SetProjParm('longitude_of_center',clong)

        # Lambert Azimuthal Equal Area / ZEA projection
        elif ( wcsproj == 'ZEA' ):
            gdalproj = 'Lambert_Azimuthal_Equal_Area'
            srs.SetProjection(gdalproj)
            clong = self.__header['CRVAL1']
            srs.SetProjParm('longitude_of_center',clong)
            #clat = self.__header['XXXXX']
            #srs.SetProjParm('latitude_of_center',clat)

        # Lambert Conformal Conic 1SP / COO projection
        elif ( wcsproj == 'COO' ):
            gdalproj = 'Lambert_Conformal_Conic_1SP'
            srs.SetProjection(gdalproj)
            clong = self.__header['CRVAL1']
            srs.SetProjParm('longitude_of_center',clong)
            #clat = self.__header['XXXXX']
            #srs.SetProjParm('latitude_of_center',clat)
            scale = self.__header['XXXXX']
            if scale is not None:
                srs.SetProjParm('scale_factor',scale)
            else: #set default of 1.0
                srs.SetProjParm('scale_factor',1.0)

        # Equirectangular / CAR projection
        elif ( wcsproj == 'CAR' ):
            gdalproj = 'Equirectangular'
            srs.SetProjection(gdalproj)
            cmer = self.__header['CRVAL1']
            srs.SetProjParm('central_meridian',cmer)
            # The standard_parallel_1 defines where the local radius is calculated
            # not the center of Y Cartesian system (which is latitude_of_origin)
            # But FITS WCS only supports projections on the sphere
            # we assume here that the local radius is the one computed at the projection center
            spar = self.__header['CRVAL2']
            srs.SetProjParm('standard_parallel_1',spar)
            olat = self.__header['CRVAL2']
            srs.SetProjParm('latitude_of_origin',olat)

	#Here we are using Mercator not Transverse Mercator but
	#There is a change FITS might be Hotine Merc or Trans Merc
        # Mercator / MER projection
        elif ( wcsproj == 'MER' ):  
            gdalproj = 'Mercator'
            srs.SetProjection(gdalproj)
            cmer = self.__header['CRVAL1']
            srs.SetProjParm('central_meridian',cmer)
            #olat = self.__header['XXXXX']
            #srs.SetProjParm('latitude_of_origin',olat)
            if olat is not None:
                srs.SetProjParm('scale_factor',olat)
            else: #set default of 0.0
                srs.SetProjParm('scale_factor',0.0)

        # Orthographic / SIN projection
        elif ( wcsproj == 'SIN' ):
            gdalproj = 'Orthographic'
            srs.SetProjection(gdalproj)
            cmer = self.__header['CRVAL1']
            srs.SetProjParm('central_meridian',cmer)
            #olat = self.__header['XXXXX']
            #srs.SetProjParm('latitude_of_origin',olat)

        # Point Perspective / AZP projection
        elif ( wcsproj == 'AZP' ):
            gdalproj = 'perspective_point_height'
            srs.SetProjection(gdalproj)
            # appears to need height... maybe center lon/lat

        # Polar Stereographic / STG projection
        elif ( wcsproj == 'STG' ):
            gdalproj = 'Polar_Stereographic'
            srs.SetProjection(gdalproj)
            cmer = self.__header['CRVAL1']
            srs.SetProjParm('central_meridian',cmer)
            olat = self.__header['CRVAL2']
            srs.SetProjParm('latitude_of_origin',olat)

        else:
            print "Unknown projection"
            return 'fail'

        projname = gdalproj + '_' + target
        srs.SetAttrValue('projcs',projname)
        # false easting and northing not used for planetary
        srs.SetProjParm('false_easting',0.0)
        srs.SetProjParm('false_northing',0.0)

        wkt = srs.ExportToWkt()
        print 'projection:\n'+wkt
        dst_ds.SetProjection(wkt)

        # Adding SimpleSource info
        band = dst_ds.GetRasterBand(1)
        band.SetNoDataValue(nodata)
        band.SetScale(bscale)
        band.SetOffset(bzero)
        
        # Close properly the dataset
        dst_ds = None
        #src_ds = None

