#!/usr/bin/env python
#/******************************************************************************
# * $Id$
# *
# * Project: GDAL Utilities
# * Purpose: Create a "geo"FITS compatible image from a GDAL supported image.
# * Authors:  Trent Hare, U.S. Geological Survey, <thare@usgs.gov>
# *           Chiara Marmo Paris-Sud, France
# * Date:    June 27, 2017
# * version: 0.1
# *
# * Port from gdalinfo.py whose author is Even Rouault
# ******************************************************************************
# * Copyright (c) 2010, Even Rouault
# * Copyright (c) 1998, Frank Warmerdam
# *
# * Permission is hereby granted, free of charge, to any person obtaining a
# * copy of this software and associated documentation files (the "Software"),
# * to deal in the Software without restriction, including without limitation
# * the rights to use, copy, modify, merge, publish, distribute, sublicense,
# * and/or sell copies of the Software, and to permit persons to whom the
# * Software is furnished to do so, subject to the following conditions:
# *
# * The above copyright notice and this permission notice shall be included
# * in all copies or substantial portions of the Software.
# *
# * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# * DEALINGS IN THE SOFTWARE.
# ****************************************************************************/

import sys
import math
import os
import os.path
import re
import numpy as np

import astropy
from astropy.io import fits
from astropy import wcs

try:
    from osgeo import gdal
    from osgeo import osr
except:
    import gdal
    import osr

#/************************************************************************/
#/*                               Usage()                                */
#/************************************************************************/

def Usage(theApp):
    print( '\nUsage: gdal2fits.py in.tif output.fits') # % theApp)
    print( '   optional: to print out image information also send -debug')
    print( '   optional: to get lonsys=360, send -force360')
    print( '   optional: to computer min and maximum for FITS -computeMinMax')
    print( '   optional: to override the center Longitude, send -centerLon 180')
    print( '   optional: to set or override scaler and offset send -base 17374000 and/or -multiplier 0.5')
    sys.exit(1)


def EQUAL(a, b):
    return a.lower() == b.lower()


#/************************************************************************/
#/*                                main()                                */
#/************************************************************************/

def main( argv = None ):

    debug = False
    inFilename = None
    inProjection = None
    bShowFileList = True
    bComputeMinMax = False
    dst_fits = None
    bands = 1
    centLat = 0
    centLon = 0
    centerLon = False
    TMscale = 1.0
    UpperLeftCornerX = 0.0
    UpperLeftCornerY = 0.0
    falseEast = 0
    falseNorth = 0
    force360 = False
    base = None
    multiplier = None

    if argv is None:
        argv = sys.argv

    argv = gdal.GeneralCmdLineProcessor( argv )

    if argv is None:
        return 1

    nArgc = len(argv)

#/* -------------------------------------------------------------------- */
#/*      Parse arguments.                                                */
#/* -------------------------------------------------------------------- */
    i = 1
    while i < nArgc:

        if EQUAL(argv[i], "--utility_version"):
            print("%s is running against GDAL %s" % \
                   (argv[0], gdal.VersionInfo("RELEASE_NAME")))
            return 0
        elif EQUAL(argv[i], "-debug"):
            debug = True
        elif EQUAL(argv[i], "-computeMinMax"):
            bComputeMinMax = True
        elif EQUAL(argv[i], "-force360"):
            force360 = True
        elif EQUAL(argv[i], "-centerLon"):
            i = i + 1
            centerLon = float(argv[i])
        elif EQUAL(argv[i], "-base"):
            i = i + 1
            base = float(argv[i])
        elif EQUAL(argv[i], "-multiplier"):
            i = i + 1
            multiplier = float(argv[i])
        elif argv[i][0] == '-':
            return Usage(argv[0])
        elif inFilename is None:
            inFilename = argv[i]
        elif dst_fits is None:
            dst_fits = argv[i]
        else:
            return Usage(argv[0])

        i = i + 1

    if inFilename is None:
        return Usage(argv[0])
    if dst_fits is None:
        return Usage(argv[0])

#/* -------------------------------------------------------------------- */
#/*      Open dataset.                                                   */
#/* -------------------------------------------------------------------- */
    inDataset = gdal.Open( inFilename, gdal.GA_ReadOnly )
    if inDataset is None:
        print("gdalinfo failed - unable to open '%s'." % inFilename )
        sys.exit(1)

    # Open the output file.

#/* -------------------------------------------------------------------- */
#/*      Report general info.                                            */
#/* -------------------------------------------------------------------- */
    hDriver = inDataset.GetDriver();
    if debug:
        print( "Driver: %s/%s" % ( \
                hDriver.ShortName, \
                hDriver.LongName ))

    papszFileList = inDataset.GetFileList();
    if papszFileList is None or len(papszFileList) == 0:
        print( "Files: none associated" )
    else:
        if debug:
            print( "Files: %s" % papszFileList[0] )
            if bShowFileList:
                for i in range(1, len(papszFileList)):
                    print( "       %s" % papszFileList[i] )

    if debug:
        print( "Size is %d, %d" % (inDataset.RasterXSize, inDataset.RasterYSize))


#/* -------------------------------------------------------------------- */
#/*      Report projection.                                              */
#/* -------------------------------------------------------------------- */
    inProjection = inDataset.GetProjectionRef()
    if inProjection is not None:

        hSRS = osr.SpatialReference()
        if hSRS.ImportFromWkt(inProjection ) == gdal.CE_None:
            pszPrettyWkt = hSRS.ExportToPrettyWkt(False)
            if debug:
                print( "Coordinate System is:\n%s" % pszPrettyWkt )
        else:
            if debug:
                print( "Coordinate System is `%s'" % inProjection )


        hSRS = osr.SpatialReference()
        if hSRS.ImportFromWkt(inProjection) == gdal.CE_None:
            pszPrettyWkt = hSRS.ExportToPrettyWkt(False)

            #print( "Coordinate System is:\n%s" % pszPrettyWkt )
            mapProjection = "None"
            #Extract projection information
            target = hSRS.GetAttrValue("DATUM",0)
            target = target.replace("D_","").replace("_2000","").replace("GCS_","")

            semiMajor = hSRS.GetSemiMajor() 
            semiMinor = hSRS.GetSemiMinor()
            # if image is in degrees (deg/pix) then force meters (can be removed)
            if (inProjection[0:6] == "GEOGCS"):
                mapProjection = "CAR"
                centLon = hSRS.GetProjParm('central_meridian')

            if (inProjection[0:6] == "PROJCS"):
                mapProjection = hSRS.GetAttrValue("PROJECTION",0)

                if EQUAL(mapProjection,"Sinusoidal"):
                    mapProjection = "SFL"
                    centLon = hSRS.GetProjParm('central_meridian')

                if EQUAL(mapProjection,"Equirectangular"):
                    mapProjection = "CAR"
                    centLat = hSRS.GetProjParm('standard_parallel_1')
                    centLon = hSRS.GetProjParm('central_meridian')

                #Transverse Mercator maybe not supported in FITS....?
                #if EQUAL(mapProjection,"Transverse_Mercator"):
                #    mapProjection = "MER"
                #    centLat = hSRS.GetProjParm('standard_parallel_1')
                #    centLon = hSRS.GetProjParm('central_meridian')
                #    TMscale = hSRS.GetProjParm('scale_factor')
                #    #Need to research when TM actually applies false values
                #    #but planetary is almost always 0.0
                #    falseEast =  hSRS.GetProjParm('false_easting')
                #    falseNorth =  hSRS.GetProjParm('false_northing')

                if EQUAL(mapProjection,"Orthographic"):
                    mapProjection = "SIN"
                    centLat = hSRS.GetProjParm('standard_parallel_1')
                    centLon = hSRS.GetProjParm('central_meridian')

                #Mercator NOT supported in FITS....?
                if (EQUAL(mapProjection,"Mercator_1SP") or EQUAL(mapProjection,"Mercator")):
                    mapProjection = "MER" # a guess
                    centLat = hSRS.GetProjParm('standard_parallel_1')
                    centLon = hSRS.GetProjParm('central_meridian')

                if EQUAL(mapProjection,"Polar_Stereographic"):
                    mapProjection = "STG"
                    centLat = hSRS.GetProjParm('latitude_of_origin')
                    centLon = hSRS.GetProjParm('central_meridian')

                if EQUAL(mapProjection,"Stereographic_South_Pole"):
                    mapProjection = "STG"
                    centLat = hSRS.GetProjParm('latitude_of_origin')
                    centLon = hSRS.GetProjParm('central_meridian')
                
                if EQUAL(mapProjection,"Stereographic_North_Pole"):
                    mapProjection = "STG"
                    centLat = hSRS.GetProjParm('latitude_of_origin')
                    centLon = hSRS.GetProjParm('central_meridian')
        else:
            print( "Warning - Currently we can't parse this type of projection" )
            print( "Coordinate System is `%s'" % inProjection )
            target = "n/a"
            #sys.exit(1)
    else:
        print( "Warning - No Coordinate System defined:\n" )
        target = "n/a"
        #sys.exit(1)
        
#/* -------------------------------------------------------------------- */
#/*      Report Geotransform.                                            */
#/* -------------------------------------------------------------------- */
    adfGeoTransform = inDataset.GetGeoTransform(can_return_null = True)
    if adfGeoTransform is not None:
        #figure it out for fits
        UpperLeftCornerX = adfGeoTransform[0] - falseEast
        UpperLeftCornerY = adfGeoTransform[3] - falseNorth

        if adfGeoTransform[2] == 0.0 and adfGeoTransform[4] == 0.0:
            if debug:
                print( "Origin = (%.15f,%.15f)" % ( \
                        adfGeoTransform[0], adfGeoTransform[3] ))

                print( "Pixel Size = (%.15f,%.15f)" % ( \
                        adfGeoTransform[1], adfGeoTransform[5] ))

        else:
            if debug:
                print( "GeoTransform =\n" \
                        "  %.16g, %.16g, %.16g\n" \
                        "  %.16g, %.16g, %.16g" % ( \
                        adfGeoTransform[0], \
                        adfGeoTransform[1], \
                        adfGeoTransform[2], \
                        adfGeoTransform[3], \
                        adfGeoTransform[4], \
                        adfGeoTransform[5] ))

        #Using a very simple method to calculate cellsize from degree to meters. 
        #Warning: might not always be good.
        if (inProjection[0:6] == "GEOGCS"):
            #convert degrees/pixel to m/pixel 
             mapres = 1 / adfGeoTransform[1]
             mres = adfGeoTransform[1] * (semiMajor * math.pi / 180.0)
        else:
            #convert m/pixel to pixel/degree
             mapres = 1 / (adfGeoTransform[1] / (semiMajor * math.pi / 180.0))
             mres = adfGeoTransform[1] 

        #from fits2vrt notes
        # Defining Geotransform: if linear WCS is defined 
        # GeoTransform[1] = CD1_1
        # GeoTransform[2] = CD1_2
        # GeoTransform[4] = CD2_1
        # GeoTransform[5] = CD2_2
        # GeoTransform[0] and GeoTransform[3] must be computed.

             
#/* -------------------------------------------------------------------- */
#/*      Setup projected to lat/long transform if appropriate.           */
#/* -------------------------------------------------------------------- */
    if inProjection is not None and len(inProjection) > 0:
        hProj = osr.SpatialReference( inProjection )
        if hProj is not None:
            hLatLong = hProj.CloneGeogCS()

        if hLatLong is not None:
            gdal.PushErrorHandler( 'CPLQuietErrorHandler' )
            hTransform = osr.CoordinateTransformation( hProj, hLatLong )
            gdal.PopErrorHandler()
            if gdal.GetLastErrorMsg().find( 'Unable to load PROJ.4 library' ) != -1:
                hTransform = None

#/* -------------------------------------------------------------------- */
#/*      Report corners.                                                 */
#/* -------------------------------------------------------------------- */
    if debug:
        print( "Corner Coordinates:" )
        GDALInfoReportCorner( inDataset, hTransform, "Upper Left", \
                              0.0, 0.0 );
        GDALInfoReportCorner( inDataset, hTransform, "Lower Left", \
                              0.0, inDataset.RasterYSize);
        GDALInfoReportCorner( inDataset, hTransform, "Upper Right", \
                              inDataset.RasterXSize, 0.0 );
        GDALInfoReportCorner( inDataset, hTransform, "Lower Right", \
                              inDataset.RasterXSize, \
                              inDataset.RasterYSize );
        GDALInfoReportCorner( inDataset, hTransform, "Center", \
                              inDataset.RasterXSize/2.0, \
                              inDataset.RasterYSize/2.0 );

    #Get bounds -- Do not need
    #ulx = GDALGetLon( inDataset, hTransform, 0.0, 0.0 );
    #uly = GDALGetLat( inDataset, hTransform, 0.0, 0.0 );
    #lrx = GDALGetLon( inDataset, hTransform, inDataset.RasterXSize, \
    #                      inDataset.RasterYSize );
    #lry = GDALGetLat( inDataset, hTransform, inDataset.RasterXSize, \
    #                      inDataset.RasterYSize );
   
    #If user updates the centerLon, set here
    if (centerLon):
       centLon = centerLon
 
    #Calculate Simple Cylindrical X,Y in meters from bounds if not projected.
    #Needs more testing.                     
    if (inProjection[0:6] == "GEOGCS"):
        #note that: mres = adfGeoTransform[1] * (semiMajor * math.pi / 180.0)
        UpperLeftCornerX = semiMajor * (ulx - centLon) * math.pi / 180.0
        UpperLeftCornerY = semiMajor * uly * math.pi / 180.0
        

    #Initialize output image and header
    tofits = open(dst_fits,'wt')


#/* ==================================================================== */
#/*      Loop over bands to write image.                                                */
#/* ==================================================================== */
    bands = inDataset.RasterCount
    for i in range(1, inDataset.RasterCount + 1):

        iBand = inDataset.GetRasterBand(i)
        (nBlockXSize, nBlockYSize) = iBand.GetBlockSize()
        if debug:
                print( "Band %d Block=%dx%d Type=%s, ColorInterp=%s" % ( i, \
                       nBlockXSize, nBlockYSize, \
                       gdal.GetDataTypeName(iBand.DataType), \
                       gdal.GetColorInterpretationName( \
                       iBand.GetRasterColorInterpretation()) ))

                if iBand.GetDescription() is not None \
                                          and len(iBand.GetDescription()) > 0 :
                    print( "  Description = %s" % iBand.GetDescription() )

        dfMin = iBand.GetMinimum()
        dfMax = iBand.GetMaximum()
        if dfMin is not None or dfMax is not None or bComputeMinMax:
                    line =  "  "
                    if dfMin is not None:
                        line = line + ("Min=%.3f " % dfMin)
                    if dfMax is not None:
                        line = line + ("Max=%.3f " % dfMax)

                    if bComputeMinMax:
                        gdal.ErrorReset()
                        adfCMinMax = iBand.ComputeRasterMinMax(False)
                        dfMin = adfCMinMax[0]
                        dfMax = adfCMinMax[1]
                        if gdal.GetLastErrorType() == gdal.CE_None:
                            line = line + ( "  Computed Min/Max=%.3f,%.3f" % ( \
                                      dfMin, dfMax ))
                    if debug:
                        print( line )


        dfNoData = iBand.GetNoDataValue()
        if dfNoData is not None:
            if debug:
                if dfNoData != dfNoData:
                    print( "  NoData Value=nan" )
                else:
                    print( "  NoData Value=%.18g" % dfNoData )


        if debug:               
            if iBand.GetScale() != 1.0 or iBand.GetOffset() != 0.0:
                print( "  Offset: %.15g,   Scale:%.15g" % \
                         ( iBand.GetOffset(), iBand.GetScale()))

        #Load band into numpy array for writing using Astropy
        raster_data = iBand.ReadAsArray(0, 0, iBand.XSize, iBand.YSize)

        #Note we are currently loading full band. This can be changed to per line
        #for i in range(iBand.YSize):
        #   scanLine = iBand.ReadAsArray(0, i, iBand.XSize, 1, iBand.XSize, 1)

   
        #get the datatype
        if EQUAL(gdal.GetDataTypeName(iBand.DataType), "Float32"):
            fbittype = -32
        elif EQUAL(gdal.GetDataTypeName(iBand.DataType), "Float64"):
            fbittype = -64
        elif EQUAL(gdal.GetDataTypeName(iBand.DataType), "INT64"):
            fbittype = 64
        elif EQUAL(gdal.GetDataTypeName(iBand.DataType), "INT32"):
            fbittype = 32
        elif EQUAL(gdal.GetDataTypeName(iBand.DataType), "INT16"):
            fbittype = 16
        elif EQUAL(gdal.GetDataTypeName(iBand.DataType), "UINT16"):
            fbittype = 16
            bzero = -32768
            bscale = 1
        elif EQUAL(gdal.GetDataTypeName(iBand.DataType), "Byte"):
            fbittype = 8
        else:
            print( "  %s: Not supported pixel type. Please convert to 8, 16 Int, or 32 Float" % gdal.GetDataTypeName(iBand.DataType))
            sys.exit(1)

        if debug:
            print "GDAL type: %s" % gdal.GetDataTypeName(iBand.DataType)
            print "FITS type: %s" % str(fbittype)

        #Create FITS using this information
        # filename : dst_fits
        # NAXIS1   : inDataset.RasterXSize
        # NAXIS2   : inDataset.RasterYSize
        # NAXIS3   : inDataset.RasterCount
        # BITPIX   : fbittype
        # ByteOrder  = MSB  - but need to write out this correctly
        # BZERO section - allows for user to set
        
        #if base is None: 
        #    BZERO   : iBand.GetOffset()
        #else:
        #    BZERO   : base
        
        # BSCALE section - allows for user to set
        #if multiplier is None: 
        #    BSCALE  : iBand.GetScale()
        #else:
        #    BSCALE  : multiplier
        
        # OBJECT   : target
        # CTYPE1   : mapProjection
        # A_RADIUS : semiMajor
        # B_RADIUS : SemiMajor
        # C_RADIUS : SemiMinor
        
        #Need to inverse these to get values (from fits2vrt.py)
        #topleftx = header['CRVAL1'+altkey] + geot1 * (0.5 - header['CRPIX1'+altkey]) + geot2 * (0.5 - header['CRPIX2'+altkey])
        #toplefty = header['CRVAL2'+altkey] + geot5 * (0.5 - header['CRPIX2'+altkey]) + geot4 * (0.5 - header['CRPIX1'+altkey])

        # CRPIX1   : UpperLeftCornerX (BUT calc to pixel space)
        # CRPIX2   : UpperLeftCornerY (BUT calc to pixel space)
        
        
        if ((centLon < 0) and force360):
           centLon = centLon + 360
        # CRVAL1   : centLon  # not sure this is correct
        # CRVAL2   : centLat  # not sure this is correct
        # CRPIX1   : need to calc
        # CRPIX2   : need to calc
        
        #####################################################
        # Below this line, none of these are probably needed
        # pixel size in degrees : mapres 
        # MinimumLatitude : lry
        # MaximumLatitude : uly

        #push into 360 domain
        if (force360):
          if (ulx < 0):
             ulx = ulx + 360
          if (lrx < 0):
             lrx = lrx + 360
        # MinimumLongitude : urx)
        # MaximumLongitude : lrx)

        # UpperLeftCornerX in meters : UpperLeftCornerX
        # UpperLeftCornerY in meters : UpperLeftCornerY

        if debug:
			#testing print per band
            tofits.write('band {} read\n'.format(i))
			#test numpy, and how to Mask for stats using GDAL's nodata
            nodatamask = raster_data == dfNoData
            raster_data[nodatamask] = dfNoData
			#or another way
			#raster_data[raster_data == dfNoData] = np.nan
            tofits.write('band min: {}\n'.format(raster_data.min()))
            tofits.write('band max: {}\n'.format(raster_data.max()))
			

    #end read/write band loop
    tofits.write("complete\n")
    #tofits.close()
    tofits = None
    raster_data = None
    iBand = None
    inDataset = None
    return 0

#/************************************************************************/
#/*                        GDALInfoReportCorner()                        */
#/************************************************************************/

def GDALInfoReportCorner( inDataset, hTransform, corner_name, x, y ):

    line = "%-11s " % corner_name

#/* -------------------------------------------------------------------- */
#/*      Transform the point into georeferenced coordinates.             */
#/* -------------------------------------------------------------------- */
    adfGeoTransform = inDataset.GetGeoTransform(can_return_null = True)
    if adfGeoTransform is not None:
        dfGeoX = adfGeoTransform[0] + adfGeoTransform[1] * x \
            + adfGeoTransform[2] * y
        dfGeoY = adfGeoTransform[3] + adfGeoTransform[4] * x \
            + adfGeoTransform[5] * y

    else:
        line = line + ("(%7.1f,%7.1f)" % (x, y ))
        print(line)
        return False

#/* -------------------------------------------------------------------- */
#/*      Report the georeferenced coordinates.                           */
#/* -------------------------------------------------------------------- */
    if abs(dfGeoX) < 181 and abs(dfGeoY) < 91:
        line = line + ( "(%12.7f,%12.7f) " % (dfGeoX, dfGeoY ))

    else:
        line = line + ( "(%12.3f,%12.3f) " % (dfGeoX, dfGeoY ))

#/* -------------------------------------------------------------------- */
#/*      Transform to latlong and report.                                */
#/* -------------------------------------------------------------------- */
    if hTransform is not None:
        pnt = hTransform.TransformPoint(dfGeoX, dfGeoY, 0)
        if pnt is not None:
            line = line + ( "(%s," % gdal.DecToDMS( pnt[0], "Long", 2 ) )
            line = line + ( "%s)" % gdal.DecToDMS( pnt[1], "Lat", 2 ) )

    print(line)

    return True

#/************************************************************************/
#/*                        GDALGetLon()                              */
#/************************************************************************/
def GDALGetLon( inDataset, hTransform, x, y ):

#/* -------------------------------------------------------------------- */
#/*      Transform the point into georeferenced coordinates.             */
#/* -------------------------------------------------------------------- */
    adfGeoTransform = inDataset.GetGeoTransform(can_return_null = True)
    if adfGeoTransform is not None:
        dfGeoX = adfGeoTransform[0] + adfGeoTransform[1] * x \
            + adfGeoTransform[2] * y
        dfGeoY = adfGeoTransform[3] + adfGeoTransform[4] * x \
            + adfGeoTransform[5] * y
    else:
        return 0.0

#/* -------------------------------------------------------------------- */
#/*      Transform to latlong and report.                                */
#/* -------------------------------------------------------------------- */
    if hTransform is not None:
        pnt = hTransform.TransformPoint(dfGeoX, dfGeoY, 0)
        if pnt is not None:
          return pnt[0]
    return dfGeoX


#/************************************************************************/
#/*                        GDALGetLat()                              */
#/************************************************************************/

def GDALGetLat( inDataset, hTransform, x, y ):
#/* -------------------------------------------------------------------- */
#/*      Transform the point into georeferenced coordinates.             */
#/* -------------------------------------------------------------------- */
    adfGeoTransform = inDataset.GetGeoTransform(can_return_null = True)
    if adfGeoTransform is not None:
        dfGeoX = adfGeoTransform[0] + adfGeoTransform[1] * x \
            + adfGeoTransform[2] * y
        dfGeoY = adfGeoTransform[3] + adfGeoTransform[4] * x \
            + adfGeoTransform[5] * y
    else:
        return 0.0

#/* -------------------------------------------------------------------- */
#/*      Transform to latlong and report.                                */
#/* -------------------------------------------------------------------- */
    if hTransform is not None:
        pnt = hTransform.TransformPoint(dfGeoX, dfGeoY, 0)
        if pnt is not None:
          return pnt[1]
    return dfGeoY

if __name__ == '__main__':
    version_num = int(gdal.VersionInfo('VERSION_NUM'))
    if version_num < 1800: # because of GetGeoTransform(can_return_null)
        print('ERROR: Python bindings of GDAL 1.8.0 or later required')
        sys.exit(1)

    sys.exit(main(sys.argv))
        
