from osgeo import gdal
import numpy

dst_filename = "D:\PROJECT\FOREST 2020\TRAINING\PyQgis\RESULT\Landsat8\L301017\Band6R_IC.tif"
#output to special GDAL "in memory" (/vsimem) path just for testing
#dst_filename = '/vsimem/test.tif'

#Raster size
raster_list=glob.glob('D:\PROJECT\FOREST 2020\TRAINING\PyQgis\DATA\Landsat8\NEW\*.tif')
read=[]
for i in raster_list:
    band=gdal.Open(i)
    read.append(band.GetRasterBand(1).ReadAsArray())

nrows=read[0].shape[0]
ncols=read[0].shape[1]
nbands=len(read)

#min & max random values of the output raster
zmin=0
zmax=12345

## See http://gdal.org/python/osgeo.gdal_array-module.html#codes
## for mapping between gdal and numpy data types
gdal_datatype = gdal.GDT_UInt32
np_datatype = numpy.uint32
driver = gdal.GetDriverByName( "GTiff" )
dst_ds = driver.Create( dst_filename, ncols, nrows, nbands, gdal_datatype )

## These are only required if you wish to georeference (http://en.wikipedia.org/wiki/Georeference)
## your output geotiff, you need to know what values to input, don't just use the ones below
#Coordinates of the lower left corner of the image
#in same units as spatial reference
#xllcorner=147.2
#yllcorner=-34.54

#Cellsize in same units as spatial reference
#cellsize=0.01

#dst_ds.SetGeoTransform( [ xllcorner, cellsize, 0, yllcorner, 0, -cellsize ] )
#srs = osr.SpatialReference()
#srs.SetWellKnownGeogCS("WGS84")
#dst_ds.SetProjection( srs.ExportToWkt() )

raster = numpy.random.randint(zmin,zmax, (nbands, nrows, ncols)).astype(np_datatype )
for band in range(nbands):
    dst_ds.GetRasterBand(band+1).WriteArray( raster[band, :, :] )

# Once we're done, close properly the dataset
dst_ds = None