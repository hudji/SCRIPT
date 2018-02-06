from osgeo import gdal
import numpy as np
from numpy import *
import glob
#import codecs
raster_list=glob.glob('C:/Users/hudjimartsu/Documents/PROJECT/FOREST 2020/TRAINING/PyQgis/DATA/source/*.TIF')
read=[]
for i in raster_list:
    band=gdal.Open(i)
    read.append(band.ReadAsArray())
#g=band[3]
#g=band[4]
red=array(read[3], dtype=float32)
nir=array(read[4], dtype=float)
check = np.logical_and ( red > 1, nir > 1 )
ndvi = np.where ( check,  (nir - red ) / ( nir + red ), -999 )
geo = band.GetGeoTransform()
proj = band.GetProjection()
shape = red.shape
driver = gdal.GetDriverByName("GTiff")
dst_ds = driver.Create("C:/Users/hudjimartsu/Documents/PROJECT/FOREST 2020/TRAINING/PyQgis/DATA/landsat/ndvi4.tif", shape[1], shape[0], 1, gdal.GDT_Float32)
#dst_ds = driver.Create("C:/Users/hudjimartsu/Documents/PROJECT/FOREST 2020/TRAINING/PyQgis/DATA/landsat/ndvi2.tif", shape [1], shape [0], 1, gdal.GDT_Float32)
dst_ds.SetGeoTransform(geo)
dst_ds.SetProjection(proj)
dst_ds.GetRasterBand(1).WriteArray(ndvi)
dst_ds = None  # save, close"""

#print 'nilai', band
print 'angka', shape

