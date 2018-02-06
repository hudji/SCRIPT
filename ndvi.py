from osgeo import gdal
import numpy as np
from numpy import *
g = gdal.Open("D:/PROJECT/FOREST 2020/TRAINING/PyQgis/DATA/source/LC81220652016230LGN00_B4.TIF")
red = g.ReadAsArray()
g = gdal.Open("D:/PROJECT/FOREST 2020/TRAINING/PyQgis/DATA/source/LC81220652016230LGN00_B5.TIF")
nir = g.ReadAsArray()

#
red = array(red, dtype = float)
nir = array(nir, dtype = float)
check = np.logical_and ( red > 1, nir > 1 )
ndvi = np.where ( check,  (nir - red ) / ( nir + red )*100, -999 ) 
geo = g.GetGeoTransform()  
proj = g.GetProjection()   
shape = red.shape        
driver = gdal.GetDriverByName("GTiff")
dst_ds = driver.Create("D:/PROJECT/FOREST 2020/TRAINING/PyQgis/DATA/RESULT/ndvi3.tif", shape[1], shape[0], 1, gdal.GDT_Float32)
dst_ds.SetGeoTransform(geo)
dst_ds.SetProjection(proj)
dst_ds.GetRasterBand(1).WriteArray(ndvi)
dst_ds = None  # save, close
print 'nilai',shape