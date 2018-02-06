from osgeo import gdal
import numpy as np
from numpy import *
g = gdal.Open("C:/Users/hudjimartsu/Documents/PROJECT/FOREST 2020/TRAINING/PyQgis/DATA/landsat/LT52240631988227CUB02_B4.TIF")
red = g.ReadAsArray()
g = gdal.Open("C:/Users/hudjimartsu/Documents/PROJECT/FOREST 2020/TRAINING/PyQgis/DATA/landsat/LT52240631988227CUB02_B5.TIF")
nir = g.ReadAsArray()
red = array(red, dtype = float)
nir = array(nir, dtype = float)
check = np.logical_and ( red > 1, nir > 1 )
ndvi = np.where ( check,  (nir - red ) / ( nir + red )*100, -999 ) 
geo = g.GetGeoTransform()  
proj = g.GetProjection()   
shape = red.shape        
driver = gdal.GetDriverByName("GTiff")
dst_ds = driver.Create("C:/Users/hudjimartsu/Documents/PROJECT/FOREST 2020/TRAINING/PyQgis/DATA/landsat/ndvi.tif", shape[1], shape[0], 1, gdal.GDT_Float32)
dst_ds.SetGeoTransform(geo)
dst_ds.SetProjection(proj)
dst_ds.GetRasterBand(1).WriteArray(ndvi)
dst_ds = None  # save, close
print 'nilai',shape