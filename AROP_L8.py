import datetime
import math
from datetime import datetime, date
import numpy as np
import glob
import os
from osgeo import gdal
from scipy.stats import linregress
import pandas as pd
# from dict import dict
import numexpr


#Load data raster
path='D:/FORESTS2020/DATA/LANDSAT/LC08_L1TP_123064_20170827_20170913_01_T1.tar/LC08_L1TP_123064_20170827_20170913_01_T1/'
raster_list=glob.glob(path+ '*.TIF')
warp=[]
for i in raster_list:
    band_warp=gdal.Open(i)
    warp.append(band_warp.GetRasterBand(1).ReadAsArray())

band_base=gdal.Open("D:/FORESTS2020/DATA/LANDSAT/GEOCOVER/CLIP/G120065.TIF")
base= band_base.GetRasterBand(1).ReadAsArray()

#lat long value
# get columns and rows of your image from gdalinfo

def pixel2coord(x, y):
    xp = a * x + b * y + xoff
    yp = d * x + e * y + yoff
    return(xp, yp)

###Base Image
xoff, a, b, yoff, d, e = band_base.GetGeoTransform()
ulx_base= xoff
uly_base=yoff
rows=base.shape[0]
colms=base.shape[1]
coordinate_base=[]
for row in  range(0,rows):
  for col in  range(0,colms):
      coordinate_base.append(pixel2coord(col,row))
coor_2=np.array(coordinate_base, dtype=float)
long=coor_2[:,0]
lat=coor_2[:,1]
long_base=long.reshape(rows,colms)
lat_base=lat.reshape(rows,colms)
del xoff, a, b, yoff, d, e


###Warp Image
xoff, a, b, yoff, d, e = band_warp.GetGeoTransform()
ulx_warp= xoff
uly_warp=yoff
rowz=warp[7].shape[0]
colmz=warp[7].shape[1]
coordinate_warp=[]
for rown in  range(0,rowz):
  for coln in  range(0,colmz):
      coordinate_warp.append(pixel2coord(coln,rown))
coor_w=np.array(coordinate_warp, dtype=float)
long=coor_w[:,0]
lat=coor_w[:,1]
long_warp=long.reshape(rowz,colmz)
lat_warp=lat.reshape(rowz,colmz)

print (ulx_warp)
print (ulx_base)
### Automated Registration

wx= ((long_base*28.5)+(ulx_base-ulx_warp))/30
wy= ((lat_base*28.5)+(uly_base-uly_warp))/30

