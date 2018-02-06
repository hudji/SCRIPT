# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 13:41:34 2017

@author: sahid
"""

# NDVI Python Script
#
# GNU GENERAL PUBLIC LICENSE
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# Created by Alexandros Falagas.
#
from osgeo import gdal
# this allows GDAL to throw Python Exceptions
gdal.UseExceptions()
from gdalconst import *
import numpy as np
import sys
red=gdal.Open(r"D:/PROJECT/FOREST 2020/TRAINING/PyQgis/DATA/source/LC81220652016230LGN00_B4.TIF", GA_ReadOnly)
r=np.array(red.GetRasterBand(1).ReadAsArray(), dtype=float)
nir=gdal.Open(r"D:/PROJECT/FOREST 2020/TRAINING/PyQgis/DATA/source/LC81220652016230LGN00_B5.TIF", GA_ReadOnly)
n=np.array(nir.GetRasterBand(1).ReadAsArray(), dtype=float)
geotr=red.GetGeoTransform()
proj=red.GetProjection()
tableshape=r.shape
np.seterr(divide='ignore', invalid='ignore') #Ignore the divided by zero or Nan appears
ndvi=(n-r)/(n+r).astype(np.float64) # The NDVI formula
driver=gdal.GetDriverByName('GTiff')
dst_ds=driver.Create(r"D:/PROJECT/FOREST 2020/TRAINING/PyQgis/DATA/RESULT/NDVI_3.TIF", tableshape[1], tableshape[0], 1, gdal.GDT_Float32)
dst_ds.SetGeoTransform(geotr)
dst_ds.SetProjection(proj)
dst_ds.GetRasterBand(1).WriteArray(ndvi)
dst_ds=None # save, close
print 'The NDVI image is saved.'