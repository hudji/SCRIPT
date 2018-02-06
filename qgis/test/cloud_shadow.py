import numpy as np
import glob
import os
from osgeo import gdal
import pandas as pd
import numexpr
from dict import dict
f = open('D:\FORESTS2020\TRAINING\PyQgis\DATA\source\mtl.txt', 'r') #open file for reading
def build_data(f):
    output = {}
    for line in f.readlines():
        if "=" in line:
            l = line.split("=")
            output[l[0].strip()] = l[1].strip()
    return output
data = build_data(f)
#Load data raster
raster_list=glob.glob('D:\FORESTS2020\TRAINING\PyQgis\DATA\source\*.tif')
read=[]
for i in raster_list:
    band=gdal.Open(i)
    read.append(band.GetRasterBand(1).ReadAsArray())

bqa= np.array(read[11], dtype=float)
Cloud_cirrus2 = numexpr.evaluate("((bqa >= 53248))") # https://landsat.usgs.gov/qualityband

#export auto
geo = band.GetGeoTransform()
proj = band.GetProjection()
shape = bqa.shape
driver = gdal.GetDriverByName("GTiff")
dst_ds = driver.Create("D:\FORESTS2020\TRAINING\PyQgis\RESULT\Landsat8\Cloud Masking\C051217\Cloud.tif", shape[1], shape[0], 1, gdal.GDT_Float32)
dst_ds.SetGeoTransform(geo)
dst_ds.SetProjection(proj)
dst_ds.GetRasterBand(1).WriteArray(Cloud_cirrus2)
dst_ds.FlushCache()
dst_ds = None  # save, close"""