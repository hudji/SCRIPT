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
band4=np.array(read[5], dtype=float)
band5=np.array(read[6], dtype=float)

reflectance_red=(float(data['REFLECTANCE_MULT_BAND_4'])*band4+float(data['REFLECTANCE_ADD_BAND_4']))/np.sin(54.68707440*3.1459/180)
reflectance_nir=(float(data['REFLECTANCE_MULT_BAND_5'])*band5+float(data['REFLECTANCE_ADD_BAND_5']))/np.sin(54.68707440*3.1459/180)
#cloud= bqa >= 53248

ndvi=numexpr.evaluate("(reflectance_nir - reflectance_red)/(reflectance_nir + reflectance_red)")
test_ndvi= numexpr.evaluate("(ndvi <= 0.1)")
Cloud_cirrus = numexpr.evaluate("((bqa >= 53248) | (bqa == 28672)) | ((bqa == 31744) | (bqa == 48128) & (ndvi<=0.1))") # https://landsat.usgs.gov/qualityband
Cloud_cirrus2 = numexpr.evaluate("((bqa >= 53248) | (bqa == 28672)) | ((bqa == 31744) | (bqa == 48128))") # https://landsat.usgs.gov/qualityband
(uniform_filter(cloudmask*2.0, size=3) >= 1.0)
#export auto
geo = band.GetGeoTransform()
proj = band.GetProjection()
shape = bqa.shape
driver = gdal.GetDriverByName("GTiff")
dst_ds = driver.Create("D:\FORESTS2020\TRAINING\PyQgis\RESULT\Landsat8\Cloud Masking\C041217\Cloud5.tif", shape[1], shape[0], 1, gdal.GDT_Float32)
dst_ds.SetGeoTransform(geo)
dst_ds.SetProjection(proj)
dst_ds.GetRasterBand(1).WriteArray(Cloud_cirrus2)
dst_ds.FlushCache()
dst_ds = None  # save, close"""