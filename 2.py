# -*- coding: utf-8 -*-
"""
Created on Thu Sep 07 10:06:36 2017

@author: hudjimartsu
"""

from osgeo import gdal
import numpy as np
from numpy import *
import glob
import pandas as pd
#from __future__ import division

#Load Metadata
f = open('D:/PROJECT/FOREST 2020/TRAINING/PyQgis/DATA/source/mtl.txt', 'r') #open file for reading
def build_data(f): #build dictionary
    output = {} #Dict
    for line in f.readlines(): #Iterates through every line in the string
        if "=" in line: #make sure line has data as wanted
            l = line.split("=") #Seperate by "=" and put into a list
            output[l[0].strip()] = l[1].strip() #First word is key, second word is value
    return output #Returns a dictionary with the key, value pairs.
data = build_data(f)

#Load data raster
raster_list=glob.glob('D:/PROJECT/FOREST 2020/TRAINING/PyQgis/DATA/source/*.TIF')
read=[]
for i in raster_list:
    band=gdal.Open(i)
    read.append(band.GetRasterBand(1).ReadAsArray())

coastal=np.array(read[0],dtype=float)
blue=np.array(read[1], dtype=float)
green=np.array(read[2],dtype=float)   
red=np.array(read[3],dtype=float)
nir=np.array(read[4],dtype=float)
swir1=np.array(read[5],dtype=float)
swir2=np.array(read[6], dtype=float)
panch=np.array(read[7],dtype=float)
cirrus=np.array(read[8],dtype=float)
tirs1=np.array(read[9],dtype=float)
tirs2=np.array(read[10],dtype=float)

#Reflectance
reflectance_blue=(float(data['REFLECTANCE_MULT_BAND_2'])*blue+float(data['REFLECTANCE_ADD_BAND_2']))/sin(54.68707440*3.1459/180)
reflectance_green=(float(data['REFLECTANCE_MULT_BAND_3'])*green+float(data['REFLECTANCE_ADD_BAND_3']))/sin(54.68707440*3.1459/180)
reflectance_red=(float(data['REFLECTANCE_MULT_BAND_4'])*red+float(data['REFLECTANCE_ADD_BAND_4']))/sin(54.68707440*3.1459/180)
reflectance_nir=(float(data['REFLECTANCE_MULT_BAND_5'])*nir+float(data['REFLECTANCE_ADD_BAND_5']))/sin(54.68707440*3.1459/180)
reflectance_swir2=(float(data['REFLECTANCE_MULT_BAND_7'])*swir2+float(data['REFLECTANCE_ADD_BAND_7']))/sin(54.68707440*3.1459/180)
reflectance_panch=(float(data['REFLECTANCE_MULT_BAND_8'])*panch+float(data['REFLECTANCE_ADD_BAND_8']))/sin(54.68707440*3.1459/180)



np.seterr(divide='ignore', invalid='ignore') #Ignore the divided by zero or Nan appears
ndvi=(n-r)/(n+r).astype(np.float32) # The NDVI formula    

#NDVI

# Step I


geo = band.GetGeoTransform()
proj = band.GetProjection()
shape = red.shape
driver = gdal.GetDriverByName("GTiff")
dst_ds = driver.Create("D:/PROJECT/FOREST 2020/TRAINING/PyQgis/DATA/RESULT/REFLECTANCE/ndvi_2.tif", shape[1], shape[0], 1, gdal.GDT_Float32)
#dst_ds = driver.Create("C:/Users/hudjimartsu/Documents/PROJECT/FOREST 2020/TRAINING/PyQgis/DATA/landsat/ndvi2.tif", shape [1], shape [0], 1, gdal.GDT_Float32)
dst_ds.SetGeoTransform(geo)
dst_ds.SetProjection(proj)
dst_ds.GetRasterBand(1).WriteArray(ndvi)
dst_ds = None  # save, close"""

#print 'nilai', band
#print 'angka', ndvi
