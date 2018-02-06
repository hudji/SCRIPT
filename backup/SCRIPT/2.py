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
f = open('C:/Users/hudjimartsu/Documents/PROJECT/FOREST 2020/TRAINING/PyQgis/DATA/source/mtl.txt', 'r') #open file for reading
def build_data(f): #build dictionary
    output = {} #Dict
    for line in f.readlines(): #Iterates through every line in the string
        if "=" in line: #make sure line has data as wanted
            l = line.split("=") #Seperate by "=" and put into a list
            output[l[0].strip()] = l[1].strip() #First word is key, second word is value
    return output #Returns a dictionary with the key, value pairs.
data = build_data(f)

#Load data raster
raster_list=glob.glob('C:/Users/hudjimartsu/Documents/PROJECT/FOREST 2020/TRAINING/PyQgis/DATA/source/*.TIF')
read=[]
for i in raster_list:
    band=gdal.Open(i)
    read.append(band.ReadAsArray())

coastal=array(read[0],dtype=float32)
blue=array(read[1], dtype=float32)
green=array(read[2],dtype=float32)   
red=array(read[3],dtype=float32)
nir=array(read[4],dtype=float32)
swir1=array(read[5],dtype=float32)
swir2=array(read[6], dtype=float32)
panch=array(read[7],dtype=float32)
cirrus=array(read[8],dtype=float32)
tirs1=array(read[9],dtype=float32)
tirs2=array(read[10],dtype=float32)

#Reflectance
reflectance_blue=(float(data['REFLECTANCE_MULT_BAND_2'])*blue+float(data['REFLECTANCE_ADD_BAND_2']))/sin(54.68707440*3.1459/180)
reflectance_green=(float(data['REFLECTANCE_MULT_BAND_3'])*green+float(data['REFLECTANCE_ADD_BAND_3']))/sin(54.68707440*3.1459/180)
reflectance_red=(float(data['REFLECTANCE_MULT_BAND_4'])*red+float(data['REFLECTANCE_ADD_BAND_4']))/sin(54.68707440*3.1459/180)
reflectance_nir=(float(data['REFLECTANCE_MULT_BAND_5'])*nir+float(data['REFLECTANCE_ADD_BAND_5']))/sin(54.68707440*3.1459/180)
reflectance_swir2=(float(data['REFLECTANCE_MULT_BAND_7'])*swir2+float(data['REFLECTANCE_ADD_BAND_7']))/sin(54.68707440*3.1459/180)

np.seterr(invalid='ignore')
#NDSI
#def ndsi():
    a=np.array(reflectance_blue - reflectance_nir)
    b=reflectance_blue + reflectance_nir
    result=a/b
 #   return result
    
ndsi=((reflectance_blue - reflectance_nir)/(reflectance_blue + reflectance_nir))

#NDVI
ndvi=(reflectance_red-reflectance_green)/(reflectance_red+reflectance_green)

# Step I


geo = band.GetGeoTransform()
proj = band.GetProjection()
shape = red.shape
driver = gdal.GetDriverByName("GTiff")
dst_ds = driver.Create("C:/Users/hudjimartsu/Documents/PROJECT/FOREST 2020/TRAINING/PyQgis/DATA/RESULT/REFLECTANCE/ref_ndsi_1.tif", shape[1], shape[0], 1, gdal.GDT_Float32)
#dst_ds = driver.Create("C:/Users/hudjimartsu/Documents/PROJECT/FOREST 2020/TRAINING/PyQgis/DATA/landsat/ndvi2.tif", shape [1], shape [0], 1, gdal.GDT_Float32)
dst_ds.SetGeoTransform(geo)
dst_ds.SetProjection(proj)
dst_ds.GetRasterBand(1).WriteArray(ndsi)
dst_ds = None  # save, close"""

#print 'nilai', band
print 'angka', ndvi
