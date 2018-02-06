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
import numexpr
#from __future__ import division

#Load Metadata
f = open('D:/PROJECT/FOREST 2020/TRAINING/PyQgis/DATA/clip/mtl.txt', 'r') #open file for reading
def build_data(f): #build dictionary
    output = {} #Dict
    for line in f.readlines(): #Iterates through every line in the string
        if "=" in line: #make sure line has data as wanted
            l = line.split("=") #Seperate by "=" and put into a list
            output[l[0].strip()] = l[1].strip() #First word is key, second word is value
    return output #Returns a dictionary with the key, value pairs.
data = build_data(f)

#Load data raster
raster_list=glob.glob('D:/PROJECT/FOREST 2020/TRAINING/PyQgis/DATA/clip/*.TIF')
read=[]
for i in raster_list:
    band=gdal.Open(i)
    read.append(band.GetRasterBand(1).ReadAsArray())

#coastal=np.array(read[0],dtype=float)
blue=np.array(read[3], dtype=float)
green=np.array(read[4],dtype=float)   
red=np.array(read[5],dtype=float)
nir=np.array(read[6],dtype=float)
swir1=np.array(read[7],dtype=float)
swir2=np.array(read[8], dtype=float)
#panch=np.array(read[9],dtype=float)
#cirrus=np.array(read[10],dtype=float)
tirs1=np.array(read[1],dtype=float)
tirs2=np.array(read[2],dtype=float)

#Reflectance
reflectance_blue=(float(data['REFLECTANCE_MULT_BAND_1'])*blue+float(data['REFLECTANCE_ADD_BAND_1']))/sin(54.68707440*3.1459/180)
reflectance_green=(float(data['REFLECTANCE_MULT_BAND_2'])*green+float(data['REFLECTANCE_ADD_BAND_2']))/sin(54.68707440*3.1459/180)
reflectance_red=(float(data['REFLECTANCE_MULT_BAND_3'])*red+float(data['REFLECTANCE_ADD_BAND_3']))/sin(54.68707440*3.1459/180)
reflectance_nir=(float(data['REFLECTANCE_MULT_BAND_5'])*nir+float(data['REFLECTANCE_ADD_BAND_5']))/sin(54.68707440*3.1459/180)
reflectance_swir2=(float(data['REFLECTANCE_MULT_BAND_7'])*swir2+float(data['REFLECTANCE_ADD_BAND_7']))/sin(54.68707440*3.1459/180)

###Step I
#NDVI Formula
n=nir
r=red
"""np.seterr(divide='ignore', invalid='ignore') #Ignore the divided by zero or Nan appears
upper=n-r
lower=n+r
mask = (upper == lower) & (upper==0)
ndvi=np.where(mask, 0, upper/lower)"""
NDVI = numexpr.evaluate("(n - r) / (n + r)")
#NDSI Formula
"""np.seterr(divide='ignore', invalid='ignore')
upper2=green-swir1
lower2=green+swir1
mask2=(upper2 == lower2) & (upper2 ==0)
ndsi=np.where(mask2, 0, upper2/lower2)"""
NDSI=numexpr.evaluate("(green - swir1) / (green + swir1)")

#Brightness Temperature
radiance_tirs1=(float(data['RADIANCE_MULT_BAND_10'])*tirs1+float(data['RADIANCE_ADD_BAND_10']))
bt=float(data['K2_CONSTANT_BAND_10'])/np.log(float(data['K1_CONSTANT_BAND_10'])/radiance_tirs1+1)
bt_n=bt-273.15
#DataFrame
"""
row = ndvi.shape[0]
col = ndvi.shape[1]
ndvi2 = ndvi.ravel()
reflectance_swir2n = reflectance_swir2.ravel()
ndsi2=ndsi.ravel()
bt_n2=bt_n.ravel()
d = {'band7':reflectance_swir2n, 'BT':bt_n2, 'NDSI':ndsi2, 'NDVI':ndvi2}
dataframe = pd.DataFrame(d)
E=dataframe.query('band7 > 0.03 and BT < 27 and NDSI < 0.8 and NDVI < 0.8')
#truee = E.index
Hasil = np.repeat('False',len(ndvi2))
Hasil[E.index] = 'True'
basic_test = Hasil.reshape(row,col)"""
idplcd = numexpr.evaluate("(NDSI < 0.8) & (NDVI < 0.8) & (reflectance_swir2 > 0.03) & (bt_n < 27)")

#Mean Vis
"""
vis=[blue, green, red]
MeanVis=(blue+green+red)/3
Mean_test=0
for a in vis:
    upper3=(a-MeanVis)
    mask3=(upper3 == MeanVis) & (upper3 ==0)
    test=abs(np.where(mask3, 0, upper3/MeanVis))
    Mean_test +=test
whiteness= Mean_test < 0.7  
whiteness2=whiteness.ravel()"""
visimean = numexpr.evaluate("(blue + green + red) / 3 ")
whiteness = numexpr.evaluate("(abs(blue - visimean) + abs(green - visimean)+ abs(red - visimean)) / visimean")
del visimean

# update idplcd
#whiteness[satu_Bv] = 0  # If one visible is saturated whiteness == 0
idplcd &= whiteness < 0.7

# Haze test
HOT = numexpr.evaluate("reflectance_blue- 0.5 * reflectance_red - 0.08")  # Haze test
idplcd &= numexpr.evaluate("(HOT > 0)")
del HOT  # need to find thick warm cloud

# Ratio4/5>0.75 cloud test
idplcd &= numexpr.evaluate("(nir / swir1) > 0.75")

#water test
WT = numexpr.evaluate("((NDVI < 0.01) & (nir < 0.11)) | ((NDVI < 0.1) & (NDVI > 0) & (nir < 0.05))")




"""
geo = band.GetGeoTransform()
proj = band.GetProjection()
shape = red.shape
driver = gdal.GetDriverByName("GTiff")
dst_ds = driver.Create("D:/PROJECT/FOREST 2020/TRAINING/PyQgis/DATA/RESULT/REFLECTANCE/pass_one.tif", shape[1], shape[0], 1, gdal.GDT_Float32)
dst_ds.SetGeoTransform(geo)
dst_ds.SetProjection(proj)
dst_ds.GetRasterBand(1).WriteArray(pass_one)
dst_ds = None  # save, close"""

#print 'nilai', band
#print 'angka', ndvi




