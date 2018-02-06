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
import numexpr,sys
from datetime import datetime, date
from scipy import ndimage
from subprocess import call
import os


#gdalwarp2=r'C:/Program Files/GDAL/gdalwarp.exe'

f = open('D:/FORESTS2020/TRAINING/PyQgis/DATA/Landsat8/mtl.txt', 'r') #open file for reading
def build_data(f): #build dictionary
    output = {} #Dict
    for line in f.readlines(): #Iterates through every line in the string
        if "=" in line: #make sure line has data as wanted
            l = line.split("=") #Seperate by "=" and put into a list
            output[l[0].strip()] = l[1].strip() #First word is key, second word is value
    return output #Returns a dictionary with the key, value pairs.
data = build_data(f)

#Load data raster
raster_list=glob.glob('D:/FORESTS2020/TRAINING/PyQgis/DATA/source/*.TIF')
read=[]
for i in raster_list:
    band=gdal.Open(i)
    read.append(band.GetRasterBand(1).ReadAsArray())

#coastal=np.array(read[0],dtype=float)
blue=np.array(read[1], dtype=float)
green=np.array(read[2],dtype=float)
red=np.array(read[3],dtype=float)
nir=np.array(read[4],dtype=float)
swir1=np.array(read[5],dtype=float)
swir2=np.array(read[6],dtype=float)
tirs1=np.array(read[0],dtype=float)
bqa=np.array(read[7],dtype=float)
#bqa[bqa > 64512] = 0
def export_array(in_array, output_path):
    """This function is used to produce output of array as a map."""
    global proj, geotrans, row, col
    proj     = band.GetProjection()
    geotrans = band.GetGeoTransform()
    row      = band.RasterYSize
    col      = band.RasterXSize
    driver   = gdal.GetDriverByName("GTiff")
    outdata  = driver.Create(output_path, col, row, 1, gdal.GDT_Float32)
    outband  = outdata.GetRasterBand(1)
    outband.SetNoDataValue(np.nan)
    outband.WriteArray(in_array)
    # Georeference the image
    outdata.SetGeoTransform(geotrans)
    # Write projection information
    outdata.SetProjection(proj)
    outdata.FlushCache()
    outdata = None

def raster_buffer(raster_filepath, dist=100):
    """This function creates a distance buffer around the given raster file with non-zero values.
    The value in output raster will have value of the cell to which it is close to."""
    d = gdal.Open(raster_filepath)
    if d is None:
        print("Error: Could not open image " + raster_filepath)
        sys.exit(1)
    inband      = d.GetRasterBand(1)
    in_array    = inband.ReadAsArray(0, 0, col, row).astype(int)
    Xcell_size  = int(abs(geotrans[1]))
    Ycell_size  = int(abs(geotrans[5]))
    cell_size   = (Xcell_size + Ycell_size) / 2
    cell_dist   = dist / cell_size
    in_array[in_array == (inband.GetNoDataValue() or 0 or -999)] = 0
    out_array   = np.zeros_like(in_array)
    temp_array  = np.zeros_like(in_array)
    i, j, h, k  = 0, 0, 0, 0
    print("Running distance buffer, please wait...")
    while (h < col):
        k = 0
        while (k < row):
            if (in_array[k][h] >= 1):
                i = h - cell_dist
                while ((i < cell_dist + h) and i < col):
                    j = k - cell_dist
                    while (j < (cell_dist + k) and j < row):
                        if (((i - h) ** 2 + (j - k) ** 2) <= cell_dist ** 2):
                            if (temp_array[j][i] == 0 or temp_array[j][i] > ((i - h) ** 2 + (j - k) ** 2)):
                                out_array[j][i] = in_array[k][h]
                                temp_array[j][i] = (i - h) ** 2 + (j - k) ** 2
                        j += 1
                    i += 1
            k += 1
        h += 1
    d, temp_array, in_array = None, None, None
    return out_array

def cos(x):
    cos= np.cos(np.deg2rad(x))
    return  cos
def sin(x):
    sin=np.sin(np.deg2rad(x))
    return sin


#Reflectance
reflectance_blue=(float(data['REFLECTANCE_MULT_BAND_2'])*blue+float(data['REFLECTANCE_ADD_BAND_2']))/sin(54.68707440)
reflectance_red=(float(data['REFLECTANCE_MULT_BAND_4'])*red+float(data['REFLECTANCE_ADD_BAND_4']))/sin(54.68707440)
reflectance_nir=(float(data['REFLECTANCE_MULT_BAND_5'])*nir+float(data['REFLECTANCE_ADD_BAND_5']))/sin(54.68707440)
reflectance_swir1=(float(data['REFLECTANCE_MULT_BAND_6'])*swir1+float(data['REFLECTANCE_ADD_BAND_6']))/sin(54.68707440)
reflectance_swir2=(float(data['REFLECTANCE_MULT_BAND_7'])*swir2+float(data['REFLECTANCE_ADD_BAND_7']))/sin(54.68707440)
reflectance_green=(float(data['REFLECTANCE_MULT_BAND_3'])*green+float(data['REFLECTANCE_ADD_BAND_3']))/sin(54.68707440)

#Brightness Temperature
radiance_tirs1=(float(data['RADIANCE_MULT_BAND_10'])*tirs1+float(data['RADIANCE_ADD_BAND_10']))
bt_n=float(data['K2_CONSTANT_BAND_10'])/np.log(float(data['K1_CONSTANT_BAND_10'])/radiance_tirs1+1)
BT=bt_n-273.15
###Step I
#Equation 1
#NDVI
NDVI=numexpr.evaluate("(reflectance_nir - reflectance_red) / (reflectance_nir + reflectance_red)")
#NDSI
NDSI=numexpr.evaluate("(reflectance_green-reflectance_swir1)/(reflectance_green+reflectance_swir1)")

#NDVI[numexpr.evaluate("(reflectance_nir + reflectance_red) == 0")] = 0.01
#NDSI[numexpr.evaluate("(reflectance_green + reflectance_swir1) == 0")] = 0.01
basic_test=numexpr.evaluate( "(reflectance_swir2 > 0.3)& (NDSI < 0.8) & (NDVI < 0.8) & (BT < 27) ")

#Equation 2
#Mean Vis
meanVIS = numexpr.evaluate("(reflectance_blue + reflectance_green + reflectance_red) / 3 ")
whiteness = numexpr.evaluate("(abs(reflectance_blue - meanVIS) + abs(reflectance_green - meanVIS)+ abs(reflectance_red - meanVIS)) / meanVIS")
basic_test&= whiteness < 0.7
del meanVIS

#Equation 3
# Haze test
hazeTest = (numexpr.evaluate("reflectance_blue- 0.5 * reflectance_red - 0.08"))> 0  # Haze test
basic_test &= numexpr.evaluate("(hazeTest > 0)")
#del HOT  # need to find thick warm cloud

#Equation 4
# Ratio4/5>0.75 cloud test
basic_test &= numexpr.evaluate("(reflectance_nir / reflectance_swir1) > 0.75") # Equation 6

#bqa
Cloud_cirrus2 = numexpr.evaluate("((bqa >= 53248))")
#prob_bqa= prob_cloud[numexpr.evaluate('prob_cloud > 65534')] = -9999 # https://landsat.usgs.gov/qualityband

#FMVsBQ
prob_cloud2= numexpr.evaluate("((prob_cloud==True)|(basic_test==True))")

#Cleaning
cleaning = ndimage.binary_erosion(prob_cloud2)

#Median filtering
cloud = ndimage.median_filter(cleaning, 5)
export_array(cleaning, "D:/FORESTS2020/TRAINING/PyQgis/RESULT/Landsat8/Cloud_Masking/C121217/cleaning_I.tif")
#Buffering
raster_buffer_array = raster_buffer('D:/FORESTS2020/TRAINING/PyQgis/RESULT/Landsat8/Cloud_Masking/C121217/proud_cloud2_I.tif',100)
export_array(raster_buffer_array, "D:/FORESTS2020/TRAINING/PyQgis/RESULT/Landsat8/Cloud_Masking/C111217/Buf100_Fi3_II.tif")
export_array(Cloud_cirrus2,"D:/FORESTS2020/TRAINING/PyQgis/RESULT/Landsat8/Cloud_Masking/C121217/basic_prob_I.tif")

