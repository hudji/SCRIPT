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
import time
start_time = time.time()
f = open('D:/FORESTS2020/TRAINING/PyQgis/DATA/Landsat8/mtl.txt', 'r') #open file for reading
def build_data(f): #build dictionary
    output = {} #Dict
    for line in f.readlines(): #Iterates through every line in the string
        if "=" in line: #make sure line has data as wanted
            l = line.split("=") #Seperate by "=" and put into a list
            output[l[0].strip()] = l[1].strip() #First word is key, second word is value
    return output #Returns a dictionary with the key, value pairs.
data = build_data(f)
print "Loading data raster..."
#Load data raster
raster_list=glob.glob('D:/FORESTS2020/TRAINING/PyQgis/DATA/Landsat8/NO SEA/FINAL/GEO/LANDSAT/TRIAL/*.TIF')
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

def export_array(in_array, output_path):
    """This function is used to produce output of array as a map."""
    global proj, geotrans, row, col
    proj     = band.GetProjection()
    geotrans = band.GetGeoTransform()
    row      = band.RasterYSize
    col      = band.RasterXSize
    driver   = gdal.GetDriverByName("GTiff")
    outdata  = driver.Create(output_path, col, row, 1)
    outband  = outdata.GetRasterBand(1)
    outband.SetNoDataValue(-9999)
    outband.WriteArray(in_array)
    # Georeference the image
    outdata.SetGeoTransform(geotrans)
    # Write projection information
    outdata.SetProjection(proj)
    outdata.FlushCache()
    outdata = None
def cos(x):
    cos= np.cos(np.deg2rad(x))
    return  cos
def sin(x):
    sin=np.sin(np.deg2rad(x))
    return sin

#Load data raster aspect, slope & sample area
raster_list_dem=glob.glob('D:/FORESTS2020/TRAINING/PyQgis/DATA/Landsat8/NO SEA/FINAL/GEO/DEM/*.tif')
read2=[]
for d in raster_list_dem:
    band2=gdal.Open(d)
    read2.append(band2.GetRasterBand(1).ReadAsArray())

dem_name=['aspect', 'sample', 'slope']
dem_dict= dict (zip(dem_name, read2))
print "Calculating Solar Position angle..."
def year_date():
    year_file=data['DATE_ACQUIRED']
    date_file=data['SCENE_CENTER_TIME']
    date_file2= date_file [1:16]
    all= year_file+" "+date_file2
    parsing = datetime.strptime(all, '%Y-%m-%d %H:%M:%S.%f')
    return parsing
dt=year_date()

def hour():
    h=dt.hour+7
    return h
def second():
    s= float(dt.microsecond)/1000000+dt.second
    return s
def leap():
    if (dt.year % 4) == 0:
        if (dt.year % 100) == 0:
            if (dt.year % 400) == 0:
               a = int(366)
            else:
                a = int(365)
        else:
            a= int(366)
    else:
        a= int(365)
    return a
def cos(x):
    cos= np.cos(np.deg2rad(x))
    return  cos
def sin(x):
    sin=np.sin(np.deg2rad(x))
    return sin
def day():
    day_date= date(dt.year, dt.month, dt.day)
    sum_of_day=int(day_date.strftime('%j'))
    return sum_of_day

gamma=((2 * math.pi) / leap()) * ((day() - 1) + (((hour()+dt.minute/60+second()/3600) - 12) / 24) )# degree


#sun declination angle
decl=0.006918 - 0.399912 * cos(gamma) + 0.070257 * sin(gamma) - 0.006758 * cos (2 * gamma)\
     + 0.000907 * sin (2 * gamma) - 0.002697 * cos (3 * gamma) + 0.00148 * sin (3 * gamma) #radians
decl_deg= (360 / (2 * math.pi)) * decl

#lat long value
# get columns and rows of your image from gdalinfo

xoff, a, b, yoff, d, e = band.GetGeoTransform()
def pixel2coord(x, y):
    xp = a * x + b * y + xoff
    yp = d * x + e * y + yoff
    return(xp, yp)
rows=read[0].shape[0]
colms=read[0].shape[1]
coordinate=[]
for row in  range(0,rows):
  for col in  range(0,colms):
      coordinate.append(pixel2coord(col,row))
coor_2=np.array(coordinate, dtype=float)
long=coor_2[:,0]
lat=coor_2[:,1]
long_n=long.reshape(rows,colms)
lat_n=lat.reshape(rows,colms)

#eqtime
eqtime = 229.18 * (0.000075 + 0.001868 * cos(gamma) - 0.032077 * sin(gamma) - 0.014615 * cos(2 * gamma) - 0.040849 * sin(2 * gamma))  # minutes
timeoff= eqtime - 4 * long_n + 60 * 7 #minutes
tst=hour() * 60 + dt.minute + second() / 60 + timeoff #minutes
ha=(tst /4)-180 #degree

#sun zenith angle
zenit1 =sin(lat_n)* sin(decl_deg) + cos (lat_n)* cos(decl_deg) * cos(ha)
zenit2=np.arccos(zenit1) #radians
zenit_angle= np.rad2deg(zenit2)

print "Calculating Reflectances...."
#Reflectance
reflectance_blue=(float(data['REFLECTANCE_MULT_BAND_2'])*blue+float(data['REFLECTANCE_ADD_BAND_2']))/cos(zenit_angle)
reflectance_red=(float(data['REFLECTANCE_MULT_BAND_4'])*red+float(data['REFLECTANCE_ADD_BAND_4']))/cos(zenit_angle)
reflectance_nir=(float(data['REFLECTANCE_MULT_BAND_5'])*nir+float(data['REFLECTANCE_ADD_BAND_5']))/cos(zenit_angle)
reflectance_swir1=(float(data['REFLECTANCE_MULT_BAND_6'])*swir1+float(data['REFLECTANCE_ADD_BAND_6']))/cos(zenit_angle)
reflectance_swir2=(float(data['REFLECTANCE_MULT_BAND_7'])*swir2+float(data['REFLECTANCE_ADD_BAND_7']))/cos(zenit_angle)
reflectance_green=(float(data['REFLECTANCE_MULT_BAND_3'])*green+float(data['REFLECTANCE_ADD_BAND_3']))/cos(zenit_angle)

#Brightness Temperature
radiance_tirs1=(float(data['RADIANCE_MULT_BAND_10'])*tirs1+float(data['RADIANCE_ADD_BAND_10']))
bt_n=float(data['K2_CONSTANT_BAND_10'])/np.log(float(data['K1_CONSTANT_BAND_10'])/radiance_tirs1+1)
BT=bt_n-273.15
###Step I
print "Calculating BASIC FMASK ..."
#Equation 1
#NDVI
NDVI=numexpr.evaluate("(reflectance_nir - reflectance_red) / (reflectance_nir + reflectance_red)")
#NDSI
NDSI=numexpr.evaluate("(reflectance_green-reflectance_swir1)/(reflectance_green+reflectance_swir1)")

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
print"Combining BQA & FMASK..."
#bqa
prob_cloud= numexpr.evaluate("((bqa >= 53248))") # https://landsat.usgs.gov/qualityband

print"Cleaning & Filtering..."
#FMVsBQ
prob_cloud2= numexpr.evaluate("((prob_cloud==True)|(basic_test==True))")
filtering = ndimage.median_filter(prob_cloud2, 3)
#Cleaning
#cleaning_open = ndimage.binary_opening(filtering)
cleaning = ndimage.binary_erosion(filtering)
filtering2 = ndimage.median_filter(cleaning, 2)

export_array(filtering2, "D:/FORESTS2020/TRAINING/PyQgis/RESULT/Landsat8/Cloud_Masking/C121217/Cloud_masking_III.tif")
os.system('gdalwarp -dstnodata -9999 -cutline D:/FORESTS2020/TRAINING/PyQgis/DATA/source/TEST/no_sea.shp D:/FORESTS2020/TRAINING/PyQgis/RESULT/Landsat8/Cloud_Masking/C121217/Cloud_masking_III.TIF D:/FORESTS2020/TRAINING/PyQgis/RESULT/Landsat8/Cloud_Masking/C121217/Cloud_masking_III_CL.TIF')
print("Finish: running time: %s seconds" % (time.time() - start_time))
