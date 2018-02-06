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
# from subprocess import call
import os
import time
start_time = time.time()
f = open('D:/FORESTS2020/REFERENSI/LAPAN/CLoud_masking/Landsat/120065/LC81200652016120LGN00_MTL.txt', 'r') #open file for reading
def build_data(f): #build dictionary
    output = {} #Dict
    for line in f.readlines(): #Iterates through every line in the string
        if "=" in line: #make sure line has data as wanted
            l = line.split("=") #Seperate by "=" and put into a list
            output[l[0].strip()] = l[1].strip() #First word is key, second word is value
    return output #Returns a dictionary with the key, value pairs.
data = build_data(f)
print "Loading Data Raster..."
#Load data raster
path='D:/FORESTS2020/REFERENSI/LAPAN/CLoud_masking/Landsat/120065/'
raster_list=glob.glob(path+ '*.TIF')
read=[]
for i in raster_list:
    band=gdal.Open(i)
    read.append(band.GetRasterBand(1).ReadAsArray())
filename=[]
for a in [os.path.basename(x) for x in glob.glob(path + '*.TIF')]:
    p=os.path.splitext(a)[0]
    filename.append(p)
my_dict= dict(zip(filename, read))

#Load data raster aspect, slope & sample area
pathname='D:/FORESTS2020/REFERENSI/LAPAN/CLoud_masking/DEM/120065/*.TIF'
raster_list_dem=glob.glob(pathname)
filename_dem=[]
for b in [os.path.basename(z) for z in glob.glob(pathname)]:
    c=os.path.splitext(b)[0]
    filename_dem.append(c)

read2=[]
for d in raster_list_dem:
    band2=gdal.Open(d)
    read2.append(band2.GetRasterBand(1).ReadAsArray())
dem_dict= dict(zip(filename_dem, read2))

def raster_buffer(data, dist=60):
    row = data.shape[0]
    col = data.shape[1]
    in_array = data.astype(int)
    Xcell_size = 30
    Ycell_size = 30
    cell_size = (Xcell_size + Ycell_size) / 2
    cell_dist = dist / cell_size
    # in_array[in_array == (inband.GetNoDataValue() or 0 or -999)]=0
    out_array = np.zeros_like(in_array)
    temp_array = np.zeros_like(in_array)
    i, j, h, k = 0, 0, 0, 0
    print("Running distance buffer...")
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

print "Calculating Solar Position Angle..."
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
print "Calculating Reflectances..."
#Reflectance
reflectance_band1=(float(data['REFLECTANCE_MULT_BAND_1'])*my_dict[filename[0][:-2]+'B1']+float(data['REFLECTANCE_ADD_BAND_1']))/cos(zenit_angle)
reflectance_band2=(float(data['REFLECTANCE_MULT_BAND_2'])*my_dict[filename[0][:-2]+'B2']+float(data['REFLECTANCE_ADD_BAND_2']))/cos(zenit_angle)
reflectance_band3=(float(data['REFLECTANCE_MULT_BAND_3'])*my_dict[filename[0][:-2]+'B3']+float(data['REFLECTANCE_ADD_BAND_3']))/cos(zenit_angle)
reflectance_band4=(float(data['REFLECTANCE_MULT_BAND_4'])*my_dict[filename[0][:-2]+'B4']+float(data['REFLECTANCE_ADD_BAND_4']))/cos(zenit_angle)
reflectance_band5=(float(data['REFLECTANCE_MULT_BAND_5'])*my_dict[filename[0][:-2]+'B5']+float(data['REFLECTANCE_ADD_BAND_5']))/cos(zenit_angle)
reflectance_band6=(float(data['REFLECTANCE_MULT_BAND_6'])*my_dict[filename[0][:-2]+'B6']+float(data['REFLECTANCE_ADD_BAND_6']))/cos(zenit_angle)
reflectance_band7=(float(data['REFLECTANCE_MULT_BAND_7'])*my_dict[filename[0][:-2]+'B7']+float(data['REFLECTANCE_ADD_BAND_7']))/cos(zenit_angle)
reflectance_band9=(float(data['REFLECTANCE_MULT_BAND_9'])*my_dict[filename[0][:-2]+'B9']+float(data['REFLECTANCE_ADD_BAND_9']))/cos(zenit_angle)

#Brightness Temperature
radiance_tirs1=(float(data['RADIANCE_MULT_BAND_10'])*my_dict[filename[0][:-2]+'B10']+float(data['RADIANCE_ADD_BAND_10']))
bt_n=float(data['K2_CONSTANT_BAND_10'])/np.log(float(data['K1_CONSTANT_BAND_10'])/radiance_tirs1+1)
BT=bt_n-273.15

print "Calculating BASIC FMASK..."
#Equation 1
#NDVI
NDVI=numexpr.evaluate("(reflectance_band5 - reflectance_band4) / (reflectance_band5 + reflectance_band4)")
#NDSI
NDSI=numexpr.evaluate("(reflectance_band3-reflectance_band6)/(reflectance_band3+reflectance_band6)")

#NDVI[numexpr.evaluate("(reflectance_nir + reflectance_red) == 0")] = 0.01
#NDSI[numexpr.evaluate("(reflectance_green + reflectance_swir1) == 0")] = 0.01
basic_test=numexpr.evaluate( "(reflectance_band7 > 0.3)& (NDSI < 0.8) & (NDVI < 0.8) & (BT < 27) ")

#Equation 2
#Mean Vis
meanVIS = numexpr.evaluate("(reflectance_band2 + reflectance_band3 + reflectance_band4) / 3 ")
whiteness = numexpr.evaluate("(abs(reflectance_band2 - meanVIS) + abs(reflectance_band3 - meanVIS)+ abs(reflectance_band4 - meanVIS)) / meanVIS")
basic_test&= whiteness < 0.7
del meanVIS

#Equation 3
# Haze test
hazeTest = (numexpr.evaluate("reflectance_band2- 0.5 * reflectance_band4 - 0.08"))> 0  # Haze test
basic_test &= numexpr.evaluate("(hazeTest > 0)")
#del HOT  # need to find thick warm cloud

#Equation 4
basic_test &= numexpr.evaluate("(reflectance_band5 / reflectance_band6) > 0.75") # Equation 6
#bqa
bqa= my_dict[filename[0][:-2]+'BQA']
prob_cloud= numexpr.evaluate("(bqa == 53248)|(bqa == 56320)| (bqa == 61440)| (bqa== 64512) ") # https://landsat.usgs.gov/qualityband

print "Combining FMASK and Band BQA..."
#FMVsBQ
prob_cloud2= numexpr.evaluate("((prob_cloud==True)|(basic_test==True))")
print "Filtering and Cleaning..."
filtering = ndimage.median_filter(prob_cloud2, 3)
#Cleaning
cleaning = ndimage.binary_erosion(filtering)
#buffering
raster_buffer_array=raster_buffer(cleaning,60)
#stampday=datetime.date.today()
#filename_today= stampday.strftime("%Y%m%d")
output=data['DATE_ACQUIRED']
export_array(raster_buffer_array, "D:/FORESTS2020/TRAINING/Python/RESULT/CLOUDMASKING/C310118/IPB/Cloud-"+ output+".TIF")
export_array(prob_cloud, "D:/FORESTS2020/TRAINING/Python/RESULT/CLOUDMASKING/C310118/IPB/bqa"+ output+".TIF")
print("Finish: running time: %s seconds" % (time.time() - start_time))