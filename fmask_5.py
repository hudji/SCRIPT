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
from datetime import datetime, date
import scipy.stats
import scipy.signal
import scipy.ndimage.morphology
#from __future__ import division

#Load Metadata
f = open('D:\PROJECT\FOREST 2020\TRAINING\PyQgis\DATA\Landsat8\mtl.txt', 'r') #open file for reading
def build_data(f): #build dictionary
    output = {} #Dict
    for line in f.readlines(): #Iterates through every line in the string
        if "=" in line: #make sure line has data as wanted
            l = line.split("=") #Seperate by "=" and put into a list
            output[l[0].strip()] = l[1].strip() #First word is key, second word is value
    return output #Returns a dictionary with the key, value pairs.
data = build_data(f)

#Load data raster
raster_list=glob.glob('D:\PROJECT\FOREST 2020\TRAINING\PyQgis\DATA\Landsat8\GEO\*.tif')
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
cirrus=np.array(read[10],dtype=float)
tirs1=np.array(read[1],dtype=float)
tirs2=np.array(read[2],dtype=float)

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
     + 0.000907 * sin (2 * gamma) - 0.002697 * cos (3 * gamma) + 0.00148 * sin (2 * gamma) #radians
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


#Reflectance
reflectance_blue=(float(data['REFLECTANCE_MULT_BAND_2'])*blue+float(data['REFLECTANCE_ADD_BAND_2']))/cos(zenit_angle)
reflectance_red=(float(data['REFLECTANCE_MULT_BAND_4'])*red+float(data['REFLECTANCE_ADD_BAND_4']))/cos(zenit_angle)
reflectance_nir=(float(data['REFLECTANCE_MULT_BAND_5'])*nir+float(data['REFLECTANCE_ADD_BAND_5']))/cos(zenit_angle)
reflectance_swir1=(float(data['REFLECTANCE_MULT_BAND_6'])*swir1+float(data['REFLECTANCE_ADD_BAND_6']))/cos(zenit_angle)
reflectance_swir2=(float(data['REFLECTANCE_MULT_BAND_7'])*swir2+float(data['REFLECTANCE_ADD_BAND_7']))/cos(zenit_angle)
reflectance_green=(float(data['REFLECTANCE_MULT_BAND_3'])*green+float(data['REFLECTANCE_ADD_BAND_3']))/cos(zenit_angle)
reflectance_cirrus=(float(data['REFLECTANCE_MULT_BAND_9'])*cirrus+float(data['REFLECTANCE_ADD_BAND_9']))/cos(zenit_angle)

#Brightness Temperature
radiance_tirs1=(float(data['RADIANCE_MULT_BAND_10'])*tirs1+float(data['RADIANCE_ADD_BAND_10']))
bt_n=float(data['K2_CONSTANT_BAND_10'])/np.log(float(data['K1_CONSTANT_BAND_10'])/radiance_tirs1+1)
BT=bt_n-273.15
###Step I
#Eq. 1
#NDVI
NDVI=numexpr.evaluate("(reflectance_nir - reflectance_red) / (reflectance_nir + reflectance_red)")
#NDSI
NDSI=numexpr.evaluate("(reflectance_green-reflectance_swir1)/(reflectance_green+reflectance_swir1)")

#NDVI[numexpr.evaluate("(reflectance_nir + reflectance_red) == 0")] = 0.01
#NDSI[numexpr.evaluate("(reflectance_green + reflectance_swir1) == 0")] = 0.01
basic_test=numexpr.evaluate( "(reflectance_swir2 > 0.3)& (NDSI < 0.8) & (NDVI < 0.8) & (BT < 27) ")

#Mean Vis
meanVIS = numexpr.evaluate("(reflectance_blue + reflectance_green + reflectance_red) / 3 ")
whiteness = numexpr.evaluate("(abs(reflectance_blue - meanVIS) + abs(reflectance_green - meanVIS)+ abs(reflectance_red - meanVIS)) / meanVIS")
basic_test&= whiteness < 0.7
del meanVIS

# Haze test
HOT = numexpr.evaluate("reflectance_blue- 0.5 * reflectance_red - 0.08")  # Haze test
basic_test &= numexpr.evaluate("(HOT > 0)")
del HOT  # need to find thick warm cloud

# Ratio4/5>0.75 cloud test
basic_test &= numexpr.evaluate("(reflectance_nir / reflectance_swir1) > 0.75") # Equation 6
#Thin_prob
Thin_prob=reflectance_cirrus/0.04

# Cirrus tests from Landsat 8
basic_test2=basic_test
basic_test2|= numexpr.evaluate("Thin_prob > 0.25")
#water test
WT = numexpr.evaluate("((NDVI < 0.01) & (reflectance_nir < 0.11)) | ((NDVI < 0.1) & (NDVI > 0) & (reflectance_nir < 0.05))") # Equation 5


#clear sky water
#clear_sky_water2= numexpr.evaluate("(WT==True) & (reflectance_swir2<0.03)")
clear_sky_water3= WT[numexpr.evaluate("(WT==True) & (reflectance_swir2<0.03)")] #Equation 7
clear_sky_water= 1*clear_sky_water3
T_water= scipy.stats.scoreatpercentile(clear_sky_water,82.5) #Equation 8
wTemp_prob= (T_water-BT)/4 # Equation 9
wTemp_prob[numexpr.evaluate('wTemp_prob < 0')] = 0
wTemp_prob[numexpr.evaluate('wTemp_prob > 0')] = 1

# Brightness test (over water)
t_bright = 0.11
Brightness_prob = reflectance_nir / t_bright
Brightness_prob[Brightness_prob > 1] = 1
Brightness_prob[Brightness_prob < 0] = 0

#Water test
#wCloud_Prob=numexpr.evaluate('wTemp_prob * Brightness_prob + 100*Thin_prob')
wCloud_Prob=numexpr.evaluate('wTemp_prob * Brightness_prob + Thin_prob') # Equation 11
water_test=numexpr.evaluate('wCloud_Prob > 0.5')
#Temperature probability for land
Clear_sky_land2= numexpr.evaluate('(basic_test== False) & (WT==False)')
Clear_sky_land= 1*Clear_sky_land2
#Temperature
#Clear_sky_land2[numexpr.evaluate('(basic_test== False) & (water_test==False)')]
#Clear_sky_land= 1*Clear_sky_land2
T_low=scipy.stats.scoreatpercentile(Clear_sky_land2,17.5)
T_high= scipy.stats.scoreatpercentile(Clear_sky_land2,82.5)   #Equation 13


# ITemperature_Prob
ITemperature_Prob= (T_high+4-BT)/(T_high+4-(T_low - 4)) # Equation 14
ITemperature_Prob[ITemperature_Prob < 0] = 0
#Variability prob
#satu_Bv = numexpr.evaluate("(satu_B1 | satu_B2 | satu_B3)")
Vari_prob = 1 - np.maximum(np.maximum(np.absolute(NDSI), np.absolute(NDVI)), whiteness) # Equation 15

# cloud over land probability
ICloud_Prob =ITemperature_Prob * Vari_prob  #Equation 16


# Land threshold
land= ICloud_Prob[Clear_sky_land==0]
Land_threshold = scipy.stats.scoreatpercentile(land, 82.5) + 0.2# dynamic threshold (land)  # Equation 17

# Potential cloud
Potential_cloud= numexpr.evaluate('(basic_test2 & (wCloud_Prob > 0.5) | & (WT == 0)) | (basic_test2 & (wCloud_Prob > Land_threshold) & (WT == 1)) | (Temp < t_templ - 3500)')
Potential_cloud= numexpr.evaluate('(basic_test2==True) & (WT==True)&(wCloud_Prob > 0.5) | (basic_test2==True) & (WT==False) & (ICloud_Prob > Land_threshold) | (ICloud_Prob >0.99) & (WT==False)|(BT<T_low -35)')




geo = band.GetGeoTransform()
proj = band.GetProjection()
shape = red.shape
driver = gdal.GetDriverByName("GTiff")
dst_ds = driver.Create("D:\PROJECT\FOREST 2020\TRAINING\PyQgis\RESULT\Landsat8\Cloud Masking\C091117\IPotential_cloud.tif", shape[1], shape[0], 1, gdal.GDT_Float32)
dst_ds.SetGeoTransform(geo)
dst_ds.SetProjection(proj)
dst_ds.GetRasterBand(1).WriteArray(Potential_cloud)
dst_ds = None  # save, close"""

#print 'nilai', band
#print 'angka', ndvi

#CSV
#corrected
ndvi_c=NDVI.ravel()
band2=reflectance_blue.ravel()


df=pd.DataFrame({'NDVI':ndvi_c,'band2':band2})
df.to_csv('D:\PROJECT\FOREST 2020\TRAINING\PyQgis\RESULT\Landsat8\Cloud Masking\C091117\sampleR.csv', index=False)



