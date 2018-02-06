import pandas as pd
import numpy as np
import numexpr
from IPython.core.hooks import late_startup_hook
from osgeo import gdal
import gdal
import glob
import math
from affine import Affine
from scipy.stats import linregress
import csv
#Load data raster
raster_list2=glob.glob('D:\PROJECT\FOREST 2020\TRAINING\SENTINEL\Data\land8\*.tif')
read=[]
for i in raster_list2:
    band=gdal.Open(i)
    read.append(band.GetRasterBand(1).ReadAsArray())

aspect=np.array(read[0], dtype=float)
band2=np.array(read[1], dtype=float)
band5=np.array(read[2], dtype=float)
band6=np.array(read[3], dtype=float)
sample=np.array(read[4], dtype=float)
slope=np.array(read[5], dtype=float)

day_of_year=228-1
y=(2*math.pi/365)*(day_of_year +((10-12)/24))

#sun declination angle
sun_dec= (0.006918 - (0.399912*np.cos(y))) + (0.070257*np.sin(y)) - (0.006758*np.cos(2*y)) + (0.000907*np.sin(2*y))-(0.002697*np.cos(3*y)) + (0.00148*np.sin (3*y))

#lat long value
# get columns and rows of your image from gdalinfo
filename = "D:\PROJECT\FOREST 2020\TRAINING\SENTINEL\Data\land8\slope.tif" #path to raster
ds=gdal.Open(filename, gdal.GA_ReadOnly)
xoff, a, b, yoff, d, e = ds.GetGeoTransform()
def pixel2coord(x, y):
    """Returns global coordinates from pixel x, y coords"""
    xp = a * x + b * y + xoff
    yp = d * x + e * y + yoff
    return(xp, yp)

rows=aspect.shape[0]
colms=aspect.shape[1]
coba=[]
for row in  range(0,rows):
  for col in  range(0,colms):
  coba.append(pixel2coord(col,row))
coba_2=np.array(coba, dtype=float)
long=coba_2[:,0]
lat=coba_2[:,1]
long_n=long.reshape(rows,colms)
lat_n=lat.reshape(rows,colms)

#time offset
toff= (229.18*(0.000075 + 0.001868*np.cos(y))) - (0.032077*np.sin(y)) - (0.014615*np.cos(2*y)) - (0.040849*np.sin(2*y))- ((4*long_n)+ (60*7))

#sun hour angle
alpha= (((60*10)+0+ (34/60)+toff)/4)-180

#sun zenith angle
zenit_angle=np.arccos((np.sin(lat_n)*(np.sin(sun_dec)))+ (np.cos(lat_n)*(np.cos(sun_dec)*(np.cos(alpha)))))

#sun azimuth angle
theta= 180 + (np.arccos(-(((np.sin(lat_n))*(np.cos(zenit_angle)))- (np.sin(sun_dec))/(np.cos(lat_n))*(np.sin(zenit_angle)))))

# IC calculation
IC = np.cos(zenit_angle)* np.cos(slope) + np.sin(zenit_angle)*np.sin(slope)*(np.cos(theta-slope))


# sample
sample_1=sample.ravel()
area= np.where((sample_1>=1),1,0)
area2=area.reshape(rows,colms)
area_true= area2.nonzero() #outputnya index row n col
a_true=area_true[0]
b_true=area_true[1]
#correction
#calculate slope regression
band2_data=band2[a_true,b_true]
IC_data=IC[a_true,b_true]
slope2=linregress(band2_data, IC_data)
a=band2_data.ravel()
b=IC_data.ravel()
#IC_11= davel[2] - (slope11[0]*(IC-1))
IC_2= band2 - (slope2[0]*(IC-zenit_angle))

#band 6
band6_data=band6[a_true,b_true]
slope6=linregress(band6_data, IC_data)
IC_6= band6 - (slope6[0]*(IC-zenit_angle))

#band 5
band5_data=band5[a_true,b_true]
slope5=linregress(band5_data, IC_data)
IC_5= band5 - (slope5[0]*(IC-zenit_angle))

#export
geo = band.GetGeoTransform()
proj = band.GetProjection()
shape = band2.shape
driver = gdal.GetDriverByName("GTiff")
dst_ds = driver.Create("D:/PROJECT/FOREST 2020/TRAINING/SENTINEL/Result/IC_2a_net.tif", shape[1], shape[0], 1, gdal.GDT_Float32)
dst_ds.SetGeoTransform(geo)
dst_ds.SetProjection(proj)
dst_ds.GetRasterBand(1).WriteArray(IC_2)
dst_ds = None  # save, close"""