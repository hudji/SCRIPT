import pandas as pd
import numpy as np
from osgeo import gdal
import glob
import math
from scipy.stats import linregress
import csv

#Load data raster
raster_list2=glob.glob('D:\PROJECT\FOREST 2020\TRAINING\PyQgis\DATA\SENTINEL\Data\S91017\clip\*.tif')
read=[]
for i in raster_list2:
    band=gdal.Open(i)
    read.append(band.GetRasterBand(1).ReadAsArray())

aspect=np.array(read[0], dtype=float)
band11=np.array(read[1], dtype=float)
band2=np.array(read[2], dtype=float)
band8=np.array(read[3], dtype=float)
sample=np.array(read[4], dtype=float)
slope=np.array(read[5], dtype=float)

day_of_year=265
y=((2 * math.pi) / 366) * ((day_of_year - 1) + ((14 - 12) / 24) ) #radians
gamma=np.radians(y) #radians

#sun declination angle
decl=0.006918 - 0.399912 * np.cos(gamma) + 0.070257 * np.sin(gamma) - 0.006758 * np.cos (2 * gamma)\
     + 0.000907 * np.sin (2 * gamma) - 0.002697 * np.cos (3 * gamma) + 0.00148 * np.sin (2 * gamma) #radians
decl_deg= (360 / (2 * math.pi)) * decl

#lat long value
# get columns and rows of your image from gdalinfo
#filename = "D:\PROJECT\FOREST 2020\TRAINING\PyQgis\DATA\SENTINEL\Data\sentinel2\clip\slope.tif" #path to raster
#ds=gdal.Open(filename, gdal.GA_ReadOnly)
xoff, a, b, yoff, d, e = band.GetGeoTransform()
def pixel2coord(x, y):
    """Returns global coordinates from pixel x, y coords"""
    xp = a * x + b * y + xoff
    yp = d * x + e * y + yoff
    return(xp, yp)
rows=aspect.shape[0]
colms=aspect.shape[1]
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
eqtime = 229.18 * (0.000075 + 0.001868 * np.cos(gamma) - 0.032077 * np.sin(gamma) - 0.014615 * np.cos(2 * gamma) - 0.040849 * np.sin(2 * gamma))  # minutes
timeoff= eqtime - 4 * long_n + 60 * 7 #minutes
tst=14 * 60 + 55 + 41.026 / 60 + timeoff #minutes
ha=(tst /4)-180 #degree
#sun zenith angle
zenit1 =np.sin(np.deg2rad(lat_n))* np.sin(np.deg2rad(decl_deg)) + np.cos (np.deg2rad(lat_n))* np.cos(np.deg2rad(decl_deg)) * np.cos(np.deg2rad(ha))
zenit2=np.arccos(zenit1) #radians
zenit_angle= np.rad2deg(zenit2) #degrees
#sun azimuth angle
theta1= -1 * ((np.sin(np.deg2rad(lat_n)) * np.cos(np.deg2rad(zenit_angle))- np.sin(np.deg2rad(decl_deg)))/(np.cos (np.deg2rad(lat_n)) * np.sin (np.deg2rad(zenit_angle))))
theta2=np.arccos(theta1) #radians
theta3=np.rad2deg(theta2)#degree
azimuth_angle=180 - theta3 #degrees

# IC calculation
delta=azimuth_angle - aspect
IC=(np.cos(np.deg2rad(zenit_angle))* np.cos (np.deg2rad(slope))) + (np.sin(np.deg2rad(zenit_angle)) * np.sin (np.deg2rad(slope)) * np.cos(np.deg2rad(delta)))#radians

# sample (harus sama ukuran dengan pixel lainnya
area_true= sample.nonzero() #outputnya index row n col
a_true=area_true[0]
b_true=area_true[1]



#correction
#calculate slope regression
cos_zenith= np.cos(np.deg2rad(zenit_angle))
# Band 2
band2_data=band2[a_true,b_true]
IC_data=IC[a_true,b_true]
a=band2_data.ravel()
b=IC_data.ravel()
slope2=linregress(b,a)
IC_2= band2 -(slope2[0]*(IC-cos_zenith))

#band 11
band11_data=band11[a_true,b_true]
band11_R=band11_data.ravel()
slope11=linregress(b, band11_R)
IC_11= band11 - (slope11[0]*(IC-cos_zenith))

#band 8
band8_data=band8[a_true,b_true]
band8_R=band8_data.ravel()
slope8=linregress(b, band8_R)
IC_8= band8 - (slope8[0]*(IC-cos_zenith))

#export
geo = band.GetGeoTransform()
proj = band.GetProjection()
shape = band2.shape
driver = gdal.GetDriverByName("GTiff")
dst_ds = driver.Create("D:\PROJECT\FOREST 2020\TRAINING\PyQgis\RESULT\Sentinel2\S291017\IC_11.tif", shape[1], shape[0], 1, gdal.GDT_Float32)
dst_ds.SetGeoTransform(geo)
dst_ds.SetProjection(proj)
dst_ds.GetRasterBand(1).WriteArray(IC_11)
dst_ds = None  # save, close"""


#CSV
#Corrected
corr_2=IC_2[a_true,b_true]
corr_2N=corr_2.ravel()
corr_11=IC_11[a_true,b_true]
corr_11N=corr_11.ravel()
corr_8=IC_8[a_true,b_true]
corr_8N=corr_8.ravel()

df=pd.DataFrame({'IC':b,'band2un':a,'band2corr':corr_2N,'band8un':band8_R,'band8corr':corr_8N,'band11un':band11_R,'band11corr':corr_11N })
df.to_csv('D:\PROJECT\FOREST 2020\TRAINING\PyQgis\RESULT\Sentinel2\S291017\sample.csv', index=False)