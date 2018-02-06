import pandas as pd
import numpy as np
import numexpr
from osgeo import gdal
import glob
import math
from scipy.stats import linregress
import csv

#Load data raster
raster_list2=glob.glob("D:\PROJECT\FOREST 2020\TRAINING\PyQgis\RAHMI's\data\New\*.tif")
read=[]
for i in raster_list2:
    band=gdal.Open(i)
    read.append(band.GetRasterBand(1).ReadAsArray())

aspect=np.array(read[0], dtype=float)
band10=np.array(read[1], dtype=float)
band2=np.array(read[2], dtype=float)
band3=np.array(read[3], dtype=float)
band4=np.array(read[4], dtype=float)
band5=np.array(read[5], dtype=float)
band6=np.array(read[6], dtype=float)
band7=np.array(read[7], dtype=float)
sample=np.array(read[8], dtype=float)
slope=np.array(read[9], dtype=float)

day_of_year=186
y=((2 * math.pi) / 365) * ((day_of_year - 1) + ((10.0952136 - 12) / 24) ) #radians
gamma=np.radians(y) #radians

#sun declination angle
decl=0.006918 - 0.399912 * np.cos(gamma) + 0.070257 * np.sin(gamma) - 0.006758 * np.cos (2 * gamma)\
     + 0.000907 * np.sin (2 * gamma) - 0.002697 * np.cos (3 * gamma) + 0.00148 * np.sin (3 * gamma) #radians
decl_deg= (360 / (2 * math.pi)) * decl

#lat long value
# get columns and rows of your image from gdalinfo
#filename = "D:\PROJECT\FOREST 2020\TRAINING\PyQgis\DATA\Landsat8\clip\slope.tif" #path to raster
#ds=gdal.Open(filename, gdal.GA_ReadOnly)
xoff, a, b, yoff, d, e = band.GetGeoTransform()
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

#eqtime
eqtime = 229.18 * (0.000075 + 0.001868 * np.cos(gamma) - 0.032077 * np.sin(gamma) - 0.014615 * np.cos(2 * gamma) - 0.040849 * np.sin(2 * gamma))  # minutes
timeoff= eqtime - 4 * long_n + 60 * 7 #minutes
tst=10 * 60 + 5 + 42.769024 / 60 + timeoff #minutes
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


# sample
sample_1=sample.ravel()
area= np.where((sample_1==1),1,0)
area2=area.reshape(rows,colms)
area_true= area2.nonzero() #outputnya index row n col
a_true=area_true[0]
b_true=area_true[1]

#correction
#calculate slope regression
cos_zenith= (np.cos(np.deg2rad(zenit_angle)))

#Reflectance
reflectance_band2=(float(data['REFLECTANCE_MULT_BAND_2'])*band2+float(data['REFLECTANCE_ADD_BAND_2']))/cos(zenit_angle)
reflectance_band5=(float(data['REFLECTANCE_MULT_BAND_5'])*band5+float(data['REFLECTANCE_ADD_BAND_5']))/cos(zenit_angle)
reflectance_band6=(float(data['REFLECTANCE_MULT_BAND_6'])*band6+float(data['REFLECTANCE_ADD_BAND_6']))/cos(zenit_angle)


#band 10
band10_data=band10[a_true,b_true]
IC_data=IC[a_true,b_true]
b=IC_data.ravel()
band10_R=band10_data.ravel()
slope10=linregress(b, band10_R)
IC_10= band10 - (slope10[0]*(IC-cos_zenith))

#band 2
band2_data=band2[a_true,b_true]
IC_data=IC[a_true,b_true]
band2_R=band2_data.ravel()
slope2=linregress(b,band2_R)
IC_2= band2 -(slope2[0]*(IC-cos_zenith))

#band 3
band3_data=band3[a_true,b_true]
IC_data=IC[a_true,b_true]
band3_R=band3_data.ravel()
slope3=linregress(b,band3_R)
IC_3= band3 -(slope3[0]*(IC-cos_zenith))

#band 4
band4_data=band4[a_true,b_true]
band4_R=band4_data.ravel()
slope4=linregress(b, band4_R)
IC_4= band4 - (slope4[0]*(IC-cos_zenith))

#band 5
band5_data=band5[a_true,b_true]
band5_R=band5_data.ravel()
slope5=linregress(b, band5_R)
IC_5= band5 - (slope5[0]*(IC-cos_zenith))

#band 6
band6_data=band6[a_true,b_true]
band6_R=band6_data.ravel()
slope6=linregress(b, band6_R)
IC_6= band6 - (slope6[0]*(IC-cos_zenith))

#band 7
band7_data=band7[a_true,b_true]
band7_R=band7_data.ravel()
slope7=linregress(b, band7_R)
IC_7= band7 - (slope7[0]*(IC-cos_zenith))


#export
geo = band.GetGeoTransform()
proj = band.GetProjection()
shape = band2.shape
driver = gdal.GetDriverByName("GTiff")
dst_ds = driver.Create("D:\PROJECT\FOREST 2020\TRAINING\PyQgis\RAHMI's\Result\R191017\C2\IC_5.tif", shape[1], shape[0], 1, gdal.GDT_Float32)
dst_ds.SetGeoTransform(geo)
dst_ds.SetProjection(proj)
dst_ds.GetRasterBand(1).WriteArray(IC_5)
dst_ds = None  # save, close"""

#CSV
#corrected
IC_10D=IC_10[a_true,b_true]
IC_10R=IC_10D.ravel()
IC_2D=IC_2[a_true,b_true]
IC_2R=IC_2D.ravel()
IC_3D=IC_3[a_true,b_true]
IC_3R=IC_3D.ravel()
IC_4D=IC_4[a_true,b_true]
IC_4R=IC_4D.ravel()
IC_5D=IC_5[a_true,b_true]
IC_5R=IC_5D.ravel()
IC_6D=IC_6[a_true,b_true]
IC_6R=IC_6D.ravel()
IC_7D=IC_7[a_true,b_true]
IC_7R=IC_7D.ravel()


df=pd.DataFrame({'IC':b,'band2un':band2_R,'band2corr':IC_2R,'band3un':band3_R,'band3corr':IC_3R, 'band4un':band4_R,'band4corr':IC_4R,'band5un':band5_R,'band5corr':IC_5R,'band6un':band6_R,'band6corr':IC_6R, 'band7un':band7_R,'band7corr':IC_7R })
df.to_csv("D:\PROJECT\FOREST 2020\TRAINING\PyQgis\RAHMI's\Result\R191017\C2\samplePoint2.csv", index=False)