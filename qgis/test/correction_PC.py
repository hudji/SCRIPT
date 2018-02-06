import pandas as pd
import numpy as np
import numexpr
from osgeo import gdal
import glob
import math
from scipy.stats import linregress
import csv
#Load data raster
raster_list=glob.glob('D:\PROJECT\FOREST 2020\TRAINING\SENTINEL\Data\olah\*.tif')
read=[]
for i in raster_list:
    band=gdal.Open(i)
    read.append(band.GetRasterBand(1).ReadAsArray())

aspect=np.array(read[0], dtype=float)
azimuth=np.array(read[1],dtype=float)
band11=np.array(read[2],dtype=float)
band2=np.array(read[3],dtype=float)
band8=np.array(read[4],dtype=float)
sample= read [5]
slope=np.array(read[6],dtype=float)
zenith=np.array(read[7],dtype=float)
data=[aspect, azimuth, band11, band2, band8,sample, slope, zenith]
davel=[]
for a in data:
    ravel=a.ravel()
    davel.append(ravel)
row = band8.shape[0]
col = band8.shape[1]
## Illumination Condition
IC= (np.cos(davel[7])* np.cos(davel[6])) + (np.sin(davel[7])* np.sin(davel[6])*np.cos(davel[1]-davel[0]))
IC_fix= IC.reshape(row,col)

# sample
area= np.where((davel[5] >=1),1,0)
area2=area.reshape(row,col)
area_true= area2.nonzero() #outputnya index row n col
a_true=area_true[0]
b_true=area_true[1]

#calculate slope regression
band11_data=band11[a_true,b_true]
IC_fix_data=IC_fix[a_true,b_true]
slope11=linregress(band11_data, IC_fix_data)
a=band11_data.ravel()
b=IC_fix_data.ravel()
#IC_11= davel[2] - (slope11[0]*(IC-1))
IC_11= davel[2] - (0.0784*(IC-1))
IC_11F=IC_11.reshape(row,col)

band8_data=band8[a_true, b_true]
slope8=linregress(band8_data, IC_fix_data)
IC_8= davel[4] - (slope8[0]*(IC-1))
IC_8F=IC_8.reshape(row,col)
a8=band8_data.ravel()

band2_data=band2[a_true, b_true]
slope2=linregress(band2_data, IC_fix_data)
IC_2= davel[3] - (slope2[0]*(IC-1))
IC_2F=IC_2.reshape(row,col)
a2=band2_data.ravel()

##Export
df=pd.DataFrame({'band11':a, 'IC':b})
df.to_csv('D:/PROJECT/FOREST 2020/TRAINING/SENTINEL/Result/slope.csv', index=False)


geo = band.GetGeoTransform()
proj = band.GetProjection()
shape = band8.shape
driver = gdal.GetDriverByName("GTiff")
dst_ds = driver.Create("D:/PROJECT/FOREST 2020/TRAINING/SENTINEL/Result/IC2.tif", shape[1], shape[0], 1, gdal.GDT_Float32)
dst_ds.SetGeoTransform(geo)
dst_ds.SetProjection(proj)
dst_ds.GetRasterBand(1).WriteArray(IC_2F)
dst_ds = None  # save, close"""