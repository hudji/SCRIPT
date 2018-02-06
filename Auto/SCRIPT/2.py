import datetime
import math
from datetime import datetime, date
import numpy as np
import glob
import os
from osgeo import gdal
from scipy.stats import linregress
import pandas as pd
from dict import dict
#Load Metadata
f = open('D:/PROJECT/FOREST 2020/TRAINING/PyQgis/DATA/Landsat8/clip/mtl.txt', 'r') #open file for reading
def build_data(f):
    output = {}
    for line in f.readlines():
        if "=" in line:
            l = line.split("=")
            output[l[0].strip()] = l[1].strip()
    return output
data = build_data(f)

#Load data raster
raster_list=glob.glob('D:\PROJECT\FOREST 2020\TRAINING\PyQgis\DATA\Landsat8\NEW\*.tif')
read=[]
for i in raster_list:
    band=gdal.Open(i)
    read.append(band.GetRasterBand(1).ReadAsArray())
data_name=['band2', 'band5', 'band6']
my_dict= dict(zip(data_name, read))

raster_list_dem=glob.glob('D:\PROJECT\FOREST 2020\TRAINING\PyQgis\DATA\Landsat8\NEW\DEM\*.tif')
read2=[]
for d in raster_list_dem:
    band=gdal.Open(d)
    read2.append(band.GetRasterBand(1).ReadAsArray())

dem_name=['aspect', 'sample', 'slope']
dem_dict= dict (zip(dem_name, read2))

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
decl=0.006918 - 0.399912 * cos(gamma) + 0.070257 * sin(gamma) - 0.006758 * cos (2 * gamma) \
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

#sun azimuth angle
theta1= -1 * ((sin(lat_n)) * cos(zenit_angle)- sin(decl_deg)/(cos (lat_n) * sin (zenit_angle)))
theta2=np.arccos(theta1) #radians
theta3=np.rad2deg(theta2)#degree
azimuth_angle=180 - theta3 #degrees

# IC calculation
delta=azimuth_angle - dem_dict['aspect']
IC=(cos(zenit_angle)* cos (dem_dict['slope'])) + (sin(zenit_angle) * sin (dem_dict['slope']) * cos(delta))#radians

# sample
area_true= dem_dict['sample'].nonzero() #outputnya index row n col
a_true=area_true[0]
b_true=area_true[1]

#correction
cos_zenith= cos(zenit_angle)

#auto
#def IC_all(my_dict):
temp={}
IC_final={}
for y in my_dict:
    val2=my_dict[y]
    temp[y]=val2[a_true,b_true].ravel()
    IC_true=IC[a_true,b_true].ravel()
    slope=linregress(IC_true, temp[y])
    IC_final[y]=my_dict[y]-(slope[0]*(IC-cos_zenith))
    #return IC_final
#IC_new=IC_final(my_dict)
#export auto
geo = band.GetGeoTransform()
proj = band.GetProjection()
shape = my_dict['band2'].shape
driver = gdal.GetDriverByName("GTiff")
dst_ds = driver.Create("D:\PROJECT\FOREST 2020\TRAINING\PyQgis\DATA\SCRIPT\Auto\Result\Band5.tif", shape[1], shape[0], 1, gdal.GDT_Float32)
dst_ds.SetGeoTransform(geo)
dst_ds.SetProjection(proj)
dst_ds.GetRasterBand(1).WriteArray(IC_final['band5'])
dst_ds.FlushCache()
dst_ds = None  # save, close"""

print "Topo done"
#CSV
#corrected
"""
def csv():
    csv2 = {}
    for c in IC_final:
        val4=IC_final[c]
        csv2[c]= val4[a_true, b_true].ravel()
        #print "hasil", csv2
        df=pd.DataFrame(csv2)
        df2=pd.DataFrame(temp)
        df3=pd.DataFrame({'IC':IC_true})
        dfn=pd.concat([df3, df, df2], axis=1)
        dfn.to_csv("D:\PROJECT\FOREST 2020\TRAINING\PyQgis\RESULT\Landsat8\L301017\sample.csv", index= False)"""
