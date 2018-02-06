import datetime
import math
from datetime import datetime, date
import numpy as np
import glob
from osgeo import gdal
from scipy.stats import linregress
import pandas as pd
import os
from dict import dict
#Load Metadata
f = open('D:/FORESTS2020/TRAINING/PyQgis/DATA/Landsat8/clip/mtl.txt', 'r') #open file for reading
def build_data(f):
    output = {}
    for line in f.readlines():
        if "=" in line:
            l = line.split("=")
            output[l[0].strip()] = l[1].strip()
    return output
data = build_data(f)

#Load data raster
# raster_list2=glob.glob('D:/FORESTS2020/TRAINING/PyQgis/DATA/Landsat8/clip/*.tif')
path='D:/FORESTS2020/TRAINING/PyQgis/DATA/Landsat8/clip/'
raster_list=glob.glob(path+ '*.tif')
read=[]
for i in raster_list:
    band=gdal.Open(i)
    read.append(band.GetRasterBand(1).ReadAsArray())
filename=[]
for a in [os.path.basename(x) for x in glob.glob(path + '*.tif')]:
    p=os.path.splitext(a)[0]
    filename.append(p)
my_dict= dict(zip(filename, read))
#
# aspect=np.array(read[0], dtype=float)
# band2=np.array(read[1], dtype=float)
# band5=np.array(read[2], dtype=float)
# band6=np.array(read[3], dtype=float)
# sample=np.array(read[4], dtype=float)
# slope=np.array(read[5], dtype=float)

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

#sun azimuth angle
theta1= -1 * ((sin(lat_n)) * cos(zenit_angle)- sin(decl_deg)/(cos (lat_n) * sin (zenit_angle)))
theta2=np.arccos(theta1) #radians
theta3=np.rad2deg(theta2)#degree
azimuth_angle=180 - theta3 #degrees

# IC calculation
delta=azimuth_angle - my_dict['aspect']
IC=(cos(zenit_angle)* cos (my_dict['slope'])) + (sin(zenit_angle) * sin (my_dict['slope']) * cos(delta))#radians

# sample
area_true= read[4].nonzero() #outputnya index row n col
a_true=area_true[0]
b_true=area_true[1]

#correction
cos_zenith= cos(zenit_angle)

#Reflectance
reflectance_band2=(float(data['REFLECTANCE_MULT_BAND_2'])*my_dict['band2']+float(data['REFLECTANCE_ADD_BAND_2']))/cos(zenit_angle)
reflectance_band5=(float(data['REFLECTANCE_MULT_BAND_5'])*my_dict['band5']+float(data['REFLECTANCE_ADD_BAND_5']))/cos(zenit_angle)
reflectance_band6=(float(data['REFLECTANCE_MULT_BAND_6'])*my_dict['band6']+float(data['REFLECTANCE_ADD_BAND_6']))/cos(zenit_angle)

reflectance_f= {'reflectance2':reflectance_band2, 'reflectance5':reflectance_band5, 'reflectance6':reflectance_band6}
#
# #band2
# band2_data=reflectance_band2[a_true,b_true]
# IC_data=IC[a_true,b_true]
# a=band2_data.ravel()
# b=IC_data.ravel()
# slope2=linregress(b,a)
# IC_2= reflectance_band2 -(slope2[0]*(IC-cos_zenith))
#
# #band 6
# band6_data=reflectance_band6[a_true,b_true]
# band6_R=band6_data.ravel()
# slope6=linregress(b, band6_R)
# IC_6= reflectance_band6 - (slope6[0]*(IC-cos_zenith))
#
# #band 5
# band5_data=reflectance_band5[a_true,b_true]
# band5_R=band5_data.ravel()
# slope5=linregress(b, band5_R)
# IC_5= reflectance_band5 - (slope5[0]*(IC-cos_zenith))

# AUTO
# # my_dict_f=dict((k, my_dict[k]) for k in ('band2', 'band5', 'band6'))
temp={}
IC_final={}
slope={}
for y in reflectance_f:
    val2=reflectance_f[y]
    temp[y]=val2[a_true,b_true].ravel()
    IC_true=IC[a_true,b_true].ravel()
    slope[y]=linregress(IC_true, temp[y])
    IC_final[y]=reflectance_f[y]-(slope[y][0]*(IC-cos_zenith))
# my_dict_f=dict((k, my_dict[k]) for k in ('band2', 'band5', 'band6'))
# temp={}
# IC_final={}
# slope={}
# for y in my_dict_f:
#     val2=my_dict_f[y]
#     temp[y]=val2[a_true,b_true].ravel()
#     IC_true=IC[a_true,b_true].ravel()
#     slope[y]=linregress(IC_true, temp[y])
#     IC_final[y]=my_dict_f[y]-(slope[y][0]*(IC-cos_zenith))

#export
geo = band.GetGeoTransform()
proj = band.GetProjection()
shape = my_dict['band2'].shape
driver = gdal.GetDriverByName("GTiff")
dst_ds = driver.Create("D:/FORESTS2020/TRAINING/PyQgis/RESULT/Landsat8/L011117/Band2R.tif", shape[1], shape[0], 1, gdal.GDT_Float32)
dst_ds.SetGeoTransform(geo)
dst_ds.SetProjection(proj)
dst_ds.GetRasterBand(1).WriteArray(reflectance_band2)
dst_ds = None  # save, close"""

#CSV
# #corrected
# IC_2D=IC_2[a_true,b_true]
# IC_2R=IC_2D.ravel()
# IC_5D=IC_5[a_true,b_true]
# IC_5R=IC_5D.ravel()
# IC_6D=IC_6[a_true,b_true]
# IC_6R=IC_6D.ravel()


df=pd.DataFrame({'IC':b,'band2un':a,'band2corr':IC_2R,'band5un':band5_R,'band5corr':IC_5R,'band6un':band6_R,'band6corr':IC_6R })
df.to_csv('D:/FORESTS2020/TRAINING/PyQgis/RESULT/Landsat8/TOPO/TP211217/sample.csv', index=False)

csv2 = {}
for c in IC_final:
        val4=IC_final[c]
        csv2[c]= val4[a_true,b_true].ravel()
        #print "hasil", csv2
        df=pd.DataFrame(csv2)
        df2=pd.DataFrame(temp)
        df3=pd.DataFrame({'IC':IC_true})
        dfn=pd.concat([df3, df, df2], axis=1)
        dfn.to_csv("D:/FORESTS2020/TRAINING/PyQgis/RESULT/Landsat8/TOPO/TP211217/sample_V.csv", index= False)