import datetime
import math
from datetime import datetime, date
import numpy as np
import glob
from osgeo import gdal

#Load data raster (Landsat 8 or Sentinel 2)
raster_list=glob.glob('Data\*.tif')
read=[]
for i in raster_list:
    band=gdal.Open(i)
    read.append(band.GetRasterBand(1).ReadAsArray())

#Load Metadata file for landsat 8
f = open('metadata file') #open file for reading
def build_data(f):
    output = {}
    for line in f.readlines():
        if "=" in line:
            l = line.split("=")
            output[l[0].strip()] = l[1].strip()
    return output
data = build_data(f)

def year_date():
    date_file=data['FILE_DATE']
    yearTahun=date_file[:10]
    time_data=date_file[-8:-1]
    all= yearTahun+" "+time_data
    dt = datetime.strptime(all, '%Y-%m-%d %H:%M:%S')
    return dt

def hour():
    if dt.hour <= 6:
        print dt.hour + 12
    else:
        print int(dt.hour)
def leap():
    if (dt.year % 4) == 0:
        if (dt.year % 100) == 0:
            if (dt.year % 400) == 0:
               print '{:d}'.format(366)
            else:
                print '{:d}'.format(365)
        else:
            print '{:d}'.format(366)
    else:
        print '{:d}'.format(365)
def cos(x):
    cos= np.cos(np.deg2rad(x))
    return  cos
def sin(x):
    sin=np.sin(np.deg2rad(x))
    return sin
def dayOfYear():
    day_date= date(dt.year, dt.month, dt.day)
    sum_of_day=int(day_date.strftime('%j'))
    return sum_of_day

gamma=((2 * math.pi) / leap()) * ((dayOfYear() - 1) + ((hour() - 12) / 24) )

#sun declination angle
decl=0.006918 - 0.399912 * cos(gamma) + 0.070257 * sin(gamma) - 0.006758 * cos (2 * gamma)\
     + 0.000907 * sin (2 * gamma) - 0.002697 * cos (3 * gamma) + 0.00148 * sin (2 * gamma) #radians
decl_deg= (360 / (2 * math.pi)) * decl

# get lat long value of your image each pixel
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
tst=10 * 60 + 0 + 34.019582 / 60 + timeoff #minutes
ha=(tst /4)-180 #degree

#sun zenith angle
zenit1 =sin(lat_n)* sin(decl_deg) + cos (lat_n)* cos(decl_deg) * cos(ha)
zenit2=np.arccos(zenit1) #radians
zenit_angle= np.rad2deg(zenit2)

#sun azimuth angle
theta1= -1 * (sin(lat_n)) * cos(zenit_angle)- sin(decl_deg)/(cos (lat_n) * sin (zenit_angle))
theta2=np.arccos(theta1) #radians
theta3=np.rad2deg(theta2)#degree
azimuth_angle=180 - theta3 #degrees

# IC calculation
delta=azimuth_angle - aspect
IC=(cos(zenit_angle)* cos (slope)) + (sin(zenit_angle) * sin (slope) * cos(delta))#radians

# sample
area_true= read[4].nonzero() #sample area
a_true=area_true[0]
b_true=area_true[1]

#correction band
cos_zenith= cos(zenit_angle)
band2_data=read[1][a_true,b_true]
IC_data=IC[a_true,b_true]
a=band2_data.ravel()
b=IC_data.ravel()
slope2=linregress(b,a)
IC_2= band2 -(slope2[0]*(IC-cos_zenith))

#export as Raster
geo = band.GetGeoTransform()
proj = band.GetProjection()
shape = band2.shape
driver = gdal.GetDriverByName("GTiff")
dst_ds = driver.Create("Ouput folder\IC_2.tif", shape[1], shape[0], 1, gdal.GDT_Float32)
dst_ds.SetGeoTransform(geo)
dst_ds.SetProjection(proj)
dst_ds.GetRasterBand(1).WriteArray(IC_2)
dst_ds = None  # save, close"""
