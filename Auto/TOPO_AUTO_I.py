# -*- coding: utf-8 -*-
"""
Created on Thu Sep 07 10:06:36 2017

@author: hudjimartsu
"""
from numpy import *
import numexpr,sys
import datetime
import math
from datetime import datetime, date
import numpy as np
import glob
from osgeo import gdal
from scipy.stats import linregress
from dict import dict
import os
import time
import pandas as pd
start_time = time.time()
f = open('D:/FORESTS2020/TRAINING/PyQgis/DATA/Landsat8/NO_SEA/FINAL/GEO/LANDSAT/TRIAL/mtl.txt', 'r') #open file for reading
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
path='D:/FORESTS2020/TRAINING/PyQgis/DATA/Landsat8/NO_SEA/FINAL/GEO/LANDSAT/'
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
pathname='D:/FORESTS2020/TRAINING/PyQgis/DATA/Landsat8/NO_SEA/FINAL/GEO/DEM/*.tif'
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

#load cloud
cloud= gdal.Open("D:/FORESTS2020/TRAINING/PyQgis/RESULT/Landsat8/Cloud_Masking/C191217/Cloud-2016-08-17.TIF")
band3 = cloud.GetRasterBand(1)
cloudArray = band3.ReadAsArray()

# def raster_buffer(data, dist=60):
#     row = data.shape[0]
#     col = data.shape[1]
#     in_array = data.astype(int)
#     Xcell_size = 30
#     Ycell_size = 30
#     cell_size = (Xcell_size + Ycell_size) / 2
#     cell_dist = dist / cell_size
#     # in_array[in_array == (inband.GetNoDataValue() or 0 or -999)]=0
#     out_array = np.zeros_like(in_array)
#     temp_array = np.zeros_like(in_array)
#     i, j, h, k = 0, 0, 0, 0
#     print("Running distance buffer...")
#     while (h < col):
#         k = 0
#         while (k < row):
#             if (in_array[k][h] >= 1):
#                 i = h - cell_dist
#                 while ((i < cell_dist + h) and i < col):
#                     j = k - cell_dist
#                     while (j < (cell_dist + k) and j < row):
#                         if (((i - h) ** 2 + (j - k) ** 2) <= cell_dist ** 2):
#                             if (temp_array[j][i] == 0 or temp_array[j][i] > ((i - h) ** 2 + (j - k) ** 2)):
#                                 out_array[j][i] = in_array[k][h]
#                                 temp_array[j][i] = (i - h) ** 2 + (j - k) ** 2
#                         j += 1
#                     i += 1
#             k += 1
#         h += 1
#     d, temp_array, in_array = None, None, None
#     return out_array


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

#sun azimuth angle
theta1= -1 * ((sin(lat_n)) * cos(zenit_angle)- sin(decl_deg)/(cos (lat_n) * sin (zenit_angle)))
theta2=np.arccos(theta1) #radians
theta3=np.rad2deg(theta2)#degree
azimuth_angle=180 - theta3 #degrees

# IC calculation
delta=azimuth_angle - dem_dict['aspect2']
IC=(cos(zenit_angle)* cos (dem_dict['slope2'])) + (sin(zenit_angle) * sin (dem_dict['slope2']) * cos(delta))#radians
print "Calculating Reflectances..."
# #Reflectance
# # reflectance_blue=(float(data['REFLECTANCE_MULT_BAND_2'])*my_dict[filename[0][:-2]+'B2']+float(data['REFLECTANCE_ADD_BAND_2']))/cos(zenit_angle)
# reflectance_red=(float(data['REFLECTANCE_MULT_BAND_4'])*my_dict[filename[0][:-2]+'B4']+float(data['REFLECTANCE_ADD_BAND_4']))/cos(zenit_angle)
# reflectance_nir=(float(data['REFLECTANCE_MULT_BAND_5'])*my_dict[filename[0][:-2]+'B5']+float(data['REFLECTANCE_ADD_BAND_5']))/cos(zenit_angle)
# # reflectance_swir1=(float(data['REFLECTANCE_MULT_BAND_6'])*my_dict[filename[0][:-2]+'B6']+float(data['REFLECTANCE_ADD_BAND_6']))/cos(zenit_angle)
# # reflectance_swir2=(float(data['REFLECTANCE_MULT_BAND_7'])*my_dict[filename[0][:-2]+'B7']+float(data['REFLECTANCE_ADD_BAND_7']))/cos(zenit_angle)
# # reflectance_green=(float(data['REFLECTANCE_MULT_BAND_3'])*my_dict[filename[0][:-2]+'B3']+float(data['REFLECTANCE_ADD_BAND_3']))/cos(zenit_angle)
#
#
# NDVI=numexpr.evaluate("(reflectance_nir - reflectance_red) / (reflectance_nir + reflectance_red)")
#
# sample_ndvi_I= numexpr.evaluate("(NDVI >= 0) & (NDVI < 0.1)")
# sample_ndvi_II= numexpr.evaluate("(NDVI >= 0.1) & (NDVI <= 0.5)")
# sample_ndvi_III= numexpr.evaluate("(NDVI >0.5)")
# slope40 = dem_dict['slope2'] > 35
# masking_NDVI_1= numexpr.evaluate("(sample_ndvi_I == True) & (cloudArray == 0) & (slope40 == True)")
# masking_NDVI_2= numexpr.evaluate("(sample_ndvi_II == True) & (cloudArray == 0) & (slope40 == True)")
# masking_NDVI_3= numexpr.evaluate("(sample_ndvi_III == True) & (cloudArray == 0) & (slope40 == True)")
# masking_NDVI_3[masking_NDVI_1==True]=1
# masking_NDVI_3[masking_NDVI_1==False]=0
# sample
area_true= dem_dict['sample'].nonzero() #outputnya index row n col
a_true=area_true[0]
b_true=area_true[1]
#correction
cos_zenith= cos(zenit_angle)

#auto
my_dict_f=dict((k, my_dict[k]) for k in ('LC81220652016230LGN00_B1', 'LC81220652016230LGN00_B2', 'LC81220652016230LGN00_B3', 'LC81220652016230LGN00_B4','LC81220652016230LGN00_B5', 'LC81220652016230LGN00_B6', 'LC81220652016230LGN00_B7', 'LC81220652016230LGN00_B9'))
temp={}
IC_final={}
slope={}
for y in my_dict_f:
    val2=my_dict_f[y]
    temp[y]=val2[a_true,b_true].ravel()
    IC_true=IC[a_true,b_true].ravel()
    slope[y]=linregress(IC_true, temp[y])
    IC_final[y]=my_dict[y]-(slope[y][0]*(IC-cos_zenith))
# generate random boolean mask the length of data
# use p 0.75 for False and 0.25 for True
# mask_NDVI_1 = np.random.choice([False, True], shape(masking_NDVI_1), p=[0.75, 0.25])
# export_array(raster_buffer_array, "D:/FORESTS2020/TRAINING/PyQgis/RESULT/Landsat8/Cloud_Masking/C191217/CloudLAPANcmd-"+ output+".TIF")
print(IC_final)
#export data to raster (GeoTiff)
for item in IC_final:
    print(item)
    geo = band.GetGeoTransform()
    proj = band.GetProjection()
    shape = my_dict['LC81220652016230LGN00_B2'].shape
    driver = gdal.GetDriverByName("GTiff")
    dst_ds = driver.Create("D:/FORESTS2020/TRAINING/PyQgis/RESULT/Landsat8/TOPO/TP201217/"+item+"_fix.TIF", shape[1], shape[0], 1, gdal.GDT_Float32)
    dst_ds.SetGeoTransform(geo)
    dst_ds.SetProjection(proj)
    dst_ds.GetRasterBand(1).WriteArray(IC_final[item])
    dst_ds.FlushCache()
    dst_ds = None  # save, close"""


# export_array(sample_ndvi_I, "D:/FORESTS2020/TRAINING/PyQgis/RESULT/Landsat8/TOPO/TP201217/sample_ndvi_2_I.TIF")
# export_array(sample_ndvi_II, "D:/FORESTS2020/TRAINING/PyQgis/RESULT/Landsat8/TOPO/TP201217/sample_ndvi_2_I.TIF")
# export_array(sample_ndvi_III, "D:/FORESTS2020/TRAINING/PyQgis/RESULT/Landsat8/TOPO/TP201217/sample_ndvi_3_I.TIF")
# export_array(masking_NDVI_1, "D:/FORESTS2020/TRAINING/PyQgis/RESULT/Landsat8/TOPO/TP201217/mask_NDVI40_I.TIF")
# export_array(masking_NDVI_2, "D:/FORESTS2020/TRAINING/PyQgis/RESULT/Landsat8/TOPO/TP201217/mask_NDVI40_II.TIF")
# export_array(masking_NDVI_3, "D:/FORESTS2020/TRAINING/PyQgis/RESULT/Landsat8/TOPO/TP201217/mask_NDVI40_III.TIF")
print("Finish: running time: %s seconds" % (time.time() - start_time))

csv2 = {}
for c in IC_final:
        val4=IC_final[c]
        csv2[c]= val4[a_true, b_true].ravel()
        #print "hasil", csv2
        df=pd.DataFrame(csv2)
        df2=pd.DataFrame(temp)
        df3=pd.DataFrame({'IC':IC_true})
        dfn=pd.concat([df3, df, df2], axis=1)
        dfn.to_csv("D:/FORESTS2020/TRAINING/PyQgis/RESULT/Landsat8/TOPO/TP201217/sample2.csv", index= False)
