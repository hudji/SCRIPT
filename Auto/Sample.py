import datetime
import math
from datetime import datetime, date
import numpy as np
import glob
import os
from osgeo import gdal
from scipy.stats import linregress
import time
from dict import dict
start_time = time.time()
f = open('D:/FORESTS2020/TRAINING/PyQgis/DATA/Landsat8/NO_SEA/FINAL/GEO/LANDSAT/TRIAL//mtl.txt', 'r') #open file for reading
def build_data(f): #build dictionarysss
    output = {} #Dict
    for line in f.readlines(): #Iterates through every line in the string
        if "=" in line: #make sure line has data as wanted
            l = line.split("=") #Seperate by "=" and put into a list
            output[l[0].strip()] = l[1].strip() #First word is key, second word is value
    return output #Returns a dictionary with the key, value pairs.
data = build_data(f)
print "Loading Data Raster..."
#Load data raster
path='D:/FORESTS2020/TRAINING/PyQgis/DATA/Landsat8/NO_SEA/FINAL/GEO/LANDSAT/TRIAL/'
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
folder=filename[0][3:9]
pathname='D:/FORESTS2020/TRAINING/PyQgis/DATA/Landsat8/NO_SEA/FINAL/GEO/DEM2/'+folder+'/*.tif'
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

#load sample
cloud= gdal.Open("D:/FORESTS2020/TRAINING/PyQgis/RESULT/Landsat8/Cloud_Masking/C191217/Cloud-2016-08-17.TIF")
band3 = cloud.GetRasterBand(1)
cloudArray = band3.ReadAsArray()

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

# sample



area_true= dem_dict['sample'].nonzero() #outputnya index row n col
a_true=area_true[0]
b_true=area_true[1]

#correction
cos_zenith= cos(zenit_angle)

#auto
temp={}
IC_final={}
slope={}
for y in my_dict:
    val2=my_dict[y]
    temp[y]=val2[a_true,b_true].ravel()
    IC_true=IC[a_true,b_true].ravel()
    slope[y]=linregress(IC_true, temp[y])
    IC_final[y]=my_dict[y]-(slope[y][0]*(IC-cos_zenith))


print(IC_final)
#export data to raster (GeoTiff)
for item in IC_final:
    print(item)
    geo = band.GetGeoTransform()
    proj = band.GetProjection()
    shape = my_dict['LC81220652016230LGN00_B2'].shape
    driver = gdal.GetDriverByName("GTiff")
    dst_ds = driver.Create("D:/FORESTS2020/TRAINING/PyQgis/RESULT/Landsat8/TOPO/TP141217/"+item+"_fix.TIF", shape[1], shape[0], 1, gdal.GDT_Float32)
    dst_ds.SetGeoTransform(geo)
    dst_ds.SetProjection(proj)
    dst_ds.GetRasterBand(1).WriteArray(IC_final[item])
    dst_ds.FlushCache()
    dst_ds = None  # save, close"""
print("Finish: running time: %s seconds" % (time.time() - start_time))

export_array(sample_slope2, "D:/FORESTS2020/TRAINING/PyQgis/RESULT/Landsat8/TOPO/TP201217/slope40.TIF")