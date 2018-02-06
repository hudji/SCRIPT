import datetime
import math
from datetime import datetime, date
import numpy as np
import glob
import os
from osgeo import gdal
from scipy.stats import linregress
import pandas as pd
# from dict import dict
import numexpr
#Load Metadata
path_f = 'D:/FORESTS2020/DATA/LANDSAT/PATH ROW/125061/*.txt'#open file for reading
glob_f= glob.glob(path_f)
f=open(glob_f[1])
def build_data(f):
    output = {}
    for line in f.readlines():
        if "=" in line:
            l = line.split("=")
            output[l[0].strip()] = l[1].strip()
    return output
data = build_data(f)

#Load data raster
print "Loading Data Raster..."
#Load data raster
path='D:/FORESTS2020/DATA/LANDSAT/PATH ROW/125061/'
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

#Load data raster aspect & slope
pathname='D:/FORESTS2020/DATA/LANDSAT/DEM/'
raster_list_dem=glob.glob(pathname+filename[0][10:16]+'/*.TIF')
filename_dem=[]
for b in [os.path.basename(z) for z in glob.glob(pathname+filename[0][10:16]+'/*.TIF')]:
    c=os.path.splitext(b)[0]
    filename_dem.append(c)

read2=[]
for d in raster_list_dem:
    band2=gdal.Open(d)
    read2.append(band2.GetRasterBand(1).ReadAsArray())
dem_dict= dict(zip(filename_dem, read2))

#Load data raster sample area based on Baplan
pathsample='D:/FORESTS2020/DATA/LANDSAT/SAMPLE/'
sampledata=gdal.Open(pathsample +filename[0][10:16]+'/baplan.TIF')
sample=np.array(sampledata.GetRasterBand(1).ReadAsArray())

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
print "Calculating Solar Position..."
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
delta=azimuth_angle - dem_dict['aspect']
IC=(cos(zenit_angle)* cos (dem_dict['slope'])) + (sin(zenit_angle) * sin (dem_dict['slope']) * cos(delta))#radians

print "Calculating Reflectances..."
#Reflectance
reflectance_band1=(float(data['REFLECTANCE_MULT_BAND_1'])*my_dict[filename[0][:-2]+'B1']+float(data['REFLECTANCE_ADD_BAND_1']))/cos(zenit_angle)
reflectance_band2=(float(data['REFLECTANCE_MULT_BAND_2'])*my_dict[filename[0][:-2]+'B2']+float(data['REFLECTANCE_ADD_BAND_2']))/cos(zenit_angle)
reflectance_band3=(float(data['REFLECTANCE_MULT_BAND_3'])*my_dict[filename[0][:-2]+'B3']+float(data['REFLECTANCE_ADD_BAND_3']))/cos(zenit_angle)
reflectance_band4=(float(data['REFLECTANCE_MULT_BAND_4'])*my_dict[filename[0][:-2]+'B4']+float(data['REFLECTANCE_ADD_BAND_4']))/cos(zenit_angle)
reflectance_band5=(float(data['REFLECTANCE_MULT_BAND_5'])*my_dict[filename[0][:-2]+'B5']+float(data['REFLECTANCE_ADD_BAND_5']))/cos(zenit_angle)
reflectance_band6=(float(data['REFLECTANCE_MULT_BAND_6'])*my_dict[filename[0][:-2]+'B6']+float(data['REFLECTANCE_ADD_BAND_6']))/cos(zenit_angle)
reflectance_band7=(float(data['REFLECTANCE_MULT_BAND_7'])*my_dict[filename[0][:-2]+'B7']+float(data['REFLECTANCE_ADD_BAND_7']))/cos(zenit_angle)
reflectance_band9=(float(data['REFLECTANCE_MULT_BAND_9'])*my_dict[filename[0][:-2]+'B9']+float(data['REFLECTANCE_ADD_BAND_9']))/cos(zenit_angle)

reflectance_f= {filename[0][:-2]+'B1':reflectance_band1, filename[0][:-2]+'B2':reflectance_band2,filename[0][:-2]+'B3':reflectance_band3, filename[0][:-2]+'B4':reflectance_band4, filename[0][:-2]+'B5':reflectance_band5, filename[0][:-2]+'B6':reflectance_band6, filename[0][:-2]+'B7':reflectance_band7, filename[0][:-2]+'B9':reflectance_band9}

print "Sampling..."
# sample
NDVI=numexpr.evaluate("(reflectance_band5 - reflectance_band4) / (reflectance_band5 + reflectance_band4)")
hutan= sample == 1
sample_ndvi= numexpr.evaluate("(NDVI >0.5) & (hutan==True)")
# sample_ndvi= NDVI > 0.5
area_true= sample_ndvi.nonzero() #outputnya index row n col
a_true=area_true[0]
b_true=area_true[1]

#correction
cos_zenith= cos(zenit_angle)

#auto
#def IC_all(my_dict):
temp={}
IC_final={}
for y in reflectance_f:
        val2=reflectance_f[y]
        temp[y]=val2[a_true,b_true].ravel()
        IC_true=IC[a_true,b_true].ravel()
        slope=linregress(IC_true, temp[y])
        IC_final[y]=reflectance_f[y]-(slope[0]*(IC-cos_zenith))
print "Exporting to GeoTIFF..."
#export auto
for item in IC_final:
    geo = band.GetGeoTransform()
    proj = band.GetProjection()
    shape = my_dict[filename[0][:-2]+'B1'].shape
    driver = gdal.GetDriverByName("GTiff")
    dst_ds = driver.Create("D:/FORESTS2020/TRAINING/Python/RESULT/TOPO/TP240118/125061/" +item + "topo.tif", shape[1], shape[0], 1, gdal.GDT_Float64)
    dst_ds.SetGeoTransform(geo)
    dst_ds.SetProjection(proj)
    ds=dst_ds.GetRasterBand(1)
    ds.SetNoDataValue(-9999)
    ds.WriteArray(IC_final[item])
    # outband.SetNoDataValue(-9999)
    dst_ds.FlushCache()
    dst_ds = None  # save, close""


def export_array(in_array, output_path):
    """This function is used to produce output of array as a map."""
    global proj, geotrans, row, col
    proj     = band.GetProjection()
    geotrans = band.GetGeoTransform()
    row      = band.RasterYSize
    col      = band.RasterXSize
    driver   = gdal.GetDriverByName("GTiff")
    outdata  = driver.Create(output_path, col, row, 1, gdal.GDT_Float32)
    outband  = outdata.GetRasterBand(1)
    outband.SetNoDataValue(-9999)
    outband.WriteArray(in_array)
    # Georeference the image
    outdata.SetGeoTransform(geotrans)
    # Write projection information
    outdata.SetProjection(proj)
    outdata.FlushCache()
    outdata = None

export_array(sample_ndvi, "D:/FORESTS2020/TRAINING/Python/RESULT/TOPO/TP220118/126062/sample.TIF")

csv2 = {}
for c in IC_final:
        val4=IC_final[c]
        csv2[c]= val4[a_true,b_true].ravel()
        #print "hasil", csv2
        df=pd.DataFrame(csv2)
        df2=pd.DataFrame(temp)
        df3=pd.DataFrame({'IC':IC_true})
        dfn=pd.concat([df3, df, df2], axis=1)
        dfn.to_csv("D:/FORESTS2020/TRAINING/PyQgis/RESULT/Landsat8/TOPO/TP020118/sample12363.csv", index= False)