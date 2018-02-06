import numpy as np
from osgeo import gdal
import glob
import numexpr
from datetime import datetime, date
import math
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

raster_list=glob.glob('D:\PROJECT\FOREST 2020\TRAINING\PyQgis\DATA\Landsat8\cloud\*.tif')
read=[]
for i in raster_list:
    band=gdal.Open(i)
    read.append(band.GetRasterBand(1).ReadAsArray())


band2=np.array(read[0], dtype=float)
band4=np.array(read[1], dtype=float)
band5=np.array(read[2], dtype=float)
band6=np.array(read[3], dtype=float)
band7=np.array(read[4], dtype=float)

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
     + 0.000907 * sin (2 * gamma) - 0.002697 * cos (3 * gamma) + 0.00148 * sin (2 * gamma) #radians
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


#Reflectance
reflectance_blue=(float(data['REFLECTANCE_MULT_BAND_2'])*band2+float(data['REFLECTANCE_ADD_BAND_2']))/cos(zenit_angle)
reflectance_red=(float(data['REFLECTANCE_MULT_BAND_4'])*band4+float(data['REFLECTANCE_ADD_BAND_4']))/cos(zenit_angle)
reflectance_nir=(float(data['REFLECTANCE_MULT_BAND_5'])*band5+float(data['REFLECTANCE_ADD_BAND_5']))/cos(zenit_angle)
reflectance_swir1=(float(data['REFLECTANCE_MULT_BAND_6'])*band6+float(data['REFLECTANCE_ADD_BAND_6']))/cos(zenit_angle)
reflectance_swir2=(float(data['REFLECTANCE_MULT_BAND_7'])*band7+float(data['REFLECTANCE_ADD_BAND_7']))/cos(zenit_angle)

#calculate z
znir=(reflectance_nir-reflectance_nir.mean())/np.std(reflectance_nir)
zswir1=(reflectance_swir1-reflectance_swir1.mean())/np.std(reflectance_swir1)
zswir2=(reflectance_swir2-reflectance_swir2.mean())/np.std(reflectance_swir2)

# Sample z > 3
row=znir.shape[0]
col=znir.shape[1]
ravel_nir=znir.ravel()
ravel_swir1=zswir1.ravel()
ravel_swir2=zswir2.ravel()
test=np.column_stack((ravel_nir, ravel_swir1, ravel_swir2))
test2= np.where ((test[:,0] > 3) | (test [:,1] > 3) | (test[:,2] > 3),1,0) # training_sample
test3= np.where ((test[:,0] > 3) | (test [:,1] > 3) | (test[:,2] > 3),0,1)
training_sample_true = test2.reshape(row,col)
training_sample_false = test3.reshape(row,col)
coba_true= training_sample_true.nonzero() #outputnya index row n col
coba_false= training_sample_false.nonzero()
a_true=coba_true[0]
b_true=coba_true[1]
a_false=coba_false[0]
b_false=coba_false[1]

#NDVI
ndvi=numexpr.evaluate("(reflectance_nir - reflectance_red)/(reflectance_nir + reflectance_red)")
#NDBI
ndbi=numexpr.evaluate("(reflectance_swir1-reflectance_nir)/(reflectance_swir1+reflectance_nir)")

#BSI
bsi=numexpr.evaluate("((reflectance_swir1+reflectance_nir)-(reflectance_swir1+reflectance_blue))/((reflectance_swir1+reflectance_nir)+(reflectance_swir1+reflectance_blue))")

#calculate NIR
Matriks_true = np.zeros(band5.shape)
Matriks_true[a_true,b_true]= reflectance_nir[a_true,b_true]
Matriks_false = np.zeros(band5.shape)
Matriks_false[a_false,b_false]= reflectance_nir[a_false,b_false]
mean_true_nir=np.nanmean(np.where(Matriks_true!=0,Matriks_true,np.nan))
std_true_nir= np.nanstd(np.where(Matriks_true!=0,Matriks_true,np.nan))
mean_false_nir=np.nanmean(np.where(Matriks_false!=0,Matriks_false,np.nan))
std_false_nir= np.nanstd(np.where(Matriks_false!=0,Matriks_false,np.nan))
del Matriks_true, Matriks_false
#calculate Swir1
Matriks_true = np.zeros(band6.shape)
Matriks_true[a_true,b_true]= reflectance_swir1[a_true,b_true]
Matriks_false = np.zeros(band6.shape)
Matriks_false[a_false,b_false]= reflectance_swir1[a_false,b_false]
mean_true_swir1=np.nanmean(np.where(Matriks_true!=0,Matriks_true,np.nan))
std_true_swir1= np.nanstd(np.where(Matriks_true!=0,Matriks_true,np.nan))
mean_false_swir1=np.nanmean(np.where(Matriks_false!=0,Matriks_false,np.nan))
std_false_swir1= np.nanstd(np.where(Matriks_false!=0,Matriks_false,np.nan))
del Matriks_true, Matriks_false
#calculate swir2
Matriks_true = np.zeros(band7.shape)
Matriks_true[a_true,b_true]= reflectance_swir2[a_true,b_true]
Matriks_false = np.zeros(band7.shape)
Matriks_false[a_false,b_false]= reflectance_swir2[a_false,b_false]
mean_true_swir2=np.nanmean(np.where(Matriks_true!=0,Matriks_true,np.nan))
std_true_swir2= np.nanstd(np.where(Matriks_true!=0,Matriks_true,np.nan))
mean_false_swir2=np.nanmean(np.where(Matriks_false!=0,Matriks_false,np.nan))
std_false_swir2= np.nanstd(np.where(Matriks_false!=0,Matriks_false,np.nan))
del Matriks_true, Matriks_false
#calculate ndvi
Matriks_true = np.zeros(ndvi.shape)
Matriks_true[a_true,b_true]= ndvi[a_true,b_true]
Matriks_false = np.zeros(ndvi.shape)
Matriks_false[a_false,b_false]= ndvi[a_false,b_false]
mean_true_ndvi=np.nanmean(np.where(Matriks_true!=0,Matriks_true,np.nan))
std_true_ndvi= np.nanstd(np.where(Matriks_true!=0,Matriks_true,np.nan))
mean_false_ndvi=np.nanmean(np.where(Matriks_false!=0,Matriks_false,np.nan))
std_false_ndvi= np.nanstd(np.where(Matriks_false!=0,Matriks_false,np.nan))
del Matriks_true, Matriks_false

##classification
#Rule 1
rule1a=np.where((reflectance_nir >= (min(mean_false_nir, mean_true_nir))+(max(mean_false_nir, mean_true_nir)- min(mean_false_nir, mean_true_nir))/2), 1, 0)
rule1b=np.where((reflectance_swir1 >= (min(mean_false_swir1, mean_true_swir1))+(max(mean_false_swir1, mean_true_swir1)- min(mean_false_swir1, mean_true_swir1))/2), 1, 0)
rule1c=np.where((reflectance_swir2 >= (min(mean_false_swir2, mean_true_swir2))+(max(mean_false_swir2, mean_true_swir2)- min(mean_false_swir2, mean_true_swir2))/2), 1, 0)
rule_ndvi=np.where((ndvi < (min(mean_false_ndvi, mean_true_ndvi))+(max(mean_false_ndvi, mean_true_ndvi)- min(mean_false_ndvi, mean_true_ndvi))/2), 1, 0)
rule_bsi = np.where(( bsi > 0), 1, 0)
rule_ndbi= np.where((ndbi > 0), 1, 0)
ravel_rule1a=rule1a.ravel()
ravel_rule1b=rule1b.ravel()
ravel_rule1c=rule1c.ravel()
ravel_rule_ndvi=rule_ndvi.ravel()
ravel_rule_bsi=rule_bsi.ravel()
ravel_rule_ndbi=rule_ndbi.ravel()
ravel_rule1=np.column_stack((ravel_rule1a, ravel_rule1b, ravel_rule1c, ravel_rule_ndvi, ravel_rule_bsi, ravel_rule_ndbi))
ravel_rule1n= np.where ((ravel_rule1[:,0] == 1) | (ravel_rule1 [:,1] == 1) | (ravel_rule1[:,2] == 1) & (ravel_rule1[:,3] == 1) & (ravel_rule1[:,4] == 1) & (ravel_rule1[:,5] == 1) ,1,0)
class_rule1= ravel_rule1n.reshape(row, col)

#Rule 2
rule2a=np.where((reflectance_nir >= min((mean_false_nir+std_false_nir), (mean_true_nir-std_true_nir))+ (max((mean_false_nir+std_false_nir), (mean_true_nir-std_true_nir))- min((mean_false_nir+std_false_nir), (mean_true_nir-std_true_nir)))/2), 1, 0)
rule2b=np.where((reflectance_swir1 >= min((mean_false_swir1+std_false_swir1), (mean_true_swir1-std_true_swir1))+(max((mean_false_swir1+std_false_swir1), (mean_true_swir1-std_true_swir1))- min((mean_false_swir1+std_false_swir1), (mean_true_swir1-std_true_swir1)))/2), 1, 0)
rule2c=np.where((reflectance_swir2 >= min((mean_false_swir2+std_false_swir2), (mean_true_swir2-std_true_swir2))+(max((mean_false_swir2+std_false_swir2), (mean_true_swir2-std_true_swir2))- min((mean_false_swir2+std_false_swir2), (mean_true_swir2-std_true_swir2)))/2), 1, 0)
rule_ndvi=np.where((ndvi < (min(mean_false_ndvi, mean_true_ndvi))+(max(mean_false_ndvi, mean_true_ndvi)- min(mean_false_ndvi, mean_true_ndvi))/2), 1, 0)
ravel_rule2a=rule2a.ravel()
ravel_rule2b=rule2b.ravel()
ravel_rule2c=rule2c.ravel()
ravel_rule2=np.column_stack((ravel_rule2a, ravel_rule2b, ravel_rule2c, ravel_rule_ndvi, ravel_rule_bsi, ravel_rule_ndbi))
ravel_rule2n= np.where ((ravel_rule2[:,0] == 1) | (ravel_rule2[:,1] == 1) | (ravel_rule2[:,2] == 1) & (ravel_rule2[:,3] == 1) & (ravel_rule2[:,4] == 1) & (ravel_rule2[:,5] == 1) ,1,0)
class_rule2= ravel_rule2n.reshape(row, col)

#Rule 3
rule3a=np.where((reflectance_nir >= min((mean_false_nir+2*std_false_nir), (mean_true_nir-2*std_true_nir))+ (max((mean_false_nir+2*std_false_nir), (mean_true_nir-2*std_true_nir))- min((mean_false_nir+2*std_false_nir), (mean_true_nir-2*std_true_nir)))/2), 1, 0)
rule3b=np.where((reflectance_swir1 >= min((mean_false_swir1+2*std_false_swir1), (mean_true_swir1-2*std_true_swir1))+(max((mean_false_swir1+2*std_false_swir1), (mean_true_swir1-2*std_true_swir1))- min((mean_false_swir1+2*std_false_swir1), (mean_true_swir1-2*std_true_swir1)))/2), 1, 0)
rule3c=np.where((reflectance_swir2 >= min((mean_false_swir2+2*std_false_swir2), (mean_true_swir2-2*std_true_swir2))+(max((mean_false_swir2+2*std_false_swir2), (mean_true_swir2-2*std_true_swir2))- min((mean_false_swir2+2*std_false_swir2), (mean_true_swir2-2*std_true_swir2)))/2), 1, 0)
rule_ndvi=np.where((ndvi < (min(mean_false_ndvi, mean_true_ndvi))+(max(mean_false_ndvi, mean_true_ndvi)- min(mean_false_ndvi, mean_true_ndvi))/2), 1, 0)
ravel_rule3a=rule3a.ravel()
ravel_rule3b=rule3b.ravel()
ravel_rule3c=rule3c.ravel()
ravel_rule3=np.column_stack((ravel_rule3a, ravel_rule3b, ravel_rule3c, ravel_rule_ndvi, ravel_rule_bsi, ravel_rule_ndbi))
ravel_rule3n= np.where((ravel_rule3[:,0] == 1) | (ravel_rule3[:,1] == 1) | (ravel_rule3[:,2] == 1) & (ravel_rule3[:,3] == 1) & (ravel_rule3[:,4] == 1) & (ravel_rule3[:,5] == 1) ,1,0)
class_rule3= ravel_rule3n.reshape(row, col)

#Rule 4
rule4a=np.where((reflectance_nir >= min((mean_false_nir+3*std_false_nir), (mean_true_nir-3*std_true_nir))+ (max((mean_false_nir+3*std_false_nir), (mean_true_nir-3*std_true_nir))- min((mean_false_nir+3*std_false_nir), (mean_true_nir-3*std_true_nir)))/2), 1, 0)
rule4b=np.where((reflectance_swir1 >= min((mean_false_swir1+3*std_false_swir1), (mean_true_swir1-3*std_true_swir1))+(max((mean_false_swir1+3*std_false_swir1), (mean_true_swir1-3*std_true_swir1))- min((mean_false_swir1+3*std_false_swir1), (mean_true_swir1-3*std_true_swir1)))/2), 1, 0)
rule4c=np.where((reflectance_swir2 >= min((mean_false_swir2+3*std_false_swir2), (mean_true_swir2-3*std_true_swir2))+(max((mean_false_swir2+3*std_false_swir2), (mean_true_swir2-3*std_true_swir2))- min((mean_false_swir2+3*std_false_swir2), (mean_true_swir2-3*std_true_swir2)))/2), 1, 0)
rule_ndvi=np.where((ndvi < (min(mean_false_ndvi, mean_true_ndvi))+(max(mean_false_ndvi, mean_true_ndvi)- min(mean_false_ndvi, mean_true_ndvi))/2), 1, 0)
ravel_rule4a=rule4a.ravel()
ravel_rule4b=rule4b.ravel()
ravel_rule4c=rule4c.ravel()
ravel_rule4=np.column_stack((ravel_rule4a, ravel_rule4b, ravel_rule4c, ravel_rule_ndvi, ravel_rule_bsi, ravel_rule_ndbi))
ravel_rule4n= np.where ((ravel_rule4[:,0] == 1) | (ravel_rule4 [:,1] == 1) | (ravel_rule4[:,2] == 1) & (ravel_rule4[:,3] == 1) & (ravel_rule4[:,4] == 1) & (ravel_rule4[:,5] == 1) ,1,0)
class_rule4= ravel_rule4n.reshape(row, col)

#export
geo = band.GetGeoTransform()
proj = band.GetProjection()
shape = band2.shape
driver = gdal.GetDriverByName("GTiff")
dst_ds = driver.Create("D:\PROJECT\FOREST 2020\TRAINING\PyQgis\RESULT\Landsat8\Cloud Masking\C031117\Rule3.tif", shape[1], shape[0], 1, gdal.GDT_Float32)
dst_ds.SetGeoTransform(geo)
dst_ds.SetProjection(proj)
dst_ds.GetRasterBand(1).WriteArray(class_rule3)
dst_ds = None  # save, close"""

