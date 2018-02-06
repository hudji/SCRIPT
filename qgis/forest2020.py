#forest2020
from osgeo import gdal
import numpy as np
from numpy import *
import glob


#Load Metadata
f = open('D:/PROJECT/FOREST 2020/TRAINING/PyQgis/DATA/clip/mtl.txt', 'r') #open file for reading
def build_data(f): #build dictionary
    output = {} #Dict
    for line in f.readlines(): #Iterates through every line in the string
        if "=" in line: #make sure line has data as wanted
            l = line.split("=") #Seperate by "=" and put into a list
            output[l[0].strip()] = l[1].strip() #First word is key, second word is value
    return output #Returns a dictionary with the key, value pairs.
data = build_data(f)

#Load data raster
raster_list=glob.glob('D:/PROJECT/FOREST 2020/TRAINING/PyQgis/DATA/clip/*.TIF')
read=[]
for i in raster_list:
    band=gdal.Open(i)
    read.append(band.GetRasterBand(1).ReadAsArray())

blue=np.array(read[1], dtype=float) #band2
red=np.array(read[3],dtype=float) #band4
nir=np.array(read[4],dtype=float) #band5
swir1=np.array(read[5],dtype=float) #band6
swir2=np.array(read[6], dtype=float) #band7
bqa=np.array(read[11], dtype=float)#bqa

training= np.where((bqa > 45000), 1, 0)
#Reflectance
reflectance_blue=(float(data['REFLECTANCE_MULT_BAND_2'])*blue+float(data['REFLECTANCE_ADD_BAND_2']))/sin(54.68707440*3.1459/180)
reflectance_red=(float(data['REFLECTANCE_MULT_BAND_4'])*red+float(data['REFLECTANCE_ADD_BAND_4']))/sin(54.68707440*3.1459/180)
reflectance_nir=(float(data['REFLECTANCE_MULT_BAND_5'])*nir+float(data['REFLECTANCE_ADD_BAND_5']))/sin(54.68707440*3.1459/180)
reflectance_swir1=(float(data['REFLECTANCE_MULT_BAND_6'])*swir1+float(data['REFLECTANCE_ADD_BAND_6']))/sin(54.68707440*3.1459/180)
reflectance_swir2=(float(data['REFLECTANCE_MULT_BAND_7'])*swir2+float(data['REFLECTANCE_ADD_BAND_7']))/sin(54.68707440*3.1459/180)

#calculate z
znir=(reflectance_nir-reflectance_nir.mean())/np.std(reflectance_nir)
zswir1=(reflectance_swir1-reflectance_swir1.mean())/np.std(reflectance_swir1)
zswir2=(reflectance_swir2-reflectance_swir2.mean())/np.std(reflectance_swir2)

#Training Sample
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
upper=nir-red
lower=nir+red
mask = (upper == lower) & (upper==0)
ndvi=np.where(mask, 0, upper/lower)
del upper, lower, mask

#NDBI
upper=swir1-nir
lower=swir1+nir
mask = (upper == lower) & (upper==0)
ndbi=np.where(mask, 0, upper/lower)
del upper, lower, mask

#BSI
upper=(swir1+nir)-(swir1+blue)
lower= (swir1+nir)+(swir1+blue)
mask = (upper == lower) & (upper==0)
bsi=np.where(mask, 0, upper/lower)
del upper, lower, mask
#calculate NIR
Matriks_true = np.zeros(nir.shape)
Matriks_true[a_true,b_true]= reflectance_nir[a_true,b_true]
Matriks_false = np.zeros(nir.shape)
Matriks_false[a_false,b_false]= reflectance_nir[a_false,b_false]
mean_true_nir=np.nanmean(np.where(Matriks_true!=0,Matriks_true,np.nan))
std_true_nir= np.nanstd(np.where(Matriks_true!=0,Matriks_true,np.nan))
mean_false_nir=np.nanmean(np.where(Matriks_false!=0,Matriks_false,np.nan))
std_false_nir= np.nanstd(np.where(Matriks_false!=0,Matriks_false,np.nan))
del Matriks_true, Matriks_false
#calculate Swir1
Matriks_true = np.zeros(nir.shape)
Matriks_true[a_true,b_true]= reflectance_swir1[a_true,b_true]
Matriks_false = np.zeros(nir.shape)
Matriks_false[a_false,b_false]= reflectance_swir1[a_false,b_false]
mean_true_swir1=np.nanmean(np.where(Matriks_true!=0,Matriks_true,np.nan))
std_true_swir1= np.nanstd(np.where(Matriks_true!=0,Matriks_true,np.nan))
mean_false_swir1=np.nanmean(np.where(Matriks_false!=0,Matriks_false,np.nan))
std_false_swir1= np.nanstd(np.where(Matriks_false!=0,Matriks_false,np.nan))
del Matriks_true, Matriks_false
#calculate swir2
Matriks_true = np.zeros(nir.shape)
Matriks_true[a_true,b_true]= reflectance_swir2[a_true,b_true]
Matriks_false = np.zeros(nir.shape)
Matriks_false[a_false,b_false]= reflectance_swir2[a_false,b_false]
mean_true_swir2=np.nanmean(np.where(Matriks_true!=0,Matriks_true,np.nan))
std_true_swir2= np.nanstd(np.where(Matriks_true!=0,Matriks_true,np.nan))
mean_false_swir2=np.nanmean(np.where(Matriks_false!=0,Matriks_false,np.nan))
std_false_swir2= np.nanstd(np.where(Matriks_false!=0,Matriks_false,np.nan))
del Matriks_true, Matriks_false
#calculate ndvi
Matriks_true = np.zeros(nir.shape)
Matriks_true[a_true,b_true]= ndvi[a_true,b_true]
Matriks_false = np.zeros(nir.shape)
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
shape = red.shape
driver = gdal.GetDriverByName("GTiff")
dst_ds = driver.Create("D:/PROJECT/FOREST 2020/TRAINING/PyQgis/DATA/RESULT/REFLECTANCE/training_bqa.tif", shape[1], shape[0], 1, gdal.GDT_Float32)
dst_ds.SetGeoTransform(geo)
dst_ds.SetProjection(proj)
dst_ds.GetRasterBand(1).WriteArray(training)
dst_ds = None  # save, close"""
