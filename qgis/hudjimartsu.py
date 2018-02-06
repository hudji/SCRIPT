# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 11:11:11 2017

@author: sahid
"""

from osgeo import gdal
import numpy as np
from numpy import *
import pandas as pd
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
bqa=np.array(read[11], dtype=float) #bandbqa

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
d = {'znir':ravel_nir,'zswir1':ravel_swir1, 'zswir2':ravel_swir2}
dataframe = pd.DataFrame(d)
E=dataframe.query('znir > 3 or zswir1 > 3 or zswir2 > 3')
Hasil = np.repeat(0.0,len(ravel_swir1))
Hasil[E.index] = 1
training_sample = Hasil.reshape(row,col)