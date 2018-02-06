# import glob, os, shutil
# source_dir= 'D:/FORESTS2020/DATA/SRTM/SUMATERA/SUMATRA/'
# files = glob.iglob(os.path.join(source_dir, "*.tif"))
# dest_dir = 'D:/FORESTS2020/DATA/SRTM/SUMATERA/SUMATRA/TEST/'
# for file in files:
#     if os.path.isfile(file):
#         shutil.copy2(file, dest_dir)
# import datetime
import math
from datetime import datetime, date
import numpy as np
import glob
import os
import shutil
from osgeo import gdal
from scipy.stats import linregress
import pandas as pd
# from dict import dict
import numexpr
# os.chdir("D:/FORESTS2020/TRAINING/PyQgis/DATA/Landsat8/SUMSEL/GEO/FINAL/PATH ROW")
os.chdir("D:/FORESTS2020/DATA/SRTM/SUMATERA/SUMATRA/NEW/")
# Open metadata file
txt=[]
folder_name=[]
for filename in os.listdir(os.getcwd()):
    filename=os.path.join(os.getcwd(), filename)
    folder_n=os.path.basename(filename)
    folder_name.append(folder_n)
    metadata= glob.iglob(filename, '*.tif')
    txt.append(metadata)
    print metadata



RootDir1 = r'D:/FORESTS2020/DATA/SRTM/PAPUA/'
TargetFolder = r'D:/FORESTS2020/DATA/SRTM/PAPUA/PAPUA_DEM/'
for root, dirs, files in os.walk((os.path.normpath(RootDir1)), topdown=False):
        for name in files:
            if name.endswith('.tif'):
                print "Found"
                SourceFolder = os.path.join(root,name)
                shutil.copy2(SourceFolder, TargetFolder) #copies tif to new folderRootDir1 = r'D:/FORESTS2020/DATA/SRTM/SUMATERA/SUMATRA/NEW/'


RootDir1 = r'D:/FORESTS2020/DATA/SRTM/SUMATERA/SUMATRA/NEW/'
TargetFolder = r'D:/FORESTS2020/DATA/SRTM/SUMATERA/SUMATRA/TEST/'
for root, dirs, files in os.walk((os.path.normpath(RootDir1)), topdown=False):
        for name in files:
            if name.endswith('.tif'):
                print "Found"
                SourceFolder = os.path.join(root,name)
                shutil.copy2(SourceFolder, TargetFolder) #copies tif to new folderRootDir1 = r'D:/FORESTS2020/DATA/SRTM/SUMATERA/SUMATRA/NEW/'


RootDir1 = r'D:/FORESTS2020/DATA/SRTM/NUSA TENGGARA/Tiles_201208200632/'
TargetFolder = r'D:/FORESTS2020/DATA/SRTM/NUSA TENGGARA/NT/'
for root, dirs, files in os.walk((os.path.normpath(RootDir1)), topdown=False):
        for name in files:
            if name.endswith('.tif'):
                print "Found"
                SourceFolder = os.path.join(root,name)
                shutil.copy2(SourceFolder, TargetFolder) #copies tif to new folder


RootDir1 = r'D:/FORESTS2020/DATA/SRTM/SULAWESI/New folder/'
TargetFolder = r'D:/FORESTS2020/DATA/SRTM/SULAWESI/SW/'
for root, dirs, files in os.walk((os.path.normpath(RootDir1)), topdown=False):
        for name in files:
            if name.endswith('.tif'):
                print "Found"
                SourceFolder = os.path.join(root,name)
                shutil.copy2(SourceFolder, TargetFolder) #copies tif to new folder


###TEST READ SUBFOLDER
RootDir1 = r'D:/FORESTS2020/DATA/LANDSAT/CONTOH/DATA'
dest= r'D:/FORESTS2020/DATA/LANDSAT/CONTOH/HASIL/'
# TargetFolder = r'D:/FORESTS2020/DATA/SRTM/SUMATERA/SUMATRA/TEST/'
for root, dirs, files in os.walk((os.path.normpath(RootDir1)), topdown=False):
        for name in dirs:
            test = os.path.join(dest, name)
            if not os.path.exists(test):
                os.makedirs(test)
            else:
                shutil.rmtree(test)
                os.makedirs(test)
        filename = []
        for a in [os.path.basename(x) for x in glob.glob(os.path.join(root, name) +'\*.TIF')]:
            p = os.path.splitext(a)[0]
            filename.append(p)
        raster = glob.glob(os.path.join(root, name) + '\*.TIF')
        read=[]
        for y in raster:
            band=gdal.Open(y)
            read.append(band.GetRasterBand(1).ReadAsArray().astype(float))





            # if name.endswith('.tif'):
            #     print "Found"
            #     SourceFolder = os.path.join(root,name)
            #     shutil.copy2(SourceFolder, TargetFolder) #copies tif to new folder




dest= "D:/FORESTS2020/DATA/SRTM/SUMATERA/SUMATRA/TEST/"
src_files = os.listdir(os.getcwd())
for file_name in src_files:
    full_file_name = os.path.join(src_files, file_name)
    if (os.path.isfile(full_file_name)):
        shutil.copy(full_file_name, dest)



# # export folder
# # make folder
# tes="D:/FORESTS2020/DATA/LANDSAT/OLAH/DEM"
# for i in folder_name:
#     test= os.path.join(tes, i)
#     if not os.path.exists(test):
#         os.makedirs(test)
#     else:
#         shutil.rmtree(test)
#         os.makedirs(test)
#         print test
# for y in txt:
#     shutil.copy(y, os.listdir(tes))
#
# tes2="D:/FORESTS2020/DATA/LANDSAT/OLAH/SAMPLE/122062"
# shutil.copyfile(txt[0], tes2)
#
#
#
# f=open(glob_f[1])
# def build_data(f):
#     output = {}
#     for line in f.readlines():
#         if "=" in line:
#             l = line.split("=")
#             output[l[0].strip()] = l[1].strip()
#     return output
# data = build_data(f)