import datetime
import math
from datetime import datetime, date
import numpy as np
import glob
import os
import shutil
from osgeo import gdal
from scipy.stats import linregress
import pandas as pd
from dict import dict
import numexpr
# os.chdir("D:/FORESTS2020/TRAINING/PyQgis/DATA/Landsat8/SUMSEL/GEO/FINAL/PATH ROW")
os.chdir("D:/FORESTS2020/DATA/LANDSAT/OLAH/PATH ROW")
# Open metadata file
txt=[]
folder_name=[]
for filename in os.listdir(os.getcwd()):
    filename=os.path.join(os.getcwd(), filename)
    folder_n=os.path.basename(filename)
    folder_name.append(folder_n)
    metadata= glob.iglob(filename , '*.txt')
    txt.append(metadata)
    print metadata


#
# import glob, os, shutil
#
# files = glob.iglob(os.path.join(source_dir, "*.ext"))
# for file in files:
#     if os.path.isfile(file):
#         shutil.copy2(file, dest_dir)


# export folder
# make folder
tes="D:/FORESTS2020/DATA/LANDSAT/OLAH/DEM"
for i in folder_name:
    test= os.path.join(tes, i)
    if not os.path.exists(test):
        os.makedirs(test)
    else:
        shutil.rmtree(test)
        os.makedirs(test)
        print test
for y in txt:
    shutil.copy(y, os.listdir(tes))

tes2="D:/FORESTS2020/DATA/LANDSAT/OLAH/SAMPLE/122062"
shutil.copyfile(txt[0], tes2)



f=open(glob_f[1])
def build_data(f):
    output = {}
    for line in f.readlines():
        if "=" in line:
            l = line.split("=")
            output[l[0].strip()] = l[1].strip()
    return output
data = build_data(f)