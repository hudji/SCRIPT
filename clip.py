# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 13:24:33 2017

@author: sahid
"""
import os, fnmatch
CLIP='D:/PROJECT/FOREST 2020/TRAINING/PyQgis/DATA/source/clip.shp'
INPUT_FOLDER='D:/PROJECT/FOREST 2020/TRAINING/PyQgis/DATA/source/'
OUTPUT_FOLDER= 'D:/PROJECT/FOREST 2020/TRAINING/PyQgis/DATA/clip/'

def findRasters (path, filter):
    for root, dirs, files in os.walk(path):
        for file in fnmatch.filter(files, filter):
            yield file

for raster in findRasters(INPUT_FOLDER, '*.tif'):
    inRaster = INPUT_FOLDER + '/' + raster
    outRaster = OUTPUT_FOLDER + '/clip_' + raster
    cmd = 'gdalwarp -q -cutline %s -crop_to_cutline %s %s' % (CLIP, inRaster, outRaster)
    os.system(cmd)


"""
import os, fnmatch
from subprocess import call
#call(["ls", "-l"])

inFolder= 'D:/PROJECT/FOREST 2020/TRAINING/PyQgis/DATA/source/'
outFolder= 'D:/PROJECT/FOREST 2020/TRAINING/PyQgis/DATA/clip/'

os.chdir (inFolder)

def findRasters (path, filter):
    for root, dirs, files in os.walk(path, filter):
        for file in fnmatch.filter(files, filter):
            yield os.path.join (root, file)

for raster in findRasters (inFolder, '*.tif'):
    (infilepath, infilename)= os.path.split (raster)
    print infilename
    outRaster= outFolder+ 'clip_'+ infilename
    print outRaster
    warp= 'gdalwarp -dstnodata 0 -q -cutline %s -crop_to_cutline -of GTiff %s %s' % ('clip.shp', raster, outRaster)
    os.system(warp)
    #call (warp)"""