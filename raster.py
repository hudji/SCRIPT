# -*- coding: utf-8 -*-
"""
Created on Thu Sep 07 19:24:56 2017

@author: hudjimartsu
"""

import os
import tarfile
from osgeo import gdal
from osgeo import gdal_array
import numpy as np
#Create output folder
newFolder = "LandsatData"
os.makedirs(newFolder)

#Extract files
tar = tarfile.open("LC81910182016153-SC20161208043748.tar.gz", "r:gz")
tar.extractall(newFolder)
tar.close()


filepath = r"C:/Users/hudjimartsu/Documents/PROJECT/FOREST 2020/TRAINING/PyQgis/DATA/LandsatData/LC81910182016153LGN00_sr_band4.tif"

# Open the file:
raster = gdal.Open(filepath)

# Check type of the variable 'raster'
type(raster)
# Projection
raster.GetProjection()

# Dimensions
raster.RasterXSize
raster.RasterYSize

# Number of bands
raster.RasterCount

# Metadata for the raster dataset
raster.GetMetadata()
# Read the raster band as separate variable
band = raster.GetRasterBand(1)

# Check type of the variable 'band'
type(band)

# Data type of the values
gdal.GetDataTypeName(band.DataType)
# Compute statistics if needed
if band.GetMinimum() is None or band.GetMaximum()is None:
    band.ComputeStatistics(0)
    print("Statistics computed.")

# Fetch metadata for the band
band.GetMetadata()

# Print only selected metadata:
print ("[ NO DATA VALUE ] = ", band.GetNoDataValue()) # none
print ("[ MIN ] = ", band.GetMinimum())
print ("[ MAX ] = ", band.GetMaximum())

rasterArray = raster.ReadAsArray()
rasterArray = gdal_array.LoadFile(filepath)
rasterArray.min()
rasterArray.max()

# Get nodata value from the GDAL band object
nodata = band.GetNoDataValue()

#Create a masked array for making calculations without nodata values
rasterArray = np.ma.masked_equal(rasterArray, nodata)
type(rasterArray)

# Check again array statistics
rasterArray.min()
rasterArray.max()
rasternew=rasterArray*2
rasternew.max()
#print os.getcwd()