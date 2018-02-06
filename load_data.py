from osgeo import gdal
import glob
import os
in_directory=r'C:\Users\hudjimartsu\Documents\PROJECT\FOREST 2020\TRAINING\PyQgis\DATA\landsat'
file_to_process=glob.glob(os.path.join(in_directory, '*.TIF'))
for data_path in file_to_process:
    raster_dataset=gdal.Open(data_path, gdal.GA_ReadOnly)