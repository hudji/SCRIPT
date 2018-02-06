import numpy as np
import glob
from osgeo import gdal
import numexpr, sys
from scipy import ndimage
#Load data raster
raster_list=glob.glob('D:\FORESTS2020\TRAINING\PyQgis\DATA\source\*.tif')
read=[]
for i in raster_list:
    band=gdal.Open(i)
    read.append(band.GetRasterBand(1).ReadAsArray())

bqa= np.array(read[11], dtype=float)
Cloud_cirrus2 = numexpr.evaluate("((bqa >= 53248))") # https://landsat.usgs.gov/qualityband

def raster_buffer(raster_filepath, dist=100):
    """This function creates a distance buffer around the given raster file with non-zero values.
    The value in output raster will have value of the cell to which it is close to."""
    d = gdal.Open(raster_filepath)
    if d is None:
        print("Error: Could not open image " + raster_filepath)
        sys.exit(1)
    global proj, geotrans, row, col
    proj        = d.GetProjection()
    geotrans    = d.GetGeoTransform()
    row         = d.RasterYSize
    col         = d.RasterXSize
    inband      = d.GetRasterBand(1)
    in_array    = inband.ReadAsArray(0, 0, col, row).astype(int)
    Xcell_size  = int(abs(geotrans[1]))
    Ycell_size  = int(abs(geotrans[5]))
    cell_size   = (Xcell_size + Ycell_size) / 2
    cell_dist   = dist / cell_size
    in_array[in_array == (inband.GetNoDataValue() or 0 or -999)] = 0
    out_array   = np.zeros_like(in_array)
    temp_array  = np.zeros_like(in_array)
    i, j, h, k  = 0, 0, 0, 0
    print("Running distance buffer, please wait...")
    while (h < col):
        k = 0
        while (k < row):
            if (in_array[k][h] >= 1):
                i = h - cell_dist
                while ((i < cell_dist + h) and i < col):
                    j = k - cell_dist
                    while (j < (cell_dist + k) and j < row):
                        if (((i - h) ** 2 + (j - k) ** 2) <= cell_dist ** 2):
                            if (temp_array[j][i] == 0 or temp_array[j][i] > ((i - h) ** 2 + (j - k) ** 2)):
                                out_array[j][i] = in_array[k][h]
                                temp_array[j][i] = (i - h) ** 2 + (j - k) ** 2
                        j += 1
                    i += 1
            k += 1
        h += 1
    d, temp_array, in_array = None, None, None
    return out_array

def export_array(in_array, output_path):
    """This function is used to produce output of array as a map."""
    driver  = gdal.GetDriverByName("GTiff")
    outdata = driver.Create(output_path, col, row, 1)
    outband = outdata.GetRasterBand(1)
    outband.SetNoDataValue(np.nan)
    outband.WriteArray(in_array)
    # Georeference the image
    outdata.SetGeoTransform(geotrans)
    # Write projection information
    outdata.SetProjection(proj)
    outdata.FlushCache()
    outdata = None

raster_buffer_array = raster_buffer("D:/FORESTS2020/TRAINING/PyQgis/RESULT/Landsat8/Cloud Masking/C051217/Cloud_test.tif", 100)
im_med = ndimage.median_filter(raster_buffer_array, 3)
export_array(im_med, "D:/FORESTS2020/TRAINING/PyQgis/RESULT/Landsat8/Cloud Masking/C051217/Buffer100_filter3_I.tif")
export_array(im_med, 'D:/FORESTS2020/TRAINING/PyQgis/RESULT/Landsat8/Cloud Masking/C051217/' + Cloud_cirrus2 +'_coba.tif')
export_array(Cloud_cirrus2, "D:/FORESTS2020/TRAINING/PyQgis/RESULT/Landsat8/Cloud Masking/C051217/Cloud3.tif")
