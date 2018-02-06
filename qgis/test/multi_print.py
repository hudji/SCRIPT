import gdal
import numpy as np
from affine import Affine
import os
import glob
path=("D:\PROJECT\FOREST 2020\TRAINING\PyQgis\RESULT\Landsat8\L301017")
def save_multiband (output_name, dataset, raster_data, driver,
                    NaN_Value, nband, arr_projection=None):
    if arr_projection is None:
        arr_projection = []
    else:
        if str(type(arr_projection)) == "<class 'osgeo.gdal.Dataset'>":
            srs = arr_projection.GetProjectionRef()
        elif str(type(arr_projection)) == "<type 'int'>":
            srs = osr.SpatialReference()
            srs.ImportFromEPSG(arr_projection)
            srs = srs.ExportToPrettyWkt()
        else:
            srs = arr_projection
    NaN_rast = NaN_Value
    # Open the reference dataset
    g = dataset
    # Get the Geotransform vector
    if raster_data is False:
        raster_data = g.ReadAsArray()
        rast_arr = np.copy(raster_data)
    #auxiliar
    band = 0
    if type(raster_data) == tuple:
        rast_arr = np.array(raster_data[band])
    if str(type(g)) == "<class 'osgeo.gdal.Dataset'>":
        geo_transform = g.GetGeoTransform()
        x_size = g.RasterXSize  # Raster xsize
        y_size = g.RasterYSize  # Raster ysize
        srs = g.GetProjectionRef()  # Projection
    elif str(type(g)) == "<class 'affine.Affine'>":
        geo_transform = (g[2], g[0], g[1], g[5], g[3], g[4])
        rast_arr = raster_data[band,:,:]
        x_size = int(rast_arr.shape[1])
        y_size = int(rast_arr.shape[0])
    driver = gdal.GetDriverByName(driver)
    dataset_out = driver.Create(output_name, x_size, y_size, nband,
                                gdal.GDT_Float32)
    #end auxiliar
    for band in range(nband):
        if type(raster_data) == tuple:
            rast_arr = np.array(raster_data[band])
        if str(type(g)) == "<class 'osgeo.gdal.Dataset'>":
            geo_transform = g.GetGeoTransform()
            x_size = g.RasterXSize  # Raster xsize
            y_size = g.RasterYSize  # Raster ysize
            srs = g.GetProjectionRef()  # Projection
        elif str(type(g)) == "<class 'affine.Affine'>":
            geo_transform = (g[2], g[0], g[1], g[5], g[3], g[4])
            rast_arr = raster_data[band,:,:]
            x_size = int(rast_arr.shape[1])
            y_size = int(rast_arr.shape[0])
        #PROCESS RASTERIO NUMPY
        else:
            geo_transform = (g[1][2], g[1][0], g[1][1], g[1][5], g[1][3], g[1][4])
            rast_arr = np.array(g[0])
            x_size = int(rast_arr.shape[2])
            y_size = int(rast_arr.shape[1])
        rast_arr[rast_arr == NaN_rast] = np.NaN
        dataset_out.SetGeoTransform(geo_transform)
        dataset_out.SetProjection(srs)
        dataset_out.GetRasterBand(band + 1).WriteArray(
            rast_arr.astype(np.float32))

#assuming the data had the same extent and spatial resolution
B1 = gdal.Open(r'D:\PROJECT\FOREST 2020\TRAINING\PyQgis\DATA\Landsat8\NEW\band2.tif')
B2 = gdal.Open(r'D:\PROJECT\FOREST 2020\TRAINING\PyQgis\DATA\Landsat8\NEW\band5.tif')
B3 = gdal.Open(r'D:\PROJECT\FOREST 2020\TRAINING\PyQgis\DATA\Landsat8\NEW\band6.tif')
B4 = gdal.Open(r'D:\PROJECT\FOREST 2020\TRAINING\PyQgis\DATA\Landsat8\NEW\DEM\aspect.tif')
srs = B1.GetProjectionRef()
b1, b2, b3, b4 = (B1.ReadAsArray(), B2.ReadAsArray(),
                  B3.ReadAsArray(), B4.ReadAsArray())

def GetTransform(Transform):
    transform = Transform.GetGeoTransform()
    Aff = Affine(transform[1], transform[2], transform[0], transform[4],
                 transform[5], transform[3])
    return(Aff)

B = np.stack([B1, B2, B3, B4])
save_multiband(path+'\\SR.tif', GetTransform(B1), B, "GTiff", -9999, B.ndim + 1, srs)



