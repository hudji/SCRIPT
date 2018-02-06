from osgeo import gdal
from numpy import gradient
from numpy import pi
from numpy import arctan
from numpy import arctan2
from numpy import sin
from numpy import cos
from numpy import sqrt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import subprocess
def hillshade(array, azimuth, angle_altitude):
    x, y = gradient(array)
    slope = pi/2. - arctan(sqrt(x*x + y*y))
    aspect = arctan2(-x, y)
    azimuthrad = azimuth * pi / 180.
    altituderad = angle_altitude * pi / 180.
    shaded = sin(altituderad) * sin(slope)+ cos(altituderad) * cos(slope)* cos(azimuthrad - aspect)
    return 255*(shaded + 1)/2

ds = gdal.Open('D:\PROJECT\FOREST 2020\DATA\srtm\strm_gede.tif')
arr = ds.ReadAsArray()
hs_array = hillshade(arr, 90, 45)
plt.imshow(hs_array,cmap='Greys')
plt.savefig('D:\PROJECT\FOREST 2020\Paper\UK\Sketchup\hillshade_whistler2.tif')
plt.show()
# gdal command line tool called gdaldem
# link http://www.gdal.org/gdaldem.html
# usage:
# gdaldem hillshade input_dem output_hillshade
# [-z ZFactor (default=1)] [-s scale* (default=1)]"
# [-az Azimuth (default=315)] [-alt Altitude (default=45)]
# [-alg ZevenbergenThorne] [-combined]
# [-compute_edges] [-b Band (default=1)] [-of format] [-co
#"NAME=VALUE"]* [-q]
create_hillshade = '''gdaldem hillshade -az 315 -alt 45 ds D:\PROJECT\FOREST 2020\Paper\UK\Sketchup\hillshade_whistler3.tif'''
subprocess.call(create_hillshade)


import subprocess
# SLOPE
# - To generate a slope map from any GDAL-supported
elevation raster:
# gdaldem slope input_dem output_slope_map"
# [-p use percent slope (default=degrees)]
[-s scale* (default=1)]
# [-alg ZevenbergenThorne]
# [-compute_edges] [-b Band (default=1)] [-of format]
[-co "NAME=VALUE"]* [-q]
create_slope = '''gdaldem slope ds D:\PROJECT\FOREST 2020\Paper\UK\Sketchup\slope_w-degrees.tif '''
subprocess.call(create_slope)


######
import numpy as np
import glob
from osgeo import gdal
import subprocess
raster_list_dem=glob.glob('D:\PROJECT\FOREST 2020\DATA\srtm\*.tif')
read2=[]
for d in raster_list_dem:
    band=gdal.Open(d)
    read2.append(band.GetRasterBand(1).ReadAsArray())

dem=np.array(read2[1], dtype=float)
x, y = gradient(dem)
slope = pi/2. - arctan(sqrt(x*x + y*y))


#here is your 3 by 3 array, with its indexes
for i in range(3):
    for j in range(3):
        print((i, j), end='')
    print()








#export
geo = band.GetGeoTransform()
proj = band.GetProjection()
shape = dem.shape
driver = gdal.GetDriverByName("GTiff")
dst_ds = driver.Create("D:\PROJECT\FOREST 2020\TRAINING\PyQgis\RESULT\Landsat8\L011117\slope.tif", shape[1], shape[0], 1, gdal.GDT_Float32)
dst_ds.SetGeoTransform(geo)
dst_ds.SetProjection(proj)
dst_ds.GetRasterBand(1).WriteArray(slope)
dst_ds = None  # save, close"""