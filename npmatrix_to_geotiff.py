import numpy as np
from osgeo import gdal, gdal_array, osr

def npmatrix_to_geotiff (filepath, matrices, gdal_type, transform = None, projection = None, nodata = None):
  (y_res, x_res) = matrices[0].shape
  driver = gdal.GetDriverByName('GTiff')
  image = driver.Create(filepath, x_res, y_res, len(matrices), gdal_type)
  if transform is not None:
    image.SetGeoTransform(transform)
  if projection is not None:
    image.SetProjection(projection)
  for index, matrix in enumerate(matrices):
    index = index + 1
    band = image.GetRasterBand(index)
    if nodata is not None:
      band.SetNoDataValue(nodata)
    band.WriteArray(matrix)
    band.FlushCache
  return image

# use it like this
# note that currently geotiff only supports a single nodata value for all bands
# in the future it may support more
# all the matrices are expected to be the same shape
# all the matrices must be 2 dimensional arrays
# if you want the gdal geotiff to render exactly like the numpy array construction
# you may need to flip it vertically first using [::-1]

projection = osr.SpatialReference()
projection.ImportFromEPSG(4326)

npmatrix_to_geotiff(
  './image.tif',
  [np.arange(0,10).reshape(2,5),
  gdal_array.NumericTypeCodeToGDALTypeCode(np.uint8),
  (0, 1, 0, 0, 0, 1),
  projection.ExportToWkt(),
  None
)

