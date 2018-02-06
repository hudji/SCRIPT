from osgeo import gdal
from osgeo import gdal_array
filepath = r"C:/Users/hudjimartsu/Documents/PROJECT/FOREST 2020/TRAINING/PyQgis/DATA/source/LC81220652016230LGN00_B7.TIF"

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

rasterArray = band.ReadAsArray()
rasterArray = gdal_array.LoadFile(filepath)