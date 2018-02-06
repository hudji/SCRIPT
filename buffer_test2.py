from osgeo import gdal
import numpy as np,sys
raster_filepath='D:/FORESTS2020/TRAINING/PyQgis/RESULT/Landsat8/Cloud_Masking/C121217/Cloud_masking_CL.TIF'
data=gdal.Open(raster_filepath)
read=data.GetRasterBand(1).ReadAsArray()
#def raster_buffer(raster_filepath, dist=60):
def raster_buffer(data, dist=60):
     """This function creates a distance buffer around the given raster file with non-zero values.
     The value in output raster will have value of the cell to which it is close to."""
     """d=gdal.Open(raster_filepath)
     if d is None:
         print("Error: Could not open image " + raster_filepath)
         sys.exit(1)
     global proj,geotrans,row,col
     proj=d.GetProjection()
     geotrans=d.GetGeoTransform()
     row=d.RasterYSize
     col=d.RasterXSize
     inband=d.GetRasterBand(1)
     in_array = inband.ReadAsArray(0,0,col,row).astype(int)"""
     row= data.shape[0]
     col=data.shape[1]
     in_array = data.astype(int)
     """Xcell_size=int(abs(geotrans[1]))
     Ycell_size=int(abs(geotrans[5]))"""
     Xcell_size=30
     Ycell_size=30
     cell_size = (Xcell_size+Ycell_size)/2
     cell_dist=dist/cell_size
     #in_array[in_array == (inband.GetNoDataValue() or 0 or -999)]=0
     out_array=np.zeros_like(in_array)
     temp_array=np.zeros_like(in_array)
     i,j,h,k=0,0,0,0
     print("Running distance buffer...")
     while(h<col):
         k=0
         while(k<row):
             if(in_array[k][h]>=1):
                 i=h-cell_dist
                 while((i<cell_dist+h) and i<col):
                     j=k-cell_dist
                     while(j<(cell_dist+k) and j<row):
                         if(((i-h)**2+(j-k)**2)<=cell_dist**2):
                             if(temp_array[j][i]==0 or temp_array[j][i]>((i-h)**2+(j-k)**2)):
                                 out_array[j][i]= in_array[k][h]
                                 temp_array[j][i]=(i-h)**2+(j-k)**2
                         j+=1
                     i+=1
             k+=1
         h+=1
     d,temp_array,in_array=None,None, None
     return out_array

def export_array(in_array,output_path):
    """This function is used to produce output of array as a map."""
    driver = gdal.GetDriverByName("GTiff")
    outdata = driver.Create(output_path,col,row,1)
    outband=outdata.GetRasterBand(1)
    outband.SetNoDataValue(np.nan)
    outband.WriteArray(in_array)
    # Georeference the image
    outdata.SetGeoTransform(geotrans)
    # Write projection information
    outdata.SetProjection(proj)
    outdata.FlushCache()
    outdata = None

raster_buffer_array=raster_buffer('D:/FORESTS2020/TRAINING/PyQgis/RESULT/Landsat8/Cloud_Masking/C131217/Cloud_masking.TIF',60)
raster_buffer_array=raster_buffer(read,60)
export_array(raster_buffer_array,'D:/FORESTS2020/TRAINING/PyQgis/RESULT/Landsat8/Cloud_Masking/C131217/buff60array_I.tif')
print("Done")

myRange = [8, 11]
myArr = [6, 7, 8, 9, 10, 11, 12]
myArrFiltered = [x for x in myArr if myRange[0] <= x <= myRange[1]]
print(myArrFiltered)