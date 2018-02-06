import gdal, os
from gdalconst import *
filename = "D:\PROJECT\FOREST 2020\TRAINING\PyQgis\DATA\Landsat8\clip\slope.tif" #path to raster
def findGDALCoordinates():
    data = gdal.Open(filename,GA_ReadOnly)
    if data is None:
        return []
    geoTransform = data.GetGeoTransform()
    minx = geoTransform[0]
    maxy = geoTransform[3]
    maxx = minx + geoTransform[1]*data.RasterXSize
    miny = maxy + geoTransform[5]*data.RasterYSize
    return [minx,miny,maxx,maxy]


from xml.dom import minidom
xmldoc = minidom.parse('C:\Users\sahid\Pictures\S2A_MSIL1C_20170922T025541_N0205_R032_T48MYT_20170922T031450\S2A_MSIL1C_20170922T025541_N0205_R032_T48MYT_20170922T031450.SAFE\MTD_MSIL1C.xml')
itemlist = xmldoc.getElementsByTagName("n1:General_Info")[0]
root=itemlist.getElementsByTagName("Product_Info")[0]
coba=root.getElementsByTagName("Product_Info")[0].firstChild.data

for i in xmldoc:
    dec=xmldoc.getElementsByTagName("n1:General_Info")[0]
    print (dec)

collection = xmldoc.documentElement
if collection.hasAttribute("n1:General_Info"):
   print "Root element : %s" % collection.getAttribute("n1:General_Info")


collection
for date in itemlist:
    print date
print(xmldoc.toxml())
print(itemlist)

print xmldoc.nodeName
print xmldoc.getfirstChild.tagName
