# The code for erta_ale validataion
# Jae Sung Kim
# Some code block was referered from http://pcjericks.github.io/py-gdalogr-cookbook/raster_layers.html#raster-to-vector-line
from osgeo import gdal, osr, ogr
import sys, os, json, csv
import numpy as np
from affine import Affine
import copy
from osgeo.gdalconst import *
import subprocess
import itertools
from math import sqrt, ceil, pow
import fiona
from shapely.geometry import shape, LineString, Point, MultiLineString
from shapely.wkt import dumps, loads
# register all of the GDAL drivers
gdal.AllRegister()

def raster2array(rasterfn):
	raster = gdal.Open(rasterfn)
	band = raster.GetRasterBand(1)
	return band.ReadAsArray()

def getPixelCoordinate(input_raster, x, y):
	input = gdal.Open(input_raster)
	trs = input.GetGeoTransform()

	forward_transform = Affine.from_gdal(trs[0], trs[1], trs[2], trs[3], trs[4], trs[5])
	reverse_transform = ~forward_transform
	px, py = reverse_transform * (float(x), float(y))

	return int(round(px,1)), int(round(py,1))


def pixelOffset2coord(rasterfn,xOffset,yOffset):
    raster = gdal.Open(rasterfn)
    geotransform = raster.GetGeoTransform()
    originX = geotransform[0]
    originY = geotransform[3]
    pixelWidth = geotransform[1]
    pixelHeight = geotransform[5]
    coordX = originX+pixelWidth*xOffset
    coordY = originY+pixelHeight*yOffset
    return coordX, coordY


def array2shp(array,outSHPfn,rasterfn,pixelValue):

    # max distance between points
    raster = gdal.Open(rasterfn)
    geotransform = raster.GetGeoTransform()
    pixelWidth = geotransform[1]
    maxDistance = ceil(sqrt(2*pixelWidth*pixelWidth))
    print(maxDistance)

    # array2dict
    count = 0
    roadList = np.where(array == pixelValue)
    multipoint = ogr.Geometry(ogr.wkbMultiLineString)
    pointDict = {}
    for indexY in roadList[0]:
        indexX = roadList[1][count]
        Xcoord, Ycoord = pixelOffset2coord(rasterfn,indexX,indexY)
        pointDict[count] = (Xcoord, Ycoord)
        count += 1

    # dict2wkbMultiLineString
    multiline = ogr.Geometry(ogr.wkbMultiLineString)
    for i in itertools.combinations(pointDict.values(), 2):
        point1 = ogr.Geometry(ogr.wkbPoint)
        point1.AddPoint(i[0][0],i[0][1])
        point2 = ogr.Geometry(ogr.wkbPoint)
        point2.AddPoint(i[1][0],i[1][1])

        distance = point1.Distance(point2)

        if distance < maxDistance:
            line = ogr.Geometry(ogr.wkbLineString)
            line.AddPoint(i[0][0],i[0][1])
            line.AddPoint(i[1][0],i[1][1])
            multiline.AddGeometry(line)

    # wkbMultiLineString2shp
    shpDriver = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(outSHPfn):
        shpDriver.DeleteDataSource(outSHPfn)
    outDataSource = shpDriver.CreateDataSource(outSHPfn)
    outLayer = outDataSource.CreateLayer(outSHPfn, geom_type=ogr.wkbMultiLineString )
    featureDefn = outLayer.GetLayerDefn()
    outFeature = ogr.Feature(featureDefn)
    outFeature.SetGeometry(multiline)
    outLayer.CreateFeature(outFeature)
    line_WKT = multiline.ExportToWkt()
    return [line_WKT,maxDistance]


resolution = 71


input_dem = 'ASTER_DEM_resampled.tif'
cmd_0 = "gdalwarp {0} {1} -tr {2} {2} -r bilinear -overwrite".format(sys.argv[1], input_dem, resolution)
subprocess.call(cmd_0,shell=True)	
#input_dem=sys.argv[1]
input_fill = 'dem_fill.tif'
slope_raster = 'slope.tif'
input_raster = 'dir_2.tif'

cmd_1 = "/home/parallels/Documents/TauDEM-Develop/bin/pitremove -z {0} -fel {1}".format(input_dem, input_fill)
subprocess.call(cmd_1,shell=True)
cmd_2 = "/home/parallels/Documents/TauDEM-Develop/bin/d8flowdir -fel {0} -sd8 {1} -p {2}".format(input_fill, slope_raster, input_raster)
subprocess.call(cmd_2,shell=True)

input_x = sys.argv[2]
input_y = sys.argv[3]
output_raster = sys.argv[4]+'.tif'

[x_l, y_l] = getPixelCoordinate(input_raster, input_x, input_y)


rasterArray = raster2array(input_raster)
row=len(rasterArray)
col=len(rasterArray[0])
out_raster=copy.deepcopy(rasterArray)

c=0

while x_l<col-1 and x_l>0 and y_l>0 and y_l<row-1:
	

	h=rasterArray[int(y_l),int(x_l)]
	
	
	if h==1:
		
		x_l=x_l+1
		y_l=y_l
	
	elif h==2:
		x_l=x_l+1
		y_l=y_l-1
	
	elif h==3:
		x_l=x_l
		y_l=y_l-1
	
	elif h==4:
		x_l=x_l-1
		y_l=y_l-1
	
	elif h==5:
		x_l=x_l-1
		y_l=y_l
	
	elif h==6:
		x_l=x_l-1
		y_l=y_l+1

	elif h==7:
		x_l=x_l
		y_l=y_l+1

	elif h==8:
		x_l=x_l+1
		y_l=y_l+1
	else:
		print(h)

	out_raster[y_l][x_l]=255

#create output image
inDs = gdal.Open(input_raster)
driver = inDs.GetDriver()
outDs = driver.Create(output_raster, col, row, 1, GDT_Int16)
outDs.SetGeoTransform(inDs.GetGeoTransform())
outDs.SetProjection(inDs.GetProjection())

outBand = outDs.GetRasterBand(1)
outBand.WriteArray(out_raster)
outBand.FlushCache()
 
outDs = None

array = raster2array(output_raster)
output_shp = sys.argv[4]+'.shp'
output_prj = sys.argv[4]+'.prj'
result = array2shp(array,output_shp,output_raster,255)

print(output_shp)
print(output_raster)

spatialRef = osr.SpatialReference()
spatialRef.ImportFromEPSG(32637)
spatialRef.MorphToESRI()
file = open(output_prj, 'w')
file.write(spatialRef.ExportToWkt())
file.close()


