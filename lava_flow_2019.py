# The code for erta_ale 
# Jae Sung Kim
# Some code block was referered from http://pcjericks.github.io/py-gdalogr-cookbook/raster_layers.html#raster-to-vector-line
# python3 lava_flow_2019.py dem_2016_32637.tif 683250.0 1502460.0 verify 2019_main_c.shp



from osgeo import gdal, osr, ogr
import sys, os, json, csv, shlex
import numpy as np
from affine import Affine
import copy
from osgeo.gdalconst import *
import subprocess, time
import itertools
from math import sqrt, ceil, pow
import fiona
from shapely.geometry import shape, LineString, Point, MultiLineString
from shapely.wkt import dumps, loads

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

    raster = gdal.Open(rasterfn)
    geotransform = raster.GetGeoTransform()
    pixelWidth = geotransform[1]
    #maxDistance = ceil(sqrt(2*pixelWidth*pixelWidth))
    pixelHeight = geotransform[5]
    maxDistance = ceil(np.sqrt(pixelWidth**2+pixelHeight**2))
    

    count = 0
    roadList = np.where(array == pixelValue)
    multipoint = ogr.Geometry(ogr.wkbMultiLineString)
    pointDict = {}
    for indexY in roadList[0]:
        indexX = roadList[1][count]
        Xcoord, Ycoord = pixelOffset2coord(rasterfn,indexX,indexY)
        pointDict[count] = (Xcoord, Ycoord)
        count += 1

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

resolution = list(range(30,121,1))
input_layer = ogr.Open(sys.argv[5])
i_layer = input_layer.GetLayer()
pt_list = []

for feature in i_layer:
    geom = feature.GetGeometryRef()
    n_pt = geom.GetPointCount()
    for i in range(n_pt):
    	pt_co=[geom.GetPoint(i)[0], geom.GetPoint(i)[1]]
    	if pt_co not in pt_list:
    		pt_list.append(pt_co)

dist_max = {}
dist_list = {}
max_dist_list = {}

for i in range(len(resolution)): 
	print("      ")
	print("resultion "+str(resolution[i])+" is running")
	if i == 0:
		input_dem = sys.argv[1]
	else:
		input_dem = 'ASTER_DEM_Clipped'+str(i)+'.tif'
		cmd_0 = "gdalwarp {0} {1} -tr {2} {2} -r bilinear -overwrite".format(sys.argv[1], input_dem, resolution[i])
		args_0 = shlex.split(cmd_0)
		p1=subprocess.Popen(args_0)	
		p1.wait()
	
	
	input_fill = 'dem_fill_'+str(i)+'.tif'
	slope_raster = 'slope_'+str(i)+'.tif'
	input_raster = 'dir_'+str(i)+'.tif'

	cmd_1 = "/home/parallels/Documents/TauDEM-Develop/bin/pitremove -z {0} -fel {1}".format(input_dem, input_fill)
	args_1 = shlex.split(cmd_1)
	p2=subprocess.Popen(args_1)
	p2.wait()

	cmd_2 = "/home/parallels/Documents/TauDEM-Develop/bin/d8flowdir -fel {0} -sd8 {1} -p {2}".format(input_fill, slope_raster, input_raster)
	args_2 = shlex.split(cmd_2)
	p3=subprocess.Popen(args_2)
	p3.wait()	

	input_x = sys.argv[2]
	input_y = sys.argv[3]
	output_raster = sys.argv[4]+str(i)+'.tif'

	[x_l, y_l] = getPixelCoordinate(input_raster, input_x, input_y)


	rasterArray = raster2array(input_raster)
	row=len(rasterArray)
	col=len(rasterArray[0])
	out_raster=copy.deepcopy(rasterArray)

	c=0
	
	# JK array value chages in loop. needs to be updated.
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
			print("else")
			#print(h)

		out_raster[y_l][x_l]=255

	
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
	output_shp = sys.argv[4]+str(i)+'.shp'
	output_prj = sys.argv[4]+str(i)+'.prj'
	result = array2shp(array,output_shp,output_raster,255)
	line_WKT = result[0]
	max_dist_list[resolution[i]] = result[1]
	line_shapely = loads(line_WKT)
	
	sum_dist = 0
	dist_array=[]
	for j in range(len(pt_list)): 
		pt = Point(pt_list[j][0],pt_list[j][1])
		dist=pt.distance(line_shapely)
		
		sum_dist = sum_dist + dist
		dist_array.append(dist)


	dist_ave = sum_dist / len(pt_list)
	rmse=np.sqrt(sum([(x-dist_ave)**2 for x in dist_array])/len(pt_list))
	dist_list[resolution[i]] = rmse
	dist_max[resolution[i]] = max(dist_array)

	spatialRef = osr.SpatialReference()
	spatialRef.ImportFromEPSG(32637)
	spatialRef.MorphToESRI()
	file = open(output_prj, 'w')
	file.write(spatialRef.ExportToWkt())
	file.close()

for i in range(len(resolution)): 
	
	if i == 0:
		input_dem = sys.argv[1]
		input_fill = 'dem_fill_'+str(i)+'.tif'
		if os.path.exists(input_fill):
			os.remove(input_fill)
		slope_raster = 'slope_'+str(i)+'.tif'
		if os.path.exists(slope_raster):
			os.remove(slope_raster)
		input_raster = 'dir_'+str(i)+'.tif'
		if os.path.exists(input_raster):
			os.remove(input_raster)
		erta_ale_tif='erta_ale'+str(i)+'.tif'
		if os.path.exists(erta_ale_tif):
			os.remove(erta_ale_tif)
	else:
		input_dem = 'ASTER_DEM_Clipped'+str(i)+'.tif'
		if os.path.exists(input_dem):
			os.remove(input_dem)
		input_fill = 'dem_fill_'+str(i)+'.tif'
		if os.path.exists(input_fill):
			os.remove(input_fill)
		slope_raster = 'slope_'+str(i)+'.tif'
		if os.path.exists(slope_raster):
			os.remove(slope_raster)
		input_raster = 'dir_'+str(i)+'.tif'
		if os.path.exists(input_raster):
			os.remove(input_raster)
		erta_ale_tif='erta_ale'+str(i)+'.tif'
		if os.path.exists(erta_ale_tif):
			os.remove(erta_ale_tif)
		
min_dist=min(dist_list.values())
for key, value in dist_list.items(): 
	if value==min_dist:
		opt_resolution=key
		print("resolution is"+str(key)+"m: average distance is "+str(value)+"(m)")	

f = open('lava_distance_2019.csv', 'w')
csv_file = csv.writer(f)
csv_file.writerow(["optimal resolution, average distance(m)"])
csv_file.writerow([str(opt_resolution)+","+str(min_dist)])
csv_file.writerow(["resolution, average distance(m)"])

for item in dist_list:
	csv_file.writerow([str(item)+","+str(dist_list[item])])

f.close()


