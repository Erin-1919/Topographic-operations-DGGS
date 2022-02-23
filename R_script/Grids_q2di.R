##################################################
### Generate hex grids within area of interest ###
##################################################

library(dggridR)
library(rgdal)
library(rgeos)
library(dplyr)
library(geosphere)

# read arguments - resolution level and study area
sys.argv = commandArgs(trailingOnly = TRUE)
dggs_res = as.numeric(sys.argv[1]) # 20-24
study_area = sys.argv[2] # Calgary / Canmore / BuffaloLake

# define orientation
vlat = 31
vlon = -56
vazimuth = 0

# define functions
# generate hex grids with q2di index attached
generate_grids = function(area,resolution){
  study_area = readOGR(dsn="Data",layer=area)
  study_area = bbox(study_area[1,])
  # get study area info
  minx = study_area[1,1]
  miny = study_area[2,1]
  maxx = study_area[1,2]
  maxy = study_area[2,2]
  # construct a look-up table, storing resolution and corresponding cell size and vertical resolution
  res_list = c(19,20,21,22,23,24)
  cell_size_list = c(0.0009,0.0008,0.0003,0.0003,0.0001,0.0001)  
  look_up = data.frame("res_list" = res_list,"cell_size_list" = cell_size_list)
  # look up cell size and vertical resolution
  dggs_cellsize = look_up$cell_size_list[look_up$res_list == resolution]
  # define DGGS
  v_lat = vlat
  v_lon = vlon
  azimuth = vazimuth
  DGG = dgconstruct(projection="ISEA", aperture=3, topology="HEXAGON", res=resolution, azimuth_deg=azimuth, pole_lat_deg=v_lat, pole_lon_deg=v_lon)
  # grid
  grid = dgrectgrid(DGG, minlat = miny, minlon = minx, maxlat = maxy, maxlon = maxx, cellsize = dggs_cellsize, frame=FALSE)
  grid_centroids = centroid(grid)
  grid = as(grid, "SpatialPolygonsDataFrame")
  grid@data$i = dgGEO_to_Q2DI(DGG,grid_centroids[,1], grid_centroids[,2])$i
  grid@data$j = dgGEO_to_Q2DI(DGG,grid_centroids[,1], grid_centroids[,2])$j
  writeOGR(grid, sprintf("Result/Area_%s_hex_%d.shp",area,resolution), layer = "cell", driver='ESRI Shapefile')
}

generate_grids(sprintf("Area_%s",study_area),dggs_res)

