##################################################
### Generate centroids within area of interest ###
##################################################

library(dggridR)
library(rgdal)
library(rgeos)
library(dplyr)
library(geosphere)

# read arguments - resolution level and study area
sys.argv = commandArgs(trailingOnly = TRUE)
dggs_res = as.numeric(sys.argv[1])
study_area = sys.argv[2]

# define orientation
vlat = 31
vlon = -56
vazimuth = 0

# define functions
generate_centroids = function(area,resolution){
  study_area = readOGR(dsn="Data",layer=area)
  study_area = bbox(study_area[1,])
  # get study area info
  minx = study_area[1,1]
  miny = study_area[2,1]
  maxx = study_area[1,2]
  maxy = study_area[2,2]
  # construct a look-up table, storing resolution and corresponding cell size
  res_list = c(19,20,21,22,23,24)
  cell_size_list = c(0.0009,0.0008,0.0003,0.0003,0.0001,0.0001)
  look_up = data.frame("res_list" = res_list,"cell_size_list" = cell_size_list)
  # look up cell size
  dggs_cellsize = look_up$cell_size_list[look_up$res_list == resolution]
  # define DGGS
  v_lat = vlat
  v_lon = vlon
  azimuth = vazimuth
  DGG = dgconstruct(projection="ISEA", aperture=3, topology="HEXAGON", res=resolution, azimuth_deg=azimuth, pole_lat_deg=v_lat, pole_lon_deg=v_lon)
  # centroids
  coords = matrix(c(minx, miny, maxx, miny, maxx, maxy, maxx, miny, minx, miny), ncol = 2, byrow = TRUE)
  regbox = Polygon(coords)
  regbox = SpatialPolygons(list(Polygons(list(regbox), ID = "a")))
  centroids = sp::makegrid(regbox, cellsize = dggs_cellsize)
  centroids$quad = dgGEO_to_Q2DI(DGG,centroids$x1, centroids$x2)$quad
  centroids$i = dgGEO_to_Q2DI(DGG,centroids$x1, centroids$x2)$i
  centroids$j = dgGEO_to_Q2DI(DGG,centroids$x1, centroids$x2)$j
  centroids = centroids[!duplicated(centroids[c("quad","i","j")]),]
  centroids$lon_c = dgQ2DI_to_GEO(DGG,centroids$quad,centroids$i,centroids$j)$lon_deg
  centroids$lat_c = dgQ2DI_to_GEO(DGG,centroids$quad,centroids$i,centroids$j)$lat_deg
  print (unique(centroids[c("quad")]))
  centroids = subset(centroids, select = -c(x1,x2,quad))
  write.csv(centroids,sprintf("Result/%s_centroids_%d.csv",area,resolution), row.names = FALSE)
}

generate_centroids(sprintf("Area_%s",study_area),dggs_res)
