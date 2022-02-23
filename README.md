# Topographic operations in hexagonal Discrete Global Grid Systems -- Source Code

This work aims to develop analytical functions based on terrain data in a pure environment of ISEA3H DGGS. The developed operations included descriptive statistics, topographic analysis, and topographic indices, and were classified as three catogories, namely local, zonal, and focal. The open-sourced library [*dggridR*](https://github.com/r-barnes/dggridR) was used to complete conversion between geographic locations and ISEA3H DGGS cell indices. The experiment was carried out using a hybrid of Python 3.7.7 and R 3.6.2 environments. The code used to conduct the experiment are available in the folder [*R_script*](https://github.com/Erin-1919/Topographic-operations-DGGS/tree/main/R_script) and [*Python_script*](https://github.com/Erin-1919/Topographic-operations-DGGS/tree/main/Python_script).

## Manuscript Information
### Title of Manuscript
Topographic operations in hexagonal Discrete Global Grid Systems

### Keywords
Discrete Global Grid Systems; topographical analysis; multi-resolution; map algebra

### DOI
todo

### Authors
Mingke Li, Heather McGrath, and Emmanuel Stefanakis

### Corresponding Author
[Mingke Li](https://erin-1919.github.io/) (mingke.li@ucalgary.ca)

[ORCID](https://orcid.org/0000-0001-6310-4964)

### Abstract
Discrete Global Grid Systems (DGGS) have been increasingly adopted as the framework of multi-source geospatial data. Previous research largely studied the mathematical foundation, developed open-sourced DGGS libraries, and explored their application as integration platforms. This study investigated the analytical operations in a pure hexagonal DGGS environment, including descriptive statistics, topographic analysis, and topographic indices. Experiments across multiple resolutions were carried out in three areas with various roughness in Alberta, Canada. With five algorithms proposed to calculate slope gradient and aspect, the cell-based, pair-wise comparison showed strong positive relationships between the gradient resulted from various algorithms. Resolutions influenced the detection of elevation changes and the rate of changes, and the degree of such influence also depended on the terrain roughness. This research sets the stage for the analytical development of general DGGS and helps to bridge the gap between the existing DGGS implementations and DGGS-driven decision-making in the real world. 

### Code Repository
https://github.com/Erin-1919/Topographic-operations-DGGS

## Libraries used
*Python*
 - numpy
 - math
 - scipy
 - warnings
 - rasterio
 - gdal
 - pandas
 - anytree
 - statistics
 - itertools
 - multiprocess

*R*
 - dggridR
 - rgdal
 - rgeos
 - dplyr
 - geosphere

## Data availability
The original Canadian Digital Elevation Model (CDEM) data can be downloaded via the Geospatial-Data Extraction tool in [Canada's Open Government Portal](https://maps.canada.ca/czs/index-en.html), or they can be obtained through the STAC API. The 2020 Annual Crop Inventory produced by the Agriculture and Agri-Food Canada, used as an example in zonal statistics, is accessible [here](https://open.canada.ca/data/en/dataset/ba2645d5-4458-414d-b196-6303ac06c1c9). 
