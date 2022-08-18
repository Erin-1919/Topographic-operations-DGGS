# Multi-resolution topographic analysis in hexagonal Discrete Global Grid Systems -- Source Code

This work aims to develop analytical functions based on terrain data in a pure environment of ISEA3H DGGS. The developed operations included descriptive statistics, topographic analysis, and topographic indices, and were classified as three catogories, namely local, zonal, and focal. The open-sourced library [*dggridR*](https://github.com/r-barnes/dggridR) was used to complete conversion between geographic locations and ISEA3H DGGS cell indices. The experiment was carried out using a hybrid of Python 3.7.7 and R 3.6.2 environments. The code used to conduct the experiment are available in the folder [*R_script*](https://github.com/Erin-1919/Topographic-operations-DGGS/tree/main/R_script) and [*Python_script*](https://github.com/Erin-1919/Topographic-operations-DGGS/tree/main/Python_script).

## Manuscript Information
### Title of Manuscript
Multi-resolution topographic analysis in hexagonal Discrete Global Grid Systems

### Keywords
Discrete Global Grid Systems; topographical analysis; multi-resolution; map algebra

### DOI
10.1016/j.jag.2022.102985

### Authors
Mingke Li, Heather McGrath, and Emmanuel Stefanakis

### Corresponding Author
[Mingke Li](https://erin-1919.github.io/) (mingke.li@ucalgary.ca)

[ORCID](https://orcid.org/0000-0001-6310-4964)

### Abstract
Discrete Global Grid Systems (DGGS) have been increasingly adopted as a standard framework for multi-source geospatial data. Previous research largely studied the mathematical foundation of discrete global grids, developed open-source libraries, and explored their application as data integration platforms. This study investigated the multi-resolution terrain analysis in a pure hexagonal DGGS environment, including descriptive statistics, topographic parameters, and topographic indices. Experiments across multiple grid resolutions were carried out in three study areas with different terrain roughness in Alberta, Canada. Five algorithms were proposed to calculate both the slope gradient and terrain aspect. A cell-based pair-wise comparison showed a strong positive correlation between the gradient values as calculated from five algorithms. The grid resolutions as well as the terrain roughness had a clear effect on the computed slope gradient and topographic indices. This research aims to enhance the analytical functionality of hexagonal DGGS to better support decision-making in real world problems.  

### Code Repository
https://github.com/Erin-1919/Topographic-operations-DGGS

## Libraries used
*Python*
 - numpy 1.19.4
 - scipy 1.5.3
 - rasterio 1.2.1
 - gdal 3.1.4
 - pandas 1.1.4
 - anytree 2.8.0
 - multiprocess 0.70.12.2

*R*
 - dggridR 2.0.4
 - rgdal 1.5.16
 - rgeos 0.5.5
 - dplyr 1.0.2
 - geosphere 1.5.10

## Data availability
The original Canadian Digital Elevation Model (CDEM) data can be downloaded via the Geospatial-Data Extraction tool in [Canada's Open Government Portal](https://maps.canada.ca/czs/index-en.html), or they can be obtained through the STAC API. The 2020 Annual Crop Inventory produced by the Agriculture and Agri-Food Canada, used as an example in zonal statistics, is accessible [here](https://open.canada.ca/data/en/dataset/ba2645d5-4458-414d-b196-6303ac06c1c9). 
