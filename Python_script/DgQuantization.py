import pandas, numpy, rasterio
import time, sys, os, warnings
import multiprocess as mp
import DgBaseFunc
from scipy import interpolate

warnings.simplefilter('error', RuntimeWarning) 

# set resolution level and study area
dggs_res = int(sys.argv[1]) # 20-24
area = sys.argv[2] # Calgary / Canmore / BuffaloLake

# look up cellsize and vertical resolution
look_up = DgBaseFunc.look_up_table()
dggs_cellsize = look_up.loc[dggs_res,'cell_size']
vertical_res = look_up.loc[dggs_res,'verti_res']

## Define the main funcitons    
def find_neighbor(x,y,DEM_TIF):
    ''' 
    Find neighbors for interpolation --
    determine the DEM source
    find out the 4 neighbors geographic coords
    extract the elevations at these 4 coords
    convert 4 coords to array index then back to 4 coords of grid mesh centers
    '''
    x_index,y_index = rasterio.transform.rowcol(DEM_TIF.transform, x,y)
    xc,yc = rasterio.transform.xy(DEM_TIF.transform, x_index, y_index)
    if x > xc and y > yc:
        x_index_array = [x_index-1,x_index-1,x_index,x_index]
        y_index_array = [y_index,y_index+1,y_index,y_index+1]
    elif x > xc and y < yc:
        x_index_array = [x_index,x_index,x_index+1,x_index+1]
        y_index_array = [y_index,y_index+1,y_index,y_index+1]
    elif x < xc and y > yc:
        x_index_array = [x_index-1,x_index-1,x_index,x_index]
        y_index_array = [y_index-1,y_index,y_index-1,y_index]
    elif x < xc and y < yc:
        x_index_array = [x_index,x_index,x_index+1,x_index+1]
        y_index_array = [y_index-1,y_index,y_index-1,y_index]
    x_array,y_array = rasterio.transform.xy(DEM_TIF.transform, x_index_array, y_index_array)
    coords = [(lon,lat) for lon, lat in zip(x_array,y_array)]
    z_array = [elev[0] for elev in DEM_TIF.sample(coords)]
    return x_array, y_array, z_array

def dggs_elevation_cdem(x,y,interp = 'linear'):
    ''' 
    Resample CDEM -- 
    if an error is raised then return -32767 as its final elevation
    if the point or any of its neighbors has the value -32767 then return -32767 as its final elevation
    if none of its neighbor has value -32767 then interpolate elevation
    restrict the decimal places according to the look-up table defined earlier 
    '''
    DEM_TIF = CDEM_TIF
    try:
        x_array, y_array, z_array = find_neighbor(x,y,DEM_TIF)
        if -32767 in z_array:
            return -32767
        else:
            CDEM_interp = interpolate.interp2d(x_array, y_array, z_array, kind=interp)
            elevation = CDEM_interp(x,y)[0]
            elevation = round(elevation,vertical_res)
            return elevation
    except:
        return -32767

def dggs_elevation_df(dataframe):
    ''' a function on the dataframe '''
    dataframe['model_elev'] = [dggs_elevation_cdem(lon,lat) for lon, lat in zip(dataframe.lon_c, dataframe.lat_c)]
    dataframe = dataframe.drop(columns=['lon_c','lat_c'])
    return dataframe

if __name__=='__main__':
    
    input_csv_path = 'Result/Area_{}_centroids_{}.csv'.format(area,dggs_res)
    centroid_df = pandas.read_csv(input_csv_path, sep=',')

    ## Read DEMs -- CDEM 
    CDEM_TIF = rasterio.open('Data/{}.tif'.format(area))

    # record timing -- start
    start_time = time.time()

    # call the function by parallel processing
    n_cores = int(os.environ.get('SLURM_CPUS_PER_TASK',default=1))
    centroid_df_split = numpy.array_split(centroid_df, n_cores)
    pool = mp.Pool(processes = n_cores)
    centroid_df_output = pandas.concat(pool.map(dggs_elevation_df, centroid_df_split))
    pool.close()
    pool.join()

    # non-parallel
    # centroid_df_output = dggs_elevation_df(centroid_df)

    # record timing -- end
    print (dggs_res)
    print ("Processing time: %s seconds" % (time.time() - start_time))

    # save the results as csv
    output_csv_path = 'Result/{}_elev_{}.csv'.format(area,dggs_res)
    centroid_df_output.to_csv(output_csv_path, index=False)
