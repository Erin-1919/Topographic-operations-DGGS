import rasterio, gdal
import pandas, numpy
import time, sys, os
import multiprocess as mp
import DgBaseFunc

# set resolution level
dggs_res = int(sys.argv[1]) # 20-24
area = sys.argv[2] # Calgary / Canmore / BuffaloLake

# look up cell area
look_up = DgBaseFunc.look_up_table()
cell_area = look_up.loc[dggs_res,'cell_area'] * 1000000

## Define the main funcitons  
def reverse_projection(input_dem,output_dem):
    ''' A function to reversely project raster to the NAD83 CSRS) '''
    raster = rasterio.open(input_dem)
    wrap_option = gdal.WarpOptions(format = 'GTiff', 
                       outputType = gdal.GDT_Float32,
                       srcSRS = raster.meta.get('crs'),
                       dstSRS = 'EPSG:4617', # NAD83(CSRS)
                       dstNodata = raster.meta.get('nodata'),
                       creationOptions = ['COMPRESS=LZW'])
    gdal.Warp(output_dem, input_dem, options = wrap_option)
    
def reclassify_df(dataframe):
    ''' A function on the dataframe to resample raster by nearest neighbor method '''
    coords = [(lon,lat) for lon, lat in zip(dataframe.lon_c, dataframe.lat_c)]
    dataframe['vege_class'] = [DgBaseFunc.catch(lambda: vege[0]) for vege in VEGE_TIF.sample(coords)]
    return dataframe

if __name__=='__main__':
    
    input_csv_path = 'Result/Area_{}_centroids_{}.csv'.format(area,dggs_res)
    centroid_df = pandas.read_csv(input_csv_path, sep=',')

    input_csv_path = 'Result/{}_elev_{}.csv'.format(area,dggs_res)
    elev_df = pandas.read_csv(input_csv_path, sep=',')

    # inversely project vegetation class raster to NAD83 CSRS
    input_tif = 'Data/aci_2020_ab.tif'
    output_tif = 'Data/aci_2020_ab_NAD83.tif'
    reverse_projection(input_tif,output_tif)

    # record timing -- start
    start_time = time.time()

    # extract raster values by nearest method
    # call the function by parallel processing
    VEGE_TIF = rasterio.open('Data/aci_2020_ab_NAD83.tif')
    # centroid_df_output = reclassify_df(centroid_df)
    n_cores = int(os.environ.get('SLURM_CPUS_PER_TASK',default=1))
    centroid_df_split = numpy.array_split(centroid_df, n_cores)
    pool = mp.Pool(processes = n_cores)
    centroid_df_output = pandas.concat(pool.map(reclassify_df, centroid_df_split))
    pool.close()
    pool.join()

    # calculate stats according to zones
    elev_df = elev_df[elev_df['model_elev'] != -32767]
    centroid_df_join = pandas.merge(left = elev_df, right = centroid_df_output, how="inner", on=['i','j'])
    stats_df = centroid_df_join.groupby(["vege_class"]).agg(elev_mean = ('model_elev',numpy.nanmean), elev_max = ('model_elev',numpy.max), elev_min = ('model_elev',numpy.min),
                                                              elev_median = ('model_elev',numpy.nanmedian), elev_std = ('model_elev',numpy.nanstd), elev_sum = ('model_elev',numpy.sum))
    stats_df['elev_std'] = [0 if numpy.isnan(i) else i for i in stats_df['elev_std']]
    stats_df['elev_range'] = stats_df['elev_max'] - stats_df['elev_min']
    stats_df['volume'] = stats_df['elev_sum'] * cell_area / 1000000000

    # save the results as csv
    output_csv_path = 'Result/{}_elev_zonal_{}.csv'.format(area,dggs_res)
    stats_df.to_csv(output_csv_path, index=True)
    
