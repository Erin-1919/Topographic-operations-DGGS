import pandas, numpy, statistics
import time, sys, os
import multiprocess as mp
import DgBaseFunc
from itertools import product

# set resolution level, study area, and ring number
dggs_res = int(sys.argv[1]) # 20-24
ring_n = int(sys.argv[2]) # suggest 1-3
area = sys.argv[3] # Calgary / Canmore / BuffaloLake

# look up cellsize and vertical resolution
look_up = DgBaseFunc.look_up_table()
dggs_cellsize = look_up.loc[dggs_res,'cell_size']
vertical_res = look_up.loc[dggs_res,'verti_res']

## Define the main funcitons
def focal_elev_stats(coords,df,rings=1,res=dggs_res):
    '''
    remove voids in neighbors
    calculate descriptive stats -- 
    mean, max, min, median, std, range
    '''
    elev_neighbor = DgBaseFunc.neighbor_navig_by_rings(res,coords,df,rings)
    elev_neighbor = [x for x in elev_neighbor if x != -32767 and not numpy.isnan(x)]
    if len(elev_neighbor) == 0:
        return [numpy.nan] * 6
    if len(elev_neighbor) == 1:
        return elev_neighbor*4+[0]*2
    else: 
        elev_mean = round(statistics.mean(elev_neighbor),vertical_res)
        elev_max = round(max(elev_neighbor),vertical_res)
        elev_min = round(min(elev_neighbor),vertical_res)
        elev_median = round(statistics.median(elev_neighbor),vertical_res)
        elev_std = round(statistics.stdev(elev_neighbor),vertical_res)
        elev_range = round(elev_max-elev_min,vertical_res)
        return [elev_mean,elev_max,elev_min,elev_median,elev_std,elev_range]
        
def focal_elev_stats_df(dataframe,rings=1):
    dataframe[['mean','max','min','median','std','range']] = [focal_elev_stats(ij,elev_df,rings) for ij in dataframe.index.values]
    return dataframe

if __name__=='__main__':
    
    ## read the csv into a dataframe
    input_csv_path = 'Result/{}_elev_{}.csv'.format(area,dggs_res)
    elev_df = pandas.read_csv(input_csv_path, sep=',')
    elev_df = elev_df.set_index(['i', 'j'])

    # record timing -- start
    start_time = time.time()

    # call the function by parallel processing
    # need to initialize columns with python 3.7 on linux
    n_cores = int(os.environ.get('SLURM_CPUS_PER_TASK',default=1))
    elev_df_split = numpy.array_split(elev_df, n_cores)
    pool = mp.Pool(processes = n_cores)
    elev_df_output = pandas.concat(pool.starmap(focal_elev_stats_df, product(elev_df_split,[ring_n]*len(elev_df_split))))
    pool.close()
    pool.join()

    # record timing -- end
    print (dggs_res)
    print ("Processing time: %s seconds" % (time.time() - start_time))

    # save the results as csv
    output_csv_path = 'Result/{}_elev_focal_{}_ring{}.csv'.format(area,dggs_res,ring_n)
    elev_df_output.to_csv(output_csv_path, index=True)
