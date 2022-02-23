import pandas, numpy, math
import time, sys, os
import multiprocess as mp
import DgBaseFunc

# set resolution level and study area
dggs_res = int(sys.argv[1]) # 20-24
area = sys.argv[2] # Calgary / Canmore / BuffaloLake
method = sys.argv[3] # MAG, MDG, MDN, FDA, BFP

# look up cell spacing and vertical resolution
look_up = DgBaseFunc.look_up_table()
cell_spacing = look_up.loc[dggs_res,'cell_spacing'] * 1000
vertical_res = look_up.loc[dggs_res,'verti_res']

## Define the main funcitons
def TRI_calcu(coords,res=dggs_res):
    ''' Terrain roughness index (TRI) 
    TRI presents the square root of sum of squared difference between a central cell and the adjacent cells'''
    elev_neighbor = DgBaseFunc.neighbor_navig(res,coords,elev_df)
    if DgBaseFunc.edge_cell_exist(elev_neighbor):
        TRI_value = numpy.nan
    else:
        TRI_ls = [(n-elev_neighbor[0]) ** 2 for n in elev_neighbor[1:]]
        TRI_value = math.sqrt(sum(TRI_ls))
    return TRI_value
    
def TPI_calcu(coords,res=dggs_res):
    ''' Topographic position index (TPI) 
    TPI presents the difference between a central cell and the mean of its surrounding cells'''
    elev_neighbor = DgBaseFunc.neighbor_navig(res,coords,elev_df)
    if DgBaseFunc.edge_cell_exist(elev_neighbor):
        TPI_value = numpy.nan
    else:
        TPI_ls = elev_neighbor[1:]
        TPI_value = elev_neighbor[0] - (sum(TPI_ls) / len(TPI_ls))
    return TPI_value
    
def TWI_calcu(alpha,beta):
    ''' Topographic wetness index (TWI) '''
    if beta == 0:
        TWI_value = numpy.nan
    else:
        beta_rad = math.radians(beta)
        TWI_value = math.log(alpha/math.tan(beta_rad))
    return TWI_value
    
def SPI_calcu(alpha,beta):
    ''' Stream power index (SPI) '''
    beta_rad = math.radians(beta)
    SPI_value = alpha * math.tan(beta_rad)
    return SPI_value

def TopoIndex_df(dataframe):
    dataframe['TRI'] = [TRI_calcu(ij) for ij in dataframe.index.values]
    dataframe['TPI'] = [TPI_calcu(ij) for ij in dataframe.index.values]
    dataframe['TWI'] = [TWI_calcu(a,b) for a,b in zip(dataframe['contri_area'],dataframe['gradient_deg'])]
    dataframe['SPI'] = [SPI_calcu(a,b) for a,b in zip(dataframe['contri_area'],dataframe['gradient_deg'])]
    return dataframe

#############################################################################

print (dggs_res)
print (area)
print (method)

## read the csv into a dataframe
input_csv_path = 'Result/{}_flow_{}_{}.csv'.format(area,method,dggs_res)
elev_df = pandas.read_csv(input_csv_path, sep=',')
elev_df = elev_df.set_index(['i', 'j'])

## call fuctions in parallel computing mode and record timing

# record timing -- start
start_time = time.time()

elev_df['TRI'] = elev_df['TPI'] = elev_df['TWI'] = elev_df['SPI'] = numpy.nan

# call the function by parallel processing
n_cores = int(os.environ.get('SLURM_CPUS_PER_TASK',default=1))
elev_df_split = numpy.array_split(elev_df, n_cores)
pool = mp.Pool(processes = n_cores)
elev_df_output = pandas.concat(pool.map(TopoIndex_df, elev_df_split))
pool.close()
pool.join()

# record timing -- end
print ("Processing time: %s seconds" % (time.time() - start_time))

# save the results as csv
output_csv_path = 'Result/{}_TopoIndex_{}_{}.csv'.format(area,method,dggs_res)
elev_df_output.to_csv(output_csv_path, index=True)

