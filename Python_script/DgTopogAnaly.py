import pandas, math
import numpy as np
import time, sys, os, warnings, functools
import multiprocess as mp
import DgBaseFunc as dbfc

warnings.simplefilter('error', RuntimeWarning) 

# set resolution level
dggs_res = int(sys.argv[1])
area = sys.argv[2] # Calgary / Canmore / BuffaloLake
method = sys.argv[3] # MAG, MDG, MDN, FDA, BFP

# look up cell spacing and vertical resolution
look_up = dbfc.look_up_table()
cell_spacing = look_up.loc[dggs_res,'cell_spacing'] * 1000
vertical_res = look_up.loc[dggs_res,'verti_res']

# main functions
def slope_MAG(coords,df,res,cell_spacing):
    ''' 
    Absolute maximum differences between the center cell and its six neighbours
    If edge cell, assign nan values
    If flat, slope = 0; aspect = -1
    If not flat, slope angle is the normalized gradient, slope aspect is one of the restricted aspect values
    If multiple equal gradient found, then choose the first neighbour encountered clockwise from north
    If the slope is uphill, the downhill aspect is considered to be in the directly opposite direction
    '''
    elev_neighbor = dbfc.neighbor_navig(res,coords,df)
    if dbfc.edge_cell_exist(elev_neighbor):
        return [np.nan] * 2
    else:
        # gradient = [round(e-elev_neighbor[0],vertical_res) for e in elev_neighbor]
        gradient = [e-elev_neighbor[0] for e in elev_neighbor]
        gradient_abs = [abs(e) for e in gradient]
        gradient_max = max(gradient_abs)
        if gradient_max == 0:
            gradient_mag = 0
            aspect_mag = -1
        else:
            gradient_rad_mag = math.atan(gradient_max / cell_spacing)
            gradient_mag = math.degrees(gradient_rad_mag)
            aspect_index = gradient_abs.index(gradient_max)
            if gradient[aspect_index] > 0:
                aspect_mag = dbfc.aspect_restricted_oppo(res)[aspect_index-1]
            elif gradient[aspect_index] < 0:
                aspect_mag = dbfc.aspect_restricted(res)[aspect_index-1]
        return gradient_mag, aspect_mag

def slope_MDG(coords,df,res,cell_spacing):
    ''' 
    Maximum downward gradient between the center cell and its six neighbours
    Need pit-filling beforehand without saving the altered elevations
    If edge cell, assign nan values
    If flat, slope = 0; aspect = -1
    If not flat, slope angle is the normalized gradient, slope aspect is one of the restricted aspect values
    If multiple equal gradient found, then choose the first neighbour encountered clockwise from north
    '''
    elev_neighbor = dbfc.neighbor_navig(res,coords,df)
    if dbfc.edge_cell_exist(elev_neighbor):
        return [np.nan] * 2
    else:
        elev_neighbor_ = elev_neighbor[1:]
        if all(i > elev_neighbor[0] for i in elev_neighbor_):
            elev_neighbor[0] = min(elev_neighbor_)
        # gradient = [round(e-elev_neighbor[0],vertical_res) for e in elev_neighbor]
        gradient = [e-elev_neighbor[0] for e in elev_neighbor]
        gradient_min = min(gradient)
        if gradient_min >= 0:
            gradient_mdg = 0
            aspect_mdg = -1
        else:
            gradient_rad_mdg = math.atan(abs(gradient_min) / cell_spacing)
            gradient_mdg = math.degrees(gradient_rad_mdg)
            aspect_index = gradient.index(gradient_min)
            aspect_mdg = dbfc.aspect_restricted(res)[aspect_index-1]
        return gradient_mdg, aspect_mdg

def slope_MDN(coords,df,res,cell_spacing):
    ''' 
    Multiple downhill neighbours
    Need pit-filling beforehand without saving the altered elevations
    Distributing flow from a pixel amongst all of its lower elevation neighbor pixels
    If edge cell, assign nan values
    If flat, slope = 0; aspect = -1
    If not flat, slope and aspect are calculated by taking the average of all the downslope neighbors
    Aspect is represented by the mean norm vector which is a "mean" surface orientation
    '''
    elev_neighbor = dbfc.neighbor_navig(res,coords,df)
    if dbfc.edge_cell_exist(elev_neighbor):
        return [np.nan] * 2
    else:
        elev_neighbor_ = elev_neighbor[1:]
        if all(i > elev_neighbor[0] for i in elev_neighbor_):
            elev_neighbor[0] = min(elev_neighbor_)
        # gradient = [round(e-elev_neighbor[0],vertical_res) for e in elev_neighbor]
        gradient = [e-elev_neighbor[0] for e in elev_neighbor]
        if min(gradient) >= 0:
            gradient_mdn = 0
            aspect_mdn = -1
        else:
            gradient_rad_mdn_ls = [math.atan(abs(g) / cell_spacing) for g in gradient if g < 0]
            gradient_rad_mdn = sum(gradient_rad_mdn_ls) / len(gradient_rad_mdn_ls)
            gradient_mdn = math.degrees(gradient_rad_mdn)
            cell_ls = []
            i = 0
            for _ in gradient:
                if _ < 0:
                    cell_ls.append(i)
                i += 1
            if len(cell_ls) == 1:
                aspect_index = cell_ls[0]
                aspect_mdn = dbfc.aspect_restricted(res)[aspect_index-1]
            elif len(cell_ls) == 2 and abs(cell_ls[0]-cell_ls[1]) == 3:
                aspect_index_1, aspect_index_2 = cell_ls[0], cell_ls[1]
                gradient_1, gradient_2 = gradient[aspect_index_1], gradient[aspect_index_2]
                aspect_index = aspect_index_1 if gradient_1 <= gradient_2 else aspect_index_2
                aspect_mdn = dbfc.aspect_restricted(res)[aspect_index-1]
            else:
                avg_norm = dbfc.mean_norm_vector(res,elev_neighbor,*cell_ls)
                aspect_mdn = dbfc.aspect_unrestricted(avg_norm[0],avg_norm[1])
                aspect_mdn = math.degrees(aspect_mdn)
            if aspect_mdn < 0:
                aspect_mdn += 360
        return gradient_mdn, aspect_mdn

def slope_FDA(coords,df,res,cell_spacing):
    ''' 
    Finite-Difference Algorithm
    Project the non-normalized gradient in three directions to orthogonal x y axes
    If edge cell, assign nan values
    If flat, slope = 0; aspect = -1
    If not flat, calculate slope and aspect by using finite difference algorithms 
    and combine the two component partial derivatives
    '''
    elev_neighbor = dbfc.neighbor_navig(res,coords,df)
    if dbfc.edge_cell_exist(elev_neighbor):
        return [np.nan] * 2
    else:
        dzx, dzy = dbfc.first_derivative(res,elev_neighbor)
        if dzx == 0 and dzy == 0:
            gradient_fda = 0
            aspect_fda = -1
        else: 
            gradient_rad_fda = math.atan(math.sqrt(dzx**2 + dzy**2) / cell_spacing)
            gradient_fda = math.degrees(gradient_rad_fda)
            aspect_fda = dbfc.aspect_unrestricted(dzx,dzy)
            aspect_fda = math.degrees(aspect_fda)
        if aspect_fda < 0:
            aspect_fda += 360
        return gradient_fda, aspect_fda

def slope_BFP(coords,df,res,cell_spacing):
    '''
    Best fit plane method
    If edge cell, assign nan values
    If flat, slope = 0; aspect = -1
    If not flat, a surface is fitted to seven centroids by multiple linear regression models, 
    using least squares to minimize the sum of distances from the surface to the cells
    '''
    elev_neighbor = dbfc.neighbor_navig(res,coords,df)
    if dbfc.edge_cell_exist(elev_neighbor):
        return [np.nan] * 2
    else:
        norm_vec = np.array(dbfc.fit_plane_norm_vector(res,elev_neighbor))
        norm_vec_mag = np.linalg.norm(norm_vec)
        unit_norm_vec = np.array([i/norm_vec_mag for i in norm_vec])
        ref_vec = np.array([0,0,-1])
        ref_vec_proj = ref_vec - (np.dot(ref_vec, unit_norm_vec) * unit_norm_vec)
        try:
            gradient_rad_bfp = math.atan((ref_vec_proj[2]**2 / math.sqrt(ref_vec_proj[0]**2 + ref_vec_proj[1]**2)) / cell_spacing)
            gradient_bfp = math.degrees(gradient_rad_bfp)
            aspect_bfp = dbfc.aspect_unrestricted(ref_vec_proj[0],ref_vec_proj[1])
            aspect_bfp = math.degrees(aspect_bfp)
        except:
            gradient_bfp = 0
            aspect_bfp = -1
        if aspect_bfp < 0:
            aspect_bfp += 360
        return gradient_bfp, aspect_bfp

def curvature(coords,df,res,cell_spacing):
    '''
    slope rate of change of landfrom
    second dervivative of DEM
    first derivitave of slope
    '''
    elev_neighbor = dbfc.neighbor_navig(res,coords,df)
    if dbfc.edge_cell_exist(elev_neighbor):
        return np.nan
    else:
        dzx2, dzy2 = dbfc.second_derivative(res,elev_neighbor)
        curv = math.sqrt(dzx2**2 + dzy2**2) / cell_spacing**2
        return curv
  
def hillshade(coords,df,res,altitude=45,azimuth=315):
    ''' 
    A hillshade is a grayscale 3D representation of the surface, with the sun's relative position taken into account for shading. 
    This function uses the altitude and azimuth properties to specify the sun's position.
    Azimuth is the angular direction of the sun, measured from north in clockwise degrees from 0 to 360. 
    An azimuth of 90 degrees is east. The default azimuth is 315 degrees (NW).
    Altitude is the slope or angle of the illumination source above the horizon. 
    The units are in degrees, from 0 (on the horizon) to 90 (overhead). The default is 45 degrees.
    '''
    slope_deg,aspect_deg = slope_FDA(coords,df)
    slope_rad = math.radians(slope_deg)
    aspect_rad = math.radians(aspect_deg)
    if slope_rad == 0:
        hs = 255.0
    else:
        zenith_deg = 90.0 - altitude
        zenith_rad = math.radians(zenith_deg)
        azimuth_math = 360.0 - azimuth + 90.0
        if azimuth_math >= 360.0:
            azimuth_math = azimuth_math - 360.0
        azimuth_rad = math.radians(azimuth_math)
        hs = 255.0 * ((math.cos(zenith_rad) * math.cos(slope_rad)) + 
                       (math.sin(zenith_rad) * math.sin(slope_rad) * math.cos(azimuth_rad - aspect_rad)))
    return hs

def slope_aspect_df(dataframe,elev_df,method,res,cell_spacing):
    ''' calculate slope and aspect by specified method '''
    if method == 'MAG':
        dataframe[['gradient_deg','aspect_deg']] = [slope_MAG(ij,elev_df,res,cell_spacing) for ij in dataframe.index.values]
    elif method == 'MDG':
        dataframe[['gradient_deg','aspect_deg']] = [slope_MDG(ij,elev_df,res,cell_spacing) for ij in dataframe.index.values]
    elif method == 'MDN':
        dataframe[['gradient_deg','aspect_deg']] = [slope_MDN(ij,elev_df,res,cell_spacing) for ij in dataframe.index.values]
    elif method == 'FDA':
        dataframe[['gradient_deg','aspect_deg']] = [slope_FDA(ij,elev_df,res,cell_spacing) for ij in dataframe.index.values]
    elif method == 'BFP':
        dataframe[['gradient_deg','aspect_deg']] = [slope_BFP(ij,elev_df,res,cell_spacing) for ij in dataframe.index.values]
    return dataframe

def curvature_df(dataframe,elev_df,res,cell_spacing):
    dataframe['curv'] = [curvature(ij,elev_df,res,cell_spacing) for ij in dataframe.index.values]
    return dataframe

def hillshade_df(dataframe,elev_df,res):
    dataframe['hs'] = [hillshade(ij,res,elev_df) for ij in dataframe.index.values]
    return dataframe

## call fuctions in parallel computing mode and record timing
if __name__=='__main__':
    
    input_csv_path = 'Result/{}_elev_{}.csv'.format(area,dggs_res)
    elev_df = pandas.read_csv(input_csv_path, sep=',')
    elev_df = elev_df.set_index(['i', 'j'])

    ############################################################

    # record timing -- start
    start_time = time.time()

    elev_df_copy = elev_df.copy()

    # call the function by parallel processing
    slope_aspect_df_p = functools.partial(slope_aspect_df,elev_df=elev_df,method=method,res=dggs_res,cell_spacing=cell_spacing)

    n_cores = int(os.environ.get('SLURM_CPUS_PER_TASK',default=1))
    elev_df_split = np.array_split(elev_df_copy, n_cores)
    pool = mp.Pool(processes = n_cores)
    elev_df_output = pandas.concat(pool.map(slope_aspect_df_p, elev_df_split))
    pool.close()
    pool.join()

    # record timing -- end
    print (dggs_res)
    print (method)
    print ("Processing time: %s seconds" % (time.time() - start_time))

    # save the results as csv
    output_csv_path = 'Result/{}_elev_{}_{}.csv'.format(area,method,dggs_res)
    elev_df_output.to_csv(output_csv_path, index=True)

    ############################################################
    ## curvature
    # record timing -- start
    start_time = time.time()

    elev_df_copy = elev_df.copy()

    # call the function by parallel processing
    curvature_df_p = functools.partial(curvature_df,elev_df=elev_df,res=dggs_res,cell_spacing=cell_spacing)

    n_cores = int(os.environ.get('SLURM_CPUS_PER_TASK',default=1))
    elev_df_split = np.array_split(elev_df_copy, n_cores)
    pool = mp.Pool(processes = n_cores)
    elev_df_output = pandas.concat(pool.map(curvature_df_p, elev_df_split))
    pool.close()
    pool.join()

    # record timing -- end
    print (dggs_res)
    print ("Processing time: %s seconds" % (time.time() - start_time))

    # save the results as csv
    output_csv_path = 'Result/{}_curvature_{}.csv'.format(area,dggs_res)
    elev_df_output.to_csv(output_csv_path, index=True)

    ############################################################
    ## hillshade
    # record timing -- start
    start_time = time.time()

    elev_df_copy = elev_df.copy()

    # call the function by parallel processing
    hillshade_df_p = functools.partial(hillshade_df,elev_df=elev_df,res=dggs_res)

    n_cores = int(os.environ.get('SLURM_CPUS_PER_TASK',default=1))
    elev_df_split = np.array_split(elev_df_copy, n_cores)
    pool = mp.Pool(processes = n_cores)
    elev_df_output = pandas.concat(pool.map(hillshade_df_p, elev_df_split))
    pool.close()
    pool.join()

    # record timing -- end
    print (dggs_res)
    print ("Processing time: %s seconds" % (time.time() - start_time))

    # save the results as csv
    output_csv_path = 'Result/{}_hillshade_{}.csv'.format(area,dggs_res)
    elev_df_output.to_csv(output_csv_path, index=True)
