import pandas, numpy, scipy, math
import scipy.linalg

def look_up_table():
    ''' construct a look-up table, storing DGGS resolution, corresponding cell size, cell area, and vertical resolution '''
    res_list = [16,17,18,19,20,21,22,23,24,25,26,27,28,29]
    cell_size_list = [0.005,0.003,0.001,0.0009,0.0008,0.0003,0.0003,0.0001,0.0001,0.00006,0.00003,0.00002,0.00001,0.000005]
    cell_spacing_list = [1.0750880097,0.6207023518,0.3583626699,0.2069007839,0.1194542233,0.0689669280,0.0398180744,
                     0.0229889760,0.0132726915,0.0076629920,0.0044242305,0.0025543307,0.0014747435,0.0008514436]
    cell_area_list = [1.1849116724224,0.3949705574741,0.1316568524914,0.0438856174971,0.0146285391657,0.0048761797219,
                  0.0016253932406,0.0005417977469,0.0001805992490,0.0000601997497,0.0000200665832,0.0000066888611,
                  0.0000022296204,0.0000007432068]
    vertical_res_list = [0,0,0,1,1,2,2,3,3,4,4,5,5,6]
    look_up = pandas.DataFrame({'res': res_list, 'cell_size': cell_size_list, 'verti_res': vertical_res_list, 
                                'cell_area': cell_area_list, 'cell_spacing': cell_spacing_list}, index = res_list)
    return look_up

def checkIfDuplicates(listOfElems):
    ''' Check if given list contains any duplicates '''    
    setOfElems = set()
    for elem in listOfElems:
        if elem in setOfElems:
            return True
        else:
            setOfElems.add(elem)         
    return False

def common_member(a, b):
    ''' check if two lists have at-least one element common '''
    a_set = set(a)
    b_set = set(b)
    return (a_set & b_set)
    
def catch(func, handle = lambda e : e, *args, **kwargs):
    ''' Handle the try except in a general function'''
    try:
        return func(*args, **kwargs)
    except:
        return numpy.nan

def edge_cell_exist(ls):
    ''' determine if edge cell exists in the neighborhood '''
    cell_ls = [x for x in ls if x != -32767 and not numpy.isnan(x)]
    return len(cell_ls) != 7

def edge_cell(res,coords,allcell):
    neighbor = neighbor_coords(res,coords)[1:]
    neighbor_in = [i in allcell for i in neighbor]
    neighbor_filt = list(filter(lambda x: x is True,neighbor_in))
    return len(neighbor_filt) != 6
    
def edge_cell_df(df,res,allcell):
    ls = []
    for ij in df.index.values:
        if edge_cell(res,ij,allcell):
            ls.append(ij)
    df_out = df[df.index.isin(ls)]
    return df_out

def edge_cell_chain(res,coords,allcell):
    chain = [coords]
    nxt = numpy.nan
    while nxt != coords:
        c = chain[-1]
        neighbor = neighbor_coords(res,c)[1:]
        neighbor_in = [i for i in neighbor if i in allcell]
        for n in neighbor_in:
            if edge_cell(res,n,allcell) and (n not in chain[1:]):
                nxt = n
                if nxt != coords:
                    chain.append(nxt)
                break
    return chain

def neighbor_coords(res,coords):
    ''' find out the ij coordinates within the neighborhood -- 6 neighbors + center cell '''
    i,j = coords[0],coords[1]
    if res%2 == 0:
        index_i = [i,i-1,i-1,i,i+1,i+1,i]
        index_j = [j,j-1,j,j+1,j+1,j,j-1]
    elif res%2 == 1:
        index_i = [i,i-2,i-1,i+1,i+2,i+1,i-1]
        index_j = [j,j-1,j+1,j+2,j+1,j-1,j-2]
    coords_neighbor = [(ih,jh) for ih,jh in zip(index_i,index_j)]
    return coords_neighbor

def neighbor_navig_by_rings(res,coords,df,rings=1):
    ''' navigate among neighborhood defined by ring numbers (default 1), 
    return elevation values of neighbors + center cell '''
    i = 1
    coords_edge = [coords]
    coords_neighbor = [coords]
    while i <= rings:
        i += 1
        coords_neighbor_ = coords_neighbor
        for c in coords_edge:
            coords_neighbor = coords_neighbor + neighbor_coords(res,c)[1:]
        coords_neighbor = list(set(coords_neighbor))
        coords_edge = [x for x in coords_neighbor if x not in coords_neighbor_]
    elev_neighbor = [catch(lambda: df['model_elev'].loc[ij]) for ij in coords_neighbor]
    return elev_neighbor

def edge_navig_by_rings(res,coords,rings=1):
    ''' navigate among neighborhood defined by ring numbers (default 1), 
    return elevation values of neighbors + center cell '''
    i = 1
    coords_edge = [coords]
    coords_neighbor = [coords]
    while i <= rings:
        i += 1
        coords_neighbor_ = coords_neighbor
        for c in coords_edge:
            coords_neighbor = coords_neighbor + neighbor_coords(res,c)[1:]
        coords_neighbor = list(set(coords_neighbor))
        coords_edge = [x for x in coords_neighbor if x not in coords_neighbor_]
    return coords_edge

def hex_dist(coord1,coord2,res):
    ''' calculate hex distance represented by ring number '''
    i_ = coord1[0] - coord2[0]
    j_ = coord1[1] - coord2[1]
    if i_ * j_ >= 0:
        ring = max(abs(i_),abs(j_))
    else:
        ring = abs(i_) + abs(j_)
    if res%2 == 1:
        # ring = math.ceil(ring/2)
        ring = round(ring/numpy.sqrt(3))
    return ring

def hex_dist_l(coord,coords_target,res):
    ''' calculate hex distance from one coord to a list of target coords (if there are a few), 
    represented by ring number '''
    dist_ls = [hex_dist(coord,coord_,res) for coord_ in coords_target]
    min_dist = min(dist_ls)
    return min_dist

def hex_dist_m(coord,coords_target,res,thred=30):
    ''' calculate hex distance from one coord to a list of target coords (if there are a lot), 
    represented by ring number '''
    i = 0
    coords_edge = [coord]
    coords_neighbor = [coord]
    while not common_member(coords_target, coords_edge):
        i += 1
        coords_neighbor_ = coords_neighbor
        for c in coords_edge:
            coords_neighbor = coords_neighbor + neighbor_coords(res,c)[1:]
        coords_neighbor = list(set(coords_neighbor))
        coords_edge = [x for x in coords_neighbor if x not in coords_neighbor_]
        if i > thred:
            i = numpy.nan
            break
    return i

def neighbor_navig(res,coords,df):
    ''' extract elevation values among the neighborhood -- 1 center + 6 neighbors '''
    coords_neighbor = neighbor_coords(res,coords)
    elev_neighbor = [catch(lambda: df['model_elev'].loc[ij]) for ij in coords_neighbor]
    return elev_neighbor

def neighbor_slope_navig(res,coords,df):
    ''' extract slope values among the neighborhood -- 1 center + 6 neighbors '''
    coords_neighbor = neighbor_coords(res,coords)
    slp_neighbor = [catch(lambda: df['gradient_deg'].loc[ij]) for ij in coords_neighbor]
    return slp_neighbor

def neighbor_direc_navig(res,coords,df):
    ''' extract direction codes among the neighborhood -- 1 center + 6 neighbors '''
    coords_neighbor = neighbor_coords(res,coords)
    direc_neighbor = [catch(lambda: df['direction_code'].loc[ij]) for ij in coords_neighbor]
    return direc_neighbor

def neighbor_multi_direc_navig(res,coords,df):
    ''' extract direction codes among the neighborhood -- 1 center + 6 neighbors '''
    coords_neighbor = neighbor_coords(res,coords)
    direc_multi_neighbor = [catch(lambda: df['direction_code'].loc[ij]) for ij in coords_neighbor]
    return direc_multi_neighbor

def cross_product(p0,p1,p2):
    ''' calculate the cross product of two vectors composed by three points '''
    x0, y0, z0 = p0
    x1, y1, z1 = p1
    x2, y2, z2 = p2
    ux, uy, uz = [x1-x0, y1-y0, z1-z0] # vector u
    vx, vy, vz = [x2-x0, y2-y0, z2-z0] # vector v
    u_cross_v = numpy.array([uy*vz-uz*vy, uz*vx-ux*vz, ux*vy-uy*vx])
    if u_cross_v[2] < 0:
        u_cross_v = u_cross_v * (-1)
    return u_cross_v

def aspect_restricted(res):
    ''' restricted aspect lists '''
    if res%2 == 1:
        asp_list = [30,90,150,210,270,330]
    elif res%2 == 0:
        asp_list = [0,60,120,180,240,300]   
    return asp_list

def aspect_restricted_oppo(res):
    ''' restricted aspect in opposite directions '''
    if res%2 == 1:
        asp_list = [210,270,330,30,90,150]
    elif res%2 == 0:
        asp_list = [180,240,300,0,60,120]
    return asp_list

def aspect_unrestricted(dzx,dzy):
    ''' the unrestricted direction to which that cell is oriented; algorithms from Hodgson 1998'''
    if dzx == 0:
        if dzy < 0:
            asp = math.pi
        elif dzy > 0:
            asp = 0
        else:
            asp = -1
    elif dzx > 0:
        asp = math.pi/2 - math.atan(dzy/dzx)
    else:
        asp = math.pi/2*3 - math.atan(dzy/dzx)
    return asp

def mean_norm_vector(res,z_ls, *cells):
    ''' calculate the mean of norm vectors, vector count >1 '''
    if res%2 == 1:
        x_ls = [0,math.sqrt(3)/3,math.sqrt(3)/3*2,math.sqrt(3)/3,math.sqrt(3)/-3,math.sqrt(3)/-3*2,math.sqrt(3)/-3]
        y_ls = [0,1,0,-1,-1,0,1]
    elif res%2 == 0:
        x_ls = [0,0,1,1,0,-1,-1]
        y_ls = [0,math.sqrt(3)/3*2,math.sqrt(3)/3,math.sqrt(3)/-3,math.sqrt(3)/-3*2,math.sqrt(3)/-3,math.sqrt(3)/3]
    pt_c = [x_ls[0],y_ls[0],z_ls[0]]
    points = []
    x_sum = 0
    y_sum = 0
    z_sum = 0
    for i in cells:
        points.append([x_ls[i],y_ls[i],z_ls[i]])
    for pt in points:
        if points.index(pt) == len(points) - 1:
            break
        norm = cross_product(pt_c,pt,points[points.index(pt)+1])
        x_sum = x_sum + norm[0]
        y_sum = y_sum + norm[1]
        z_sum = z_sum + norm[2]
    return numpy.array([x_sum,y_sum,z_sum])

def fit_plane_norm_vector(res,z_ls):
    ''' fit a plane by using least squares and find out its norm vector '''
    if res%2 == 1:
        x_ls = [0,math.sqrt(3)/3,math.sqrt(3)/3*2,math.sqrt(3)/3,math.sqrt(3)/-3,math.sqrt(3)/-3*2,math.sqrt(3)/-3]
        y_ls = [0,1,0,-1,-1,0,1]
    elif res%2 == 0:
        x_ls = [0,0,1,1,0,-1,-1]
        y_ls = [0,math.sqrt(3)/3*2,math.sqrt(3)/3,math.sqrt(3)/-3,math.sqrt(3)/-3*2,math.sqrt(3)/-3,math.sqrt(3)/3]
    G = numpy.c_[x_ls, y_ls, z_ls]
    A = numpy.c_[G[:,0], G[:,1], numpy.ones(G.shape[0])]
    C,_,_,_ = scipy.linalg.lstsq(A, G[:,2])
    norm_vector = numpy.array([C[0]*-1,C[1]*-1,1])
    return norm_vector

def first_derivative(res,elev_neighbor):
    ''' project the non-normalized gradient to orthogonal x y axes '''
    if res%2 == 1:
        d_zi = (elev_neighbor[3] - elev_neighbor[6]) / 2
        d_zj = (elev_neighbor[4] - elev_neighbor[1]) / 2
        d_zk = (elev_neighbor[5] - elev_neighbor[2]) / 2
        d_zx = d_zk + d_zj * math.sin(math.pi/6) - d_zi * math.sin(math.pi/6)
        d_zy = d_zi * math.cos(math.pi/6) + d_zj * math.cos(math.pi/6)
    elif res%2 == 0:
        d_zi = (elev_neighbor[4] - elev_neighbor[1]) / 2
        d_zj = (elev_neighbor[5] - elev_neighbor[2]) / 2
        d_zk = (elev_neighbor[6] - elev_neighbor[3]) / 2
        d_zx = d_zj * math.cos(math.pi/6) + d_zk * math.cos(math.pi/6)
        d_zy = d_zi + d_zj * math.sin(math.pi/6) - d_zk * math.sin(math.pi/6)
    return d_zx, d_zy
            
def second_derivative(res,elev_neighbor):
    ''' project the non-normalized second derivatives to orthogonal x y axes '''          
    if res%2 == 1:
        d_zi2 = 2 * elev_neighbor[0] - elev_neighbor[3] - elev_neighbor[6]
        d_zj2 = 2 * elev_neighbor[0] - elev_neighbor[1] - elev_neighbor[4]
        d_zk2 = 2 * elev_neighbor[0] - elev_neighbor[2] - elev_neighbor[5]
        d_zx2 = d_zk2 + d_zj2 * math.sin(math.pi/6) - d_zi2 * math.sin(math.pi/6)
        d_zy2 = d_zi2 * math.cos(math.pi/6) + d_zj2 * math.cos(math.pi/6)
    elif res%2 == 0:
        d_zi2 = 2 * elev_neighbor[0] - elev_neighbor[1] - elev_neighbor[4]
        d_zj2 = 2 * elev_neighbor[0] - elev_neighbor[2] - elev_neighbor[5]
        d_zk2 = 2 * elev_neighbor[0] - elev_neighbor[3] - elev_neighbor[6]
        d_zx2 = d_zj2 * math.cos(math.pi/6) + d_zk2 * math.cos(math.pi/6)
        d_zy2 = d_zi2 + d_zj2 * math.sin(math.pi/6) - d_zk2 * math.sin(math.pi/6)
    return d_zx2, d_zy2

def hex_dist_comb_df(dataframe,coords_target,var,res):
    dataframe[var] = numpy.nan
    dataframe[var] = [hex_dist_m(coord,coords_target,res) for coord in dataframe.index.values]
    coords_left = dataframe[numpy.isnan(dataframe[var])].index.values.tolist()
    dist_left = [hex_dist_l(coord,coords_target,res) for coord in coords_left]
    for c,d in zip(coords_left,dist_left):
        dataframe.loc[c,var] = d
        
if __name__=='__main__':
    pass
