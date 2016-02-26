# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 12:34:43 2016
@author: Evan Coopersmith

Implements the Height-Above-Nearest-Drainage (HAND) model (Nobre et al, 2011 via Renno et al, 2008)
"""

"""     |     |    
     8  |  1  |  2
        |     |
  ------------------       
        |     |    
     7  |  0  |  3
        |     |
  ------------------       
        |     |    
     6  |  5  |  4
        |     |      
        
Classification of local drain directions...N = 1, NE = 2, W = 3, ..., NW = 8.  A 'pit' = 0"""

import json
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors

_path_to_DEM = os.path.join("ARS_SM","WGW","Results","Texture_and_Topography", "ProcessedJSONs")
_file_name = 'wg2_elevation'
_path_to_HAND_output = os.path.join(_path_to_DEM, "HAND")
_json_params = ['ncols', 'nrows', 'xllcorner', 'yllcorner', 'cellsize', 'NODATA_value']

_TEST_JSON = {"xllcorner": 100, "yllcorner": 1234, "cellsize": 10, 
              'ncols':7, 'nrows':7, 'NoDATA_value':9999, 
              '1294':[36,35,36,39,42,48,50], '1284':[39,32,33,35,38,46,51], '1274':[47,40,36,35,38,43,46],
              '1264':[52,48,39,36,34,35,34], '1254':[52,50,47,38,34,33,30], '1244':[51,50,49,45,40,35,33],
              '1234':[52,51,52,50,47,41,37]}

_flow_dirs = {1:(0,-1), 2:(1,-1), 3:(1,0), 4:(1,1), 5:(0,1), 6:(-1,1), 7:(-1,0), 8:(-1,-1)}
_inv_flow_dirs = {(0,-1):1, (1,-1):2, (1,0):3, (1,1):4, (0,1):5, (-1,1):6, (-1,0):7, (-1,-1):8}

_max_A_iterations = 3 #how far out should we look to assess upstream area?
_max_sink_search_dist = 10 #how many squares out should we look before considering an INTERNAL point a sink?
_min_drainage = 100 #how much area must be drained before we consider it to be on a 'flowpath'
_HAND_thresholds = [0.01,0.03,0.09,0.27]
_cmap = matplotlib.colors.ListedColormap(['darkblue', 'blue','lawngreen','yellow','red']) #to visualize

_DEM_cellsize, _western_edge, _southern_edge = 10, 580000, 3504000 #if we are not building from scratch...this is needed information
_unit = 'm' #chosen units of elevation...mostly necessary for chart legends

def ConvertDEM_json_to_np_array(json_params, feature_json, min_northing, nrows, ncols, cellsize):
    max_northing = min_northing + cellsize*(nrows-1)
    return np.array([feature_json[str(max_northing - r*cellsize)] for r in range(nrows)])
    
def DefineLDD(DEM, nrows, ncols, flow_dirs, inv_flow_dirs):
    LDD = np.zeros(shape=(nrows,ncols))
    for row in range(nrows):
        for col in range(ncols):
            current_elev, lowest_elev = DEM[row,col], 99999
            N_i_j = GetN_i_j(row, col, nrows, ncols)  
            for i,(next_row,next_col) in enumerate(N_i_j):
                row_diff, col_diff = next_row-row, next_col-col#; print(next_row,next_col) 
                possible_flow = inv_flow_dirs[(col_diff,row_diff)]#; print(possible_flow)
                next_elev = DEM[next_row, next_col]#; print(next_elev)
                if next_elev < lowest_elev and next_elev < current_elev and LDD[next_row, next_col] != inv_flow_dirs[(row_diff*-1,col_diff*-1)]:
                    lowest_elev = next_elev; LDD[row,col] = possible_flow#; print('Best')
    return LDD

def GetF_i_j(flow_dirs, flow_dir, row, col):
    """Given a (row),(col) within the DEM, the (flow_dir), which is a number from 1-north to 8-northwest 
    working clockwise, and a dictionary (flow_dirs) describing how to locate the point to which (row,col)
    flows, return the outlet_row and outlet_col for that new point..."""
    return ((row + flow_dirs[flow_dir][1], col + flow_dirs[flow_dir][0]) if flow_dir > 0 else (-1,-1)) #if it's a pit, it doesn't flow into anywhere     

def GetF_inv_i_j(LDD, flow_dirs, row, col, N_i_j):
    """Given a (row),(col), a numpy array defining adjusted flow directions at each point (LDD), the size of 
    the DEM in question (nrows),(ncol), and the directions associated with (flow_dirs), gather all neighboring
    squares (N_i_j) that flow into row,col"""
    return [p for p in N_i_j if GetF_i_j(flow_dirs, LDD[p[0],p[1]], p[0], p[1]) == (row, col)]
    
def GetN_i_j(row, col, nrows, ncols):
    neighboring_points = [(row+1,col),(row+1,col+1),(row,col+1),(row-1,col+1),
                          (row-1,col),(row-1,col-1),(row,col-1),(row+1,col-1)]
    return [p for p in neighboring_points if IsInBox(p[0], p[1], nrows, ncols)]
    
def Generate_all_A_i_j(nrows, ncols, LDD, flow_dirs, max_iters):
    drainage_areas = np.zeros(shape=(nrows, ncols))
    for row in range(nrows):
        for col in range(ncols):
            drainage_areas[row,col] = GetA_i_j(LDD, flow_dirs, drainage_areas, row, col, nrows, ncols, max_iters)
    return drainage_areas
        
def GetA_i_j(LDD, flow_dirs, drainage_areas, row, col, nrows, ncols, max_iters, iter_count = 0):
    """This represents the approximate number of squares within the DEM that drain to the point
    (row),(col), provided we do not look farther uphill than (max_iters) squares.  (LDD) contains
    the direction in which all squares flow, with sinks resolved."""
    known_drainage = drainage_areas[row,col]
    if known_drainage > 0: #we've already done the recursion in full at this point:
        return known_drainage
    else:
        N_i_j = GetN_i_j(row, col, nrows, ncols)
        F_inv_i_j = GetF_inv_i_j(LDD, flow_dirs, row, col, N_i_j) #the points that flow into (row,col)   
        if len(F_inv_i_j) == 0: #if there aren't any points that flow into this point
            drainage_areas[row,col] = 1; return 1
        elif iter_count < max_iters:
            iter_count += 1
            A_i_j = 1 + np.sum([GetA_i_j(LDD, flow_dirs, drainage_areas, p[0], p[1], nrows, ncols, max_iters, iter_count) for p in F_inv_i_j])
            return A_i_j        
        else: #maximum recursion depth count, STOP
            return 1       
            
def Generate_all_I_i_j(nrows, ncols, LDD, flow_dirs, drainage_network, DEM, max_sink_search_dist):
    nearest_drainage = np.zeros(shape=(nrows, ncols))
    for row in range(nrows):
        print("Following paths from row", row, "to their drains...")
        for col in range(ncols):
            nearest_drainage[row,col] = GetI_i_j(LDD, flow_dirs, row, col, nrows, ncols, 
                            drainage_network, nearest_drainage, DEM, max_sink_search_dist)            
    return nearest_drainage

def GetI_i_j(LDD, flow_dirs, row, col, nrows, ncols, drainage_network, nearest_drainage, DEM, max_dist):
    """For a point (row),(col), return the point in the (drainage_network) to which it will ultimately
    flow, using the local-drain-direction (LDD) array."""    
    drain_num = nearest_drainage[row,col]
    drain_path = [] #to ensure we don't return to the same place
    while drain_num == 0:
        if (row,col) in drainage_network.keys():    
            drain_num = drainage_network[(row,col)]
        else:
            drain_path.append((row,col))
            row, col = GetF_i_j(flow_dirs, LDD[row,col], row, col)
            if (row,col) in drain_path[:-1]: #we've been here before...
                view_dist, current_elev = 2, DEM[row,col]
                while view_dist <= max_dist:
                    coord_list = GetRingAtDistance(view_dist, row, col, nrows, ncols)
                    elevs = [DEM[c] for c in coord_list]; min_elev = min(elevs)
                    if min_elev < current_elev: #we've expanded our view and found a lower point
                        row,col = coord_list[elevs.index(min_elev)]
                        break
                    elif 0 in [LDD[c] for c in coord_list]:  #we found a proximal outlet anyway
                        row,col = [c for c in coord_list if LDD[c] == 0][0]
                        break
                    view_dist += 1           
    return drain_num
        
def CalculateHAND(nrows, ncols, DEM, network_elevs, nearest_drainage):
    """Calculate the height above (nearest_drainage) for all points in the (DEM).
    (network_elevs) contains the elevations of the numbered points within the flow network."""
    HAND_estimates = np.zeros(shape=(nrows,ncols))
    for row in range(nrows):
        for col in range(ncols):
            HAND_estimates[row,col] = max(DEM[row,col] - network_elevs[nearest_drainage[(row,col)]],0)
    return HAND_estimates   

def ClassifyHAND(HAND_estimates, nrows, ncols, HAND_thresholds):
    """Given a numpy array of (HAND_estimates)...measured as elevation above the nearest drainage
    sources, return, for each row,col in (nrows),(ncols), a classification between 0 - flow-path and 
    other low-lying land, and 3 - peak as a function of the elevation difference percentiles
    defined in (HAND_thresholds)."""
    HAND_classes, max_HAND = np.zeros(shape=(nrows,ncols)), np.max(HAND_estimates)
    thresholds = [thresh * max_HAND for thresh in HAND_thresholds]
    for row in range(nrows):
        for col in range(ncols):
            HAND_classes[row,col] = np.sum([HAND_estimates[row,col] > thresh for thresh in thresholds])
    return HAND_classes, thresholds     
          
def DefineDrainageNetwork(drainage_areas, LDD, nrows, ncols, min_drainage):
    drainage_network = {} #to be filled with points along said network, each assigned a unique id_no.
    id_counter = 1 #the first point gets a tag of one
    for row in range(nrows):
        for col in range(ncols):
            if drainage_areas[row,col] >= min_drainage or LDD[row,col] == 0:
                drainage_network[(row, col)] = id_counter; id_counter += 1
    return drainage_network

def GetDrainNetworkElevs(drainage_network, DEM):
    network_elevs = {}
    for key in drainage_network.keys():
        network_elevs[drainage_network[key]] = DEM[key]
    return network_elevs
            
def RemoveSinksFromLDD(LDD, collapsed_DEM, nrows, ncols, flow_dirs, inv_flow_dirs, max_sink_search_dist):
    for row in range(nrows):
        for col in range(ncols):  
            if LDD[row,col] == 0 and not IsEdge(row, col, nrows, ncols): #this is currently a sink (and not on a boundary)
                print("Breaching sink at", row, col)#; debug_LDD = np.copy(LDD)                
                LDD, outlet_row, outlet_col = ReRouteSink(LDD, collapsed_DEM, row, col, 
                        nrows, ncols, flow_dirs, inv_flow_dirs, max_sink_search_dist)
                print("Sink routed to", outlet_row, outlet_col)
    return LDD

def ReRouteSink(LDD, collapsed_DEM, row, col, nrows, ncols, flow_dirs, inv_flow_dirs, max_sink_search_dist):
    """Given that a [(row),(col)] coordinate is deemed a sink...and not located on the DEM's edge, we must 
    re-route its flow.  We will search up to (max_sink_search_dist) squares in all directions, including the diagonal
    to locate a point of lower elevation...if we find one, direct flow there.  Else, allow a sink."""
    view_dist, point_elev = 2, collapsed_DEM[row,col]
    while view_dist <= max_sink_search_dist:
        ring_points = GetRingAtDistance(view_dist, row, col, nrows, ncols)
        ring_elevs = [collapsed_DEM[coords] for coords in ring_points]
        min_elev = min(ring_elevs)
        if min_elev < point_elev:
            min_coords = ring_points[ring_elevs.index(min_elev)]
            LDD = AlterLDDFromSink(LDD, row, col, min_coords[0], min_coords[1], 
                                   inv_flow_dirs, flow_dirs, collapsed_DEM)
            return LDD, min_coords[0], min_coords[1]
        view_dist += 1
    print("NO LOWER POINTS WITHIN", max_sink_search_dist, "SQUARES"); return LDD, row, col        
    
def AlterLDDFromSink(LDD, row, col, outlet_row, outlet_col, inv_flow_dirs, 
                     flow_dirs, collapsed_DEM):
    """If we are breaching a sink at (row),(col) by forcing flow to (outlet_row), (outlet_col),
    alter the (LDD) to direct the flow accordingly."""
    while (row,col) != (outlet_row, outlet_col):
        row_diff, col_diff = outlet_row-row, outlet_col-col
        possible_flows = [inv_flow_dirs[k] for k in inv_flow_dirs.keys() if 
                    (k[1] in [0,np.sign(row_diff)] and k[0] in [0,np.sign(col_diff)])]
        possible_new_elevs = [collapsed_DEM[GetF_i_j(flow_dirs, flow_dir, row, col)] for flow_dir in possible_flows]
        new_flow_dir = possible_flows[possible_new_elevs.index(min(possible_new_elevs))]
        LDD[row,col] = new_flow_dir
        row,col = GetF_i_j(flow_dirs, new_flow_dir, row, col)
    return LDD
        
def GetRingAtDistance(view_dist, row, col, nrows, ncols):
    north_edge, south_edge = max(0,row-view_dist), min(nrows-1, row+view_dist)
    west_edge, east_edge = max(0,col-view_dist), min(ncols-1,col+view_dist)
    north_coords = [(north_edge, i) for i in range(west_edge, east_edge + 1)]
    west_coords = [(i, west_edge) for i in range(north_edge+1, south_edge)]    
    south_coords = [(south_edge, i) for i in range(west_edge, east_edge + 1)]
    east_coords = [(i, east_edge) for i in range(north_edge+1, south_edge)]
    return south_coords + east_coords + north_coords + west_coords    

def IsEdge(row, col, nrows, ncols):
    return (True if (row == 0 or row == nrows-1 or col == 0 or col == ncols-1) else False)

def IsInBox(row, col, nrows, ncols):
    return (True if (row >= 0 and row <= nrows-1 and col >= 0 and col <= ncols-1) else False)    

def WriteNPArray(np_array, path_to_file, file_name):
    np.savetxt(os.path.join(path_to_file,file_name), np_array, newline = '\n', fmt='%.0f')
    return None
        
def GetJSON(file_path):
    with open(file_path) as inputfile:
        return json.load(inputfile)    
    
def DumpJSON(D, file_path):
    with open(file_path, 'w') as outfile:
        json.dump(D, outfile)
    return None

def GenerateSubDEM(DEM, min_northing, min_easting, nrows, ncols, cellsize, southern_edge,
    northern_edge, western_edge, eastern_edge, aggregation_fac, path_to_DEM, sub_name):
    """Generate a smaller DEM using northing/easting boundaries, and an aggregation factor
    to clump into larger squares (if necessary)."""
    #CHECK THAT OUR CHOSEN BOUNDARIES LIE WITHIN THE DEM...if not, truncate to the edge
    southern_edge = max(min_northing, southern_edge)
    northern_edge = min(min_northing + (nrows-1)*cellsize, northern_edge)
    western_edge = max(min_easting, western_edge)
    eastern_edge = min(min_easting + (ncols-1)*cellsize, eastern_edge)
    #DETERMINE SUB_ROWS OF DEM...and store the actual northing/easting of the corners
    min_row_ind, max_row_ind = nrows - int((northern_edge - min_northing)/cellsize + 0.5), nrows - int((southern_edge - min_northing)/cellsize)
    min_col_ind, max_col_ind = int((western_edge - min_easting)/cellsize), int((eastern_edge - min_easting)/cellsize + 0.5)  
    sub_DEM = DEM[min_row_ind:max_row_ind, min_col_ind:max_col_ind]; del(DEM) #save memory
    min_northing, max_northing = (nrows-max_row_ind)*cellsize + min_northing, (nrows-min_row_ind)*cellsize + min_northing
    min_easting, max_easting = min_easting + min_col_ind*cellsize, min_easting + max_col_ind*cellsize
    #Now collapse the DEM, returning the new coordinates as well        
    collapsed_DEM, min_northing, min_easting, max_northing, max_easting = CollapseDEM(sub_DEM, aggregation_fac, 
        len(sub_DEM), len(sub_DEM.T), min_northing, max_northing, min_easting, max_easting, cellsize)
    WriteNPArray(collapsed_DEM, path_to_DEM, sub_name + '.txt')
    return collapsed_DEM, min_northing, min_easting, max_northing, max_easting

def CollapseDEM(DEM, aggregation_fac, nrows, ncols, min_northing, max_northing, min_easting, max_easting, cellsize):
    """Given a numpy array for the (DEM), average into larger subsquares by an
    (aggregation_fac).  Ex: if aggregation_fac = 5, average the domain into 5x5 squares and
    decrease its size by a factor of 25 in so doing.  Extra squares in the north-east corner
    may be discarded with this procedure"""
    row_cuts, col_cuts = int(float(nrows)/aggregation_fac), int(float(ncols)/aggregation_fac)
    collapsed_DEM = np.zeros(shape=(row_cuts, col_cuts))
    for row in range(row_cuts):
        for col in range(col_cuts):
            sub_piece = DEM[row*aggregation_fac:(row+1)*aggregation_fac, col*aggregation_fac:(col+1)*aggregation_fac]
            collapsed_DEM[row,col] = int(np.mean(sub_piece))
    min_northing = min_northing + cellsize*aggregation_fac/2 
    min_easting = min_easting + cellsize*aggregation_fac/2
    max_northing = min_northing + cellsize*aggregation_fac*(row_cuts-1)
    max_easting = min_easting + cellsize*aggregation_fac*(col_cuts-1)
    return collapsed_DEM, min_northing, min_easting, max_northing, max_easting

def GetLegendCaptions(thresholds, unit, max_HAND):
    """Given the (thresholds) that delineate classifications with the HAND model, return text
    for those intervals.  A (unit) abbreviation is tacked onto the string, e.g. 'm' for meters."""
    legend, L = [], len(thresholds)
    for i,thresh in enumerate(thresholds + [-1]):
        if i == 0: #the lowest bracket - bottom is defined to be 0
            legend.append('0-' + str(int(thresh)) + unit)
        elif i == L: #the highest bracket, the top is defined as the max HAND value
            legend.append(str(int(thresholds[i-1])) + '-' + str(int(max_HAND)) + unit)
        else:
            legend.append(str(int(thresholds[i-1])) +'-'+ str(int(thresh)) + unit)
    return legend

def PlotHAND(HAND_classes, min_northing, min_easting, square_size, nrows, ncols, cmap, 
             legend_captions, path_to_images, sub_name):
	"""Visualize the results from the HAND model...first on a flat image, then with relief"""
	northing_to_easting_ratio, L, max_class = float(nrows)/ncols, len(legend_captions), np.max(HAND_classes)
	legend_ticks = [max_class/L/2 + max_class/L*i for i in range(L)] 
	plt.figure(figsize = (int(20/northing_to_easting_ratio), 20))
	#now, correct the grid to no longer reflect centers of squares, but corners of squares...
	grid_northings = [min_northing -square_size/2 + square_size*i for i in range(nrows+1)][::-1]
	grid_eastings = [min_easting -square_size/2 + square_size*i for i in range(ncols+1)]
	site_extent = [min_easting -square_size/2, min_easting + (ncols-0.5)*square_size, min_northing-square_size/2, min_northing + (nrows-0.5)*square_size]
	#generate pcolor plot 
	X,Y = np.meshgrid(grid_eastings, grid_northings) 
	plt.pcolor(X, Y, HAND_classes,cmap=cmap)
	cbar = plt.colorbar(); cbar.set_ticks(legend_ticks); cbar.set_ticklabels(legend_captions); cbar.ax.tick_params(labelsize = 24)
	plt.axis(site_extent); plt.tick_params(labelsize = 20)
	plt.savefig(os.path.join(path_to_images, sub_name + '.png'))
	plt.close()	#to avoid using additional RAM					
	return None	
 
def GetHANDClassAtPoint(HAND_classes, northing, easting, min_northing, min_easting, cellsize, nrows, ncols):
    """Given an np.array (HAND_classes), specifications of the (min_northing), (min_easting), the number of rows and columns,
    (nrows) and (ncols), and (cellsize) of that array, return the class associated with the (northing), (easting) given."""
    max_northing, max_easting = min_northing + (nrows-1)*cellsize, min_easting + (ncols-1)*cellsize 
    if northing < min_northing:
        print("Northing of", northing, "is below the southern boundary of", min_northing); northing = min_northing
    if easting < min_easting:
        print("Easting of", easting, "is beyond the western boundary of", min_easting); easting = min_easting
    if northing > max_northing:
        print("Northing of", northing, "is above the northern boundary of", max_northing); northing = max_northing
    if easting > max_easting:
        print("Easting of", easting, "is beyond the eastern boundary of", max_easting); easting = max_easting
    HAND_row, HAND_col = int(round((max_northing - northing)/cellsize, 0)), int(round((easting - min_easting)/cellsize, 0))
    return int(HAND_classes[HAND_row, HAND_col])    
    
def GetHandFileInfo(HAND_name):
    """HAND classification files are traditionally named as 'HAND_classes_minnorthing_mineasting_cellsize_nrows_ncols',
    this function should return min_northing, min_easting, cellsize, nrows, and ncols"""
    split_name = HAND_name.split('_')
    return [int(i) for i in split_name[2:]]    

if __name__ == "__main__":
    null, southern_edge, northern_edge, western_edge, eastern_edge, aggregation_fac, sub_name = sys.argv
    #Either create or read the sub_DEM used for the HAND analysis        
    sub_DEM_path = os.path.join(_path_to_DEM, sub_name + '.txt')
    if os.path.exists(sub_DEM_path):
        collapsed_DEM = np.loadtxt(sub_DEM_path); cellsize = _DEM_cellsize
    else:
        feature_json = GetJSON(os.path.join(_path_to_DEM, _file_name + '.json'))        
        min_northing, min_easting = feature_json['yllcorner'], feature_json['xllcorner']
        nrows, ncols, cellsize = feature_json['nrows'], feature_json['ncols'], feature_json['cellsize']
        DEM = ConvertDEM_json_to_np_array(_json_params, feature_json, min_northing, nrows, ncols, cellsize)
        del(feature_json) #save memory
        collapsed_DEM, min_northing, min_easting, max_northing, max_easting = GenerateSubDEM(DEM, min_northing, min_easting, 
        nrows, ncols, cellsize, southern_edge, northern_edge, western_edge, eastern_edge, aggregation_fac, _path_to_DEM, sub_name)
     
     
    nrows, ncols, cellsize = len(collapsed_DEM), len(collapsed_DEM.T), cellsize*int(aggregation_fac) 
    LDD = DefineLDD(collapsed_DEM, nrows, ncols, _flow_dirs, _inv_flow_dirs) #in which direction does each square flow?
    LDD = RemoveSinksFromLDD(LDD, collapsed_DEM, nrows, ncols, _flow_dirs, _inv_flow_dirs, _max_sink_search_dist) #let's make sure our flows EXIT the domain
    drainage_areas = Generate_all_A_i_j(nrows, ncols, LDD, _flow_dirs, _max_A_iterations) #how much area does each location drain?
    drainage_network = DefineDrainageNetwork(drainage_areas, LDD, nrows, ncols, _min_drainage) #which points are on the flow path?
    network_elevs = GetDrainNetworkElevs(drainage_network, collapsed_DEM)
    nearest_drainage = Generate_all_I_i_j(nrows, ncols, LDD, _flow_dirs, drainage_network, collapsed_DEM, _max_sink_search_dist)
    HAND_estimates = CalculateHAND(nrows, ncols, collapsed_DEM, network_elevs, nearest_drainage)
    HAND_classes, thresholds = ClassifyHAND(HAND_estimates, nrows, ncols, _HAND_thresholds)
    PlotHAND(HAND_classes, _southern_edge, _western_edge, cellsize, nrows, ncols, _cmap, 
             GetLegendCaptions(thresholds, _unit, np.max(HAND_estimates)), 
             os.path.join(_path_to_DEM,"..","HAND"), sub_name)
    
####OLD DRAINAGE NETWORK GENERATION CODE
"""            max_area_drained = 9999 #place-holder, will ultimately decrease
            if LDD[row,col] == 0: #this is an outlet, and therefore, by definition, a member of the network
                network_row, network_col = row, col                
                while max_area_drained >= min_drainage:
                    drainage_network[(network_row, network_col)] = id_counter
                    N_i_j = GetN_i_j(network_row, network_col, nrows, ncols) #which points neighbor this one on the flow network
                    id_counter += 1
                    drain_areas = [(drainage_areas[p[0],p[1]] if p not in drainage_network.keys() else -1) for p in N_i_j]
                    max_area_drained = np.max(drain_areas)
                    upstream_point = N_i_j[drain_areas.index(max_area_drained)]
                    network_row, network_col = upstream_point        """