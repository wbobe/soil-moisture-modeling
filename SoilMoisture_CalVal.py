# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 09:47:08 2015
@author: Evan Coopersmith

Functions herein allow calibration and validation of soil moisture models
"""
import Env
import os
import sys
import pandas as pd
import numpy as np
import math
import datetime as dt
from scipy.stats.stats import pearsonr
import GeneticAlgorithms_for_SM as GA
import HydroclimaticClassification as HydroClass
import ML_compendium as ML
import LocalSensors as local
import SpatialProjectionFunctions as SPF
import Precip

"""These values are default assumptions required to calibrate system memory"""
_default_triad = [0.015, 0, 0] #for estimation of system memory, a flat evapotranspiration rate
_long_beta_time = 2000 #maximum possible number of hours to consider precip history
_interval = 50 #incremental number of hours to consider when calibrating a sys memory
_thresh = 0.97 #how similar must the shorter beta series be to the longer one to be acceptable?
"""These values are for calibration of the sinusoidal loss function"""
_grow_start, _grow_end = 100, 300 #days of year, between which we'll attempt to estimate SM
_wrong_hydro_class_fac, _wrong_texture_fac, _wrong_topo_fac = 2, 2, 1 #how much are exact matches preferred?
_min_heuristic_score = 0.5 #how poorly can a site perform in validation before it is excluded
_max_similar_sites = 10 #how many sites can we use to gather parameters for a given location?
_worst_match_flex_fac = 1.2 #ensures that weaker matches still play some role
"""Necessary for SM estimation"""
_depth = 5
_ind_vars = ['BetaSeries','SMest_' + str(_depth), 'DOY']
_dep_var = 'SM'

def GetSimilarClassesAndTextures():
    return HydroClass.GetSimilarSensorMap(), HydroClass.GetSimilarTexturesMap(), HydroClass.GetSimilarTopoMap()

def GetFilesInDirectory(dir_path):
	"""Given a directory path (dir_path), return all files contained within that directory, as a list.  
	This function will avoid listing directories or files within directories..."""
	return [f for f in os.listdir(dir_path) if os.path.isfile(os.path.join(dir_path,f))]

def GetFlatFile(path_to_file, file_name):
    """Reads and returns a flat file (file_name) at (path_to_file)."""
    return pd.read_csv(os.path.join(path_to_file, file_name))

def AddEtaSeries(c_vals, times):
	"""This function reads a list of times, then calculates an Eta value from every time stamp given the 
	sinusoidal parameters contained within c_vals"""
	return [EtaSin(time, c_vals) for time in times]

def EtaSin(t, c_vec):
	"""Calculates the value of a sinusoid with a 365 day period given a value from 0 to 364.999
	Leap years are not included, though it should matter very little numerically"""
	return max(c_vec[0] + c_vec[1] * math.sin(2 * math.pi * (t + c_vec[2])/365),0)	
 
def GetBeta(p_Series, e_Series, z, n):
	"""Generates beta values as specified in Pan et al, 2011, eq 9:  The first argument takes a vector of hourly precip
	data of length n or greater.  The second argument takes an eta series, calculated via the sinusoidal evapotranspirative 
	loss function.  z denotes soil depth, n represents the length of the lookback window, in hours.  This should be the 
	minimum value needed to stabilize Beta."""
	n, p_Series = int(n), np.array(p_Series) * 1.0
	p_Series = [(p if not np.isnan(p) else 0) for p in p_Series]	
	e_Series = np.array(e_Series) * 1.0  + 0.00001 #Because vector math is going to be much, much faster
	L, BetaSeries, WetDryRatioSeries = len(e_Series), [0]*n, p_Series/e_Series
	OneMinusExpSeries = [1 - math.exp(-1 * e / z) for e in e_Series]
	WetDry_w_ExpSeries = WetDryRatioSeries * OneMinusExpSeries
	EtaOverZSeries = e_Series / z; ReverseEZSeries = EtaOverZSeries[::-1]
	CumRevEZSeries = ReverseEZSeries.cumsum(); EZSeries = CumRevEZSeries[::-1]
	null_array = np.zeros(shape = (L,n-2)); EtaRevMat = pd.DataFrame(null_array)
	EtaRevMat.index = range(L); EtaRevMat.columns = range(n-2)
	for j in range(min(n-2,L-2)):
		Zs = [0]*(j + 1) + list(EZSeries[0:(L - j - 1)])
		NewCol = EZSeries - np.array(Zs)
		EtaRevMat[j] = NewCol
	EtaRevMat[EtaRevMat > 700] = 700 #To avoid overflow.
	ExpEtaRevMat = EtaRevMat.applymap(lambda x: math.exp(x))
	for i in range(n, L):
		BetaSeries.append(sum(WetDry_w_ExpSeries[(i - n + 2):i][::-1] *
						  np.array(ExpEtaRevMat[(i - 1):i])[0])) #Converting back to an array means grabbing the first element of the "list"
	BetaSeries = np.array(BetaSeries)
	BetaSeries[n:L] = BetaSeries[n:L] + WetDry_w_ExpSeries[n:L]				  
	return BetaSeries ##Note, in the R version, the first n values appear as NA.  In python, they will simply be zeros.	

def GetBetaDiffs(p_Series, e_Series, depth, long, b_Series):
	"""In the case of an extremely wet or dry period, this function will provide the system a 'memory'
	longer than the standard beta series.  This comparison is used for error correction at a later stage."""
	b_long = GetBeta(p_Series, e_Series, depth, long)		
	b_diffs = [x - y for x,y in zip(b_long, b_Series)]
	return [max(b,0) for b in b_diffs] #we cannot have less rain over the past 2000 hours than the past 400, unless data are missing

def GetSMSeries(b_series, parameters):
	"""Applies the diagnostic soil moisture equation (Pan et al, 2012) over all values of the beta series input.
	The parameters are a list of three values: porosity, residual soil moisture, and a constant related to the 
	rate of drainage/drying.  These have already been fit by GA when this function is called."""
	return [DiagSoilEq(b, parameters) for b in b_series]
		
def DiagSoilEq(B, params):
	"""Implements the diagnostic soil moisture equation (Pan et al, 2012).  "params" will contain porosity, 
	residual soil moisture, and a drainage/drying constant, all determined via GA."""
	return params[1] + (params[0] - params[1]) * (1 - math.exp(min(min(params[2] * B, 700) * -1, 700)))

def CalculateSMEstimates(SiteData, parameters, z, sm_col, precip_col, DOY_col):
    """(SiteData) is a pandas data frame.  Its columns must contain soil moisture time series
    information (sm_col, in m3/m3), even if this information is simply a dummy value, precipitation time series
    data (precip_col, in mm), a day-of-year column (DOY_col), and a depth at which to produce estimates (z, in cm).  
    The list, (parameters), contains, in order: [vertical_shift, amplitude, horizontal_shift(days), porosity(%), 
    residual_soil_moisture(%), drainage_rate, system_memory(hrs)]"""
    EtaSeries = AddEtaSeries(parameters[:3], SiteData[DOY_col])
    BetaSeries = GetBeta(SiteData[precip_col], EtaSeries, z, parameters[-1])		
    SMSeries = GetSMSeries(BetaSeries, parameters[3:6])
    SiteData['BetaSeries'] = BetaSeries; SiteData['SMest_' + str(z)] = [s/100 for s in SMSeries]
    SiteData['ErrorSeq'] = [est - in_situ for est, in_situ in zip(SiteData['SMest_' + str(z)], SiteData[sm_col])]        
    return SiteData

def GetParameters(parameter_file, site_index):
    """Given a (parameter_file) and the (site_index) associated with the relevant set of parameters, 
    return a parameter list of the order: [vertical_shift, amplitude, horizontal_shift(days), porosity(%), 
    residual_soil_moisture(%), drainage_rate, system_memory(hrs)]."""
    row = parameter_file.loc[site_index]    
    return [row.Vert_Shift, row.Amplitude, row.Horiz_Shift, row.Porosity, row.Res_SM, row.Drainage, row.n]

def GetSimilarSensors(parameter_file, lat, lon, hydroclimatic_class, topo_class, texture_class,
                      similar_class_map, similar_texture_map, similar_topo_map):
    """Given a (parameter_file), the relevant (lat) and (lon) where estimates are required, the (hydroclimatic_class),
    (topo_class), and (texture_class), return the indices of parameter_file for site whose parameters can be used to
    generate estimates.  Additionally, return a weight associated with each set of parameters.
    (similar_class_map) and (similar_texture_map) define which other classes and textures are sufficiently similar."""        
    acceptable_textures = [texture_class] + similar_texture_map[texture_class]
    acceptable_hydro_classes = [hydroclimatic_class] + similar_class_map[hydroclimatic_class]
    acceptable_topo_classes = [topo_class] + similar_topo_map[topo_class]
    sub_parameter_file = parameter_file.loc[parameter_file.Texture.isin(acceptable_textures)]
    sub_parameter_file = sub_parameter_file.loc[sub_parameter_file.Class.isin(acceptable_hydro_classes)]
    sub_parameter_file = sub_parameter_file.loc[sub_parameter_file.Topo.isin(acceptable_topo_classes)]
    sub_parameter_file = sub_parameter_file[sub_parameter_file.Heuristic > _min_heuristic_score]
    sub_parameter_file['Dist'] = [math.sqrt((s_lat - lat)**2 + (s_lon - lon)**2) for s_lat, s_lon in 
                                  zip(sub_parameter_file.Lat, sub_parameter_file.Lon)]   
    sub_parameter_file = CalculateSiteSimilarityHeuristic(sub_parameter_file, hydroclimatic_class, topo_class, texture_class)                          
    if len(sub_parameter_file) > _max_similar_sites: sub_parameter_file = sub_parameter_file[:_max_similar_sites]
    return CalculateSiteWeights(sub_parameter_file)
    
def CalculateSiteWeights(sub_parameter_file):
    """Given the smaller list of similar sites (sub_parameter_file) to be used for producing soil moisture estimates, calculate
    an additional column determining the weight of each individual site's parameters in producing estimates."""  
    max_sim_val = np.max(sub_parameter_file.Similarity)*_worst_match_flex_fac
    weights = [max_sim_val - s for s in sub_parameter_file.Similarity]
    tot_weight = np.sum(weights)
    sub_parameter_file['Weight'] = [w/tot_weight for w in weights]
    return sub_parameter_file
    
def CalculateSiteSimilarityHeuristic(similar_sites, hydroclimatic_class, topo_class, texture_class):
    """Given a data frame, (similar_sites), the (hydroclimatic_class), (topo_class), and (texture_class)
    associated with the location where predictions are required, return the similiarity_score for 
    each site deemed a possible match."""
    #Topo class is not implemented in this version...
    similarity_score = []    
    for hydro_class, texture, topo, distance, heuristic in zip(similar_sites.Class, similar_sites.Texture, 
        similar_sites.Topo, similar_sites.Dist, similar_sites.Heuristic):
        hydro_class_fac = 1 if hydro_class == hydroclimatic_class else _wrong_hydro_class_fac
        texture_class_fac = 1 if texture == texture_class else _wrong_texture_fac
        #topo_class_fac = 1 if topo == topo_class else _wrong_topo_fac       
        similarity_score.append(hydro_class_fac*texture_class_fac*distance*(1-heuristic)**2)
    similar_sites['Similarity'] = similarity_score     
    return similar_sites.sort('Similarity')

def GetSystemMemory(default_triad, p_Series, z, long_time, interval, thresh, times):
	"""How far back must we look before we are confident that antecedent rain conditions are minimally relevant?"""
	default_etas = AddEtaSeries(default_triad, times)
	L = len(default_etas)
	BetaLong = GetBeta(p_Series, default_etas, z, long_time)
	for i in range(1, int(long_time / interval)):
		sys_mem = i * interval; print("testing", sys_mem, "as a system memory")
		TestBeta = GetBeta(p_Series, default_etas, z, sys_mem)
		r_val = pearsonr(TestBeta[long_time:L],BetaLong[long_time:L])[0] #The first value is 'rho' the second is a significance
		print(r_val)
		if(r_val > thresh): return sys_mem
	return long_time
    
def FittingEta(n, z, invalid_seq, sm_seq, time_stamps, p_series, null_b_series):
	"""Given the appropriate length of the system's memory, calculate the appropriate parameters for the EtaSeries
	using the genetic algorithm scripts developed at UIUC circa 2012"""
	return GA.FullGA(z, n, _grow_start, _grow_end, invalid_seq, "FittingEta", sm_seq, time_stamps, p_series, null_b_series)
	
def FittingDiag(n, z, v, a, h, InvalidSeq, TrueSeq, TimeStamps, PrecipSeries):
	"""Given the appropriate length of the system's memory, calculate the appropriate parameters for the Diagnostic Soil 
	Moisture Equation (Pan et al, 2012) using the genetic algorithm scripts developed at UIUC circa 2012"""
	EtaSeries = AddEtaSeries([v,a,h], TimeStamps)
	BetaSeries = GetBeta(PrecipSeries, EtaSeries, z, n)
	return GA.FullGA(z, n, _grow_start, _grow_end, InvalidSeq, "FittingDiag", TrueSeq, TimeStamps, PrecipSeries, BetaSeries)
    
def CalibrateSMParameters(SiteData, z, sm_col, p_col, DOY_col):
    """Given a (SiteData) data frame that contains soil moisture information in an (sm_col), precipitation information
    in (p_col), and day-of-year information in (DOY_col), calibrate soil moisture parameters for depth (z, cm).""" 
    if max(SiteData[sm_col] < 2): SiteData[sm_col] = [sm*100 for sm in SiteData[sm_col]] #To ensure 0-100 range for calibration
    n = GetSystemMemory(_default_triad, SiteData[p_col], z, _long_beta_time, _interval, _thresh, SiteData[DOY_col])
    v,a,h = FittingEta(n, z, SiteData.Invalid, SiteData[sm_col], SiteData[DOY_col], SiteData[p_col], SiteData.BetaSeries)
    porosity, res_sm, drainage = FittingDiag(n, z, v, a, h, SiteData.Invalid, SiteData[sm_col], SiteData[DOY_col], SiteData[p_col])
    return [v, a, h, porosity, res_sm, drainage, n]    
    
def GenerateInvalidSeq(SiteData, sm_col, p_col, DOY_col):
    """Given a (SiteData) data frame that contains soil moisture information in an (sm_col), precipitation information
    in (p_col), and day-of-year information in (DOY_col), determine which indices are not valid for calibration/validation""" 
    invalid_seq = []
    for sm, p, DOY in zip(SiteData[sm_col], SiteData[p_col], SiteData[DOY_col]):
    #criteria that renders a timestamp unusable
        invalid_seq.append(True if (not (DOY <= _grow_end and DOY >= _grow_start) or sm <= 0 or p < 0) else False) 
    return invalid_seq   
    
def GetParameters_and_Weight_of_CalSensor(ind, similar_sensors):
    """Given the index (ind) within a data frame of (similar_sensors), return the parameters
    of the diagnostic soil moisture equation along with the weight of that sensor's estimate."""    
    v, a, h = similar_sensors.loc[ind]['Vert_Shift'], similar_sensors.loc[ind]['Amplitude'], similar_sensors.loc[ind]['Horiz_Shift']
    por, res, drain = similar_sensors.loc[ind]['Porosity'], similar_sensors.loc[ind]['Res_SM'], similar_sensors.loc[ind]['Drainage']
    n, w = similar_sensors.loc[ind]['n'], similar_sensors.loc[ind]['Weight']
    return v,a,h,por,res,drain,n,w    
    
def ProduceSMFromSimilarCalSensors(similar_sensors, depth, site_data):
    """Given (similar_sensors) whose parameters can be used to produce soil moisture estimates, return an
    estimate using the weighted similarities of all suitably similar sensors."""
    sm_estimates, weights = [],[] 
    for ind in similar_sensors.index:
        v,a,h,por,res,drain,n,w = GetParameters_and_Weight_of_CalSensor(ind, similar_sensors)
        print("Estimate from parameters of:", similar_sensors.loc[ind]['Network'], similar_sensors.loc[ind]['Site_Code'], similar_sensors.loc[ind]['Site_Info'])
        recent_site_data = site_data[(n*-1 - 1):] #just enough to generate a SM estimate
        recent_site_data = CalculateSMEstimates(recent_site_data, [v,a,h,por,res,drain,n], depth, 'SM', 'P', 'DOY')
        sm_estimates.append(list(recent_site_data['SMest_' + str(depth)])[-1]); weights.append(w)
    return sm_estimates, weights    

def ProduceSMEstimateUsingInSituAndModel(similar_cal_sensors, local_dict, northing, easting, depth, site_data, current_time):
    sm_estimates = []    
    if len(local_dict) > 0: #available local sensors
        sm_estimates, weights = ProduceSMFromSimilarLocalSensors(similar_cal_sensors, local_dict, northing, easting, depth, site_data, current_time)    
    if len(sm_estimates) == 0: #no available local sensors
        sm_estimates, weights = ProduceSMFromSimilarCalSensors(similar_cal_sensors, depth, site_data)
    location_estimate = EstimateFromWeights(sm_estimates, weights)       
    print("Location estimate:", location_estimate); return location_estimate
    
def ProduceSMFromSimilarLocalSensors(similar_cal_sensors, local_dict, northing, easting, depth, site_data, current_time):
    local_sensor_inv_dist_sqs, estimates_by_sensor, cal_sensor_info = [],[],{}
    for sensor in local_dict.keys(): #for each local sensor whose value ought to be considered
        insitu = local_dict[sensor]; local_network = sensor.split("_")[0]
        local_sensor_inv_dist_sqs.append(1/(max((northing - insitu['northing'])**2 + (easting - insitu['easting'])**2, 1)))
        SM_ind = Precip.GetLastPrecipInd(insitu['site_data'], current_time, 'Year', 'DOY')        
        theta_i_tstar_s = insitu['site_data']['SM'][SM_ind] if SM_ind > 0 else -1
        #print('LOCAL INSITU: ***', theta_i_tstar_s, " *** Gathered from ", local_network, sensor)
        if theta_i_tstar_s > 0:
            estimates_by_model, model_weights = [],[]   
            for ind, cal_network, site_code in zip(similar_cal_sensors.index, similar_cal_sensors.Network, similar_cal_sensors.Site_Code):
                if site_code not in cal_sensor_info.keys(): cal_sensor_info[site_code] = {}; #print('Using cal sensor: ', site_code)
                if 'params_and_weight' not in cal_sensor_info[site_code].keys(): cal_sensor_info[site_code]['params_and_weight'] = GetParameters_and_Weight_of_CalSensor(ind, similar_cal_sensors) 
                v,a,h,por,res,drain,n,w = cal_sensor_info[site_code]['params_and_weight']; key = 'recent_site_data_' + str(n)    
                if key not in cal_sensor_info.keys(): cal_sensor_info[key] = site_data[(n*-1 -1):]                 
                recent_site_data = cal_sensor_info[key]; calc_depth = (2 if 'SCAN' in cal_network else 5)
                #print('Parameters:', v,a,h,por,res,drain,n, "Depth: ", calc_depth, "Site-P (mm): ", np.sum(recent_site_data.P))                
                adj_fac = (1 if 'ARS' in cal_network else round(1/25.4,3)); key = '_'.join(['adj_site_data', str(n), str(adj_fac)])                
                if key not in cal_sensor_info.keys(): cal_sensor_info[key] = Precip.AdjustPrecipUnit(recent_site_data, 'P', adj_fac)
                adj_site_data = cal_sensor_info[key]                
                #print('Site-P, adjusted:', np.sum(adj_site_data.P))
                key = '_'.join(["SMest", str(calc_depth), 'at-site'])
                if key not in cal_sensor_info[site_code].keys(): 
                    adj_site_data = CalculateSMEstimates(adj_site_data, [v,a,h,por,res,drain,n], calc_depth, 'SM', 'P', 'DOY')                 
                    cal_sensor_info[site_code][key] = list(adj_site_data['SMest_' + str(calc_depth)])[-1]
                theta_i_t_m = cal_sensor_info[site_code][key]
                #print('Model estimate using site precip:', theta_i_t_m)
                key = '_'.join(['Site, replaced with', sensor, str(n)])
                if key not in cal_sensor_info.keys():
                    new_site_data = local.very_deep_copy(recent_site_data) #to be manipulated in functions safely
                    cal_sensor_info[key] = ReplaceSitePrecipWSensorPrecip(new_site_data, insitu['site_data']['P'][:SM_ind], 'P', n)
                new_site_data = cal_sensor_info[key]    
                #print('LocalSensor-P:', np.sum(new_site_data.P))
                adj_fac = GetAdjFac(local_network, cal_network); key = '_'.join(['Site, replaced with', sensor, str(n), 'adj', str(adj_fac)]) 
                if key not in cal_sensor_info.keys(): cal_sensor_info[key] = Precip.AdjustPrecipUnit(new_site_data, 'P', adj_fac)
                adj_site_data = cal_sensor_info[key]                
                #print('LocalSensor-P:', np.sum(adj_site_data.P), "Depth:", calc_depth)
                key = '_'.join(['SMest', str(calc_depth), 'using_local', sensor])
                if key not in cal_sensor_info[site_code].keys(): 
                    adj_site_data = CalculateSMEstimates(adj_site_data, [v,a,h,por,res,drain,n], calc_depth, 'SM', 'P', 'DOY')
                    cal_sensor_info[site_code][key] = list(adj_site_data['SMest_' + str(calc_depth)])[-1]
                theta_i_tstar_m = cal_sensor_info[site_code][key]
                #print('Model estimate using local sensor precip:', theta_i_tstar_m)
                estimates_by_model.append(theta_i_t_m - theta_i_tstar_m); model_weights.append(w)
            model_driven_adjustment = EstimateFromWeights(estimates_by_model, model_weights)
            #print("Model suggests location is:", model_driven_adjustment, "different from insitu")
            estimates_by_sensor.append(theta_i_tstar_s + model_driven_adjustment)
        else:
            estimates_by_sensor.append(-99)
    return estimates_by_sensor, local_sensor_inv_dist_sqs    

def GetAdjFac(local_sensor_type, cal_sensor_type):
    """When using the precipitation of a (local_sensor_type) insitu sensor at the location of a site 
    where an estimate is needed, and parameters generated at a (cal_sensor_type) are employed, return a conversion
    factor betweem the local_sensor's precip and the cal_sensor parameter's expected precip unit."""
    if local_sensor_type == cal_sensor_type: #no need to adjust
        return 1
    elif local_sensor_type == 'ARS' and cal_sensor_type != 'ARS': #local sensor P is mm, but the model was calibrated based on in
        return round(1/25.4,3)
    elif local_sensor_type in ['SCAN', 'CRN'] and cal_sensor_type == 'ARS': #local senor P is inches, model calibrated in mm
        return 25.4
    else: #no need to adjust betweem 'CRN' and 'SCAN'
        return 1           
        
def ReplaceSitePrecipWSensorPrecip(new_site_data, p_series, p_col, n):
    """Given a (recent_site_data) frame from the location where SM estimates are desired, 
    return new_site_data that replaces the precipitation series (p_col) at the location itself with the
    precipitation series from the sensor (p_series) for the most recent (n) + 1 time stamps.
    This is the amount of data required for the diagnostic soil moisture equation.
    Missing precip data will be auto-filled with zeros (THIS FEATURE COULD BE CHANGED AT A LATER DATE.)."""
    L = len(new_site_data)
    if len(p_series) > n: #enough data.
       new_site_data[p_col] = list(new_site_data[p_col][0:(L-n)]) + list(p_series[-(n):])
    else: #not enough data        
       missing_p = [0 for i in range(n - len(p_series))] 
       new_site_data[p_col] = missing_p + p_series      
    return new_site_data
        
def EstimateFromWeights(sm_estimates, weights):
    """Given a list of (sm_estimates) and the corresponding (weights) associated with
    each site, produce a weighted estimate of soil moisture at the given location."""
    valid_sms, weights_of_valid_sm = [], []    
    for sm, w in zip(sm_estimates, weights):
        if sm > -1: 
            valid_sms.append(sm); weights_of_valid_sm.append(w)
    total_weight = np.sum(weights_of_valid_sm)
    return np.sum([sm * w/total_weight for sm,w in zip(valid_sms, weights_of_valid_sm)])      

def GetAllLocalSensors(sm_dict, northing, easting, current_time):
    """Given a dictionary of local soil moisture sensors (sm_dict), each of which has a key for northing, easting,
    and a soil moisture data frame with SM, DOY, and Year...at a (current_time), a datetime value, and (northing),
    (easting) coordinate, return the nearest soil moisture value from all active sensors."""
    current_SM_dict = {"sm_sensor" : [], "dist" : [], "sm_val" : []}; SM_ind = -99
    for sensor in sm_dict.keys():
        if SM_ind < -1: SM_ind = Precip.GetLastPrecipInd(sm_dict[sensor]['SM_df'], current_time, 'Year', 'DOY')  #only need this once      
        current_SM_dict['sm_sensor'].append(sensor); current_SM_dict['sm_val'].append(sm_dict[sensor]['SM_df']['SM'][SM_ind])
        current_SM_dict['dist'].append(math.sqrt((sm_dict[sensor]['Northing'] - northing)**2 + (sm_dict[sensor]['Easting'] - easting)**2))
    return pd.DataFrame(current_SM_dict)
    
def GetClosestSensorValue(active_sensors, value_col, distance_col, sensor_col):
    """Given a dictionary of active local soil moisture sensors, or other features (active_sensors), at a given time, return 
    the nearest soil moisture value as defined by (distance_col), or other value given by (value_col) from an active sensor."""
    active_sensors = active_sensors.sort(distance_col)
    for sensor, value, dist in zip(active_sensors[sensor_col], active_sensors[value_col], active_sensors[distance_col]):
        if value > 0: return sensor, value, dist
    return 'None', -99, 99999

def EvaluateModelPerformance(sm_observed, sm_modelled, max_sm, min_sm):
	"""Given two vectors for observed and modelled soil moisture, (sm_observed) and (sm_modelled), 
	crop them such that we consider only pairs where both values are non-negative. Then, return 
	the correlation between the two vectors, the RMSE between the two vectors, the optimal offset, 
	the optimal gain, and the RMSE using the offset alone, then the gain and offset in tandem.
	Finally, using the range of soil moisture (max_sm) and (min_sm), return the heuristic
	for overall performance."""
	sm_obs, sm_model = [],[] 
	for obs, model in zip(sm_observed, sm_modelled):
		if(obs > 0 and model > 0):
			sm_obs.append(obs); sm_model.append(model)
	rho, RMSE_base = GetCorr_and_RMSE(sm_obs, sm_model)
	offset = np.mean(sm_model) - np.mean(sm_obs) #best offset, with no gain involved
	rho, RMSE_off = GetCorr_and_RMSE(sm_obs + offset, sm_model)
	xi = np.array(sm_obs); A,y = np.array([xi, np.ones(len(sm_obs))]), np.array(sm_model) #prepping the regression
	if len(sm_obs) < 2: return(-1,-1,-1,-1,-1,-1,-1)
	w = np.linalg.lstsq(A.T,y)[0] #the gain, the offset as the two returned elements of this function
	if abs(w[0]) >= 100 or abs(w[1]) >= 100: return(-1,-1,-1,-1,-1,-1,-1) #fixes the degenerate case when SM from AMSR doesn't vary...thus no rho
	line = w[0]*xi+w[1];  rho, RMSE_all = GetCorr_and_RMSE(line,y)
	heuristic = rho - (RMSE_base)/(max_sm - min_sm)     
	return (rho, RMSE_base, RMSE_off, RMSE_all, offset, w, heuristic)          
 
def GetCorr_and_RMSE(v1, v2):
	"""given two vectors (v1) and (v2), return their correlation and their RMSE"""
	return pearsonr(v1,v2)[0], math.sqrt(np.mean([(a-b)**2 for a,b in zip(v1,v2)]))

if __name__ == "__main__":   
    #users input lat/lon, from there we should theoretically know: 
    #hydroclimatic_class, topo_class, and texture_class
    null, northing, easting = sys.argv 
    Env.AddEnvironmentVariables() #gather environment variables
    northing, easting, zone, hydroclimatic_class, topo_class, texture_class = 3510000, 590000, '12R', 'IAQ', 2, 6 #FOR TEST PURPOSES
        
    lat, lon = SPF.UTMtoLL(23, northing, easting, zone)
    ###Example soil moisture estimation procedure for a single site with one est of parameters
    parameter_file = GetFlatFile(os.environ['path_to_local_sensors'], os.environ['param_file_name'])
    parameter_list = GetParameters(parameter_file, 124) #random row, number 124 is in WG
    site_data = GetFlatFile(os.environ['path_to_local_sensors'], os.environ['sample_soil_moisture_file'])[3000:10000]
    site_data = CalculateSMEstimates(site_data, parameter_list, 5, 'SM', 'Precip', 'DOY') #for a single site

    #...using multiple sensors that are 'similar'
    similar_class_map, similar_texture_map, similar_topo_map = GetSimilarClassesAndTextures()    
    similar_cal_sensors = GetSimilarSensors(parameter_file, lat, lon, hydroclimatic_class, topo_class, texture_class,
                      similar_class_map, similar_texture_map, similar_topo_map)    
    sm_estimates, weights = ProduceSMFromSimilarCalSensors(similar_cal_sensors, lat, lon, _depth, site_data)
    sm_estimate = EstimateFromWeights(sm_estimates, weights)
    
    #...incorporating local data
    similar_local_sensors = local.GetSimilarLocalSensors(os.environ['path_to_local_sensors'], 
                    os.environ['local_sensor_file_name'], lat, lon, texture_class, topo_class)
    sm_estimate = ProduceSMEstimateUsingInSituAndModel(similar_cal_sensors, similar_local_sensors, lat, lon, _depth, site_data)
    
    ###Example calibration of soil moisture parameters
    site_data['Invalid'] = GenerateInvalidSeq(site_data, 'SM', 'Precip', 'DOY')
    site_data['BetaSeries'] = [0 for i in range(len(site_data))]
    parameter_list = CalibrateSMParameters(site_data, 5, 'SM', 'Precip', 'DOY')
    
    #Example improvement of SM estimates using ML: 
    #(requires reading in the dummy file and estimating SM)
    examples_to_improve = site_data[-10:] #just the last ten examples as a demo
    historical_examples = site_data[500:-10] #grabbing a hypothetical training set
    improved_predictions = ML.KNN(examples_to_improve, historical_examples, _ind_vars, _dep_var)
    
    #Evaluate a model's performance against observations
    (rho, RMSE_base, RMSE_off, RMSE_all, offset, w, heuristic) = EvaluateModelPerformance(site_data.SM, 
        site_data.SMest_5, np.max(site_data.SM), np.min(site_data.SM))
    