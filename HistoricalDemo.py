# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 19:06:29 2016
@author: Evan Coopersmith

Historical Demo module capable of generating estimates over a space incorporating:
Local Sensors, Non-Local Sensors, Calibrated Models, Topo, Texture, HydroClasses
"""
import Env
import os
import sys
import numpy as np
import SoilMoisture_CalVal as SM
import LocalSensors as Local
import Precip
import HAND
import SpatialProjectionFunctions as SPF
import pandas as pd
import datetime as dt

_hydro_class = 'IACJ' #for a demo, we'll assume hydroclimatic class does not vary within the test space
_zone = '12R' #assume a northing/easting zone that is unchanged within the context of the demo
_depth_cm = 5
_precip_sensors = ['WG0' + str(i) for i in range(1,10)] + ['WG' + str(i) for i in range(10,20)]
_min_hours = 1000 #how many hours before the start of the analysis must precipitation be available
_needed_columns = {'SMAPVEX04': ['VSM','TP_VSM_gc_avg', 'TP_VSM_ssc_avg', 'SMest_5'],
                   'SMAPVEX15': ['Gravimetric','AZ_Logger','SMest_5']}
_metrics = ['Corr', 'RMSE', 'RMSE_off', 'offset']

def GetBackgroundInformation(historical_request_file_name, output_file_name, HAND_name, texture_json_name):
    """Return a dataframe containing sensors, locations, calibrated parameters, heuristics, topo_classes, texture_classes,
    and hydro_classes.  Next, return a np.array of the HAND classifications to be used in this study.  Next, return a
    .json file containing textural classifications.  Finally, return a dataframe containing northings, eastings, and 
    datetime information for which estimates will be requested."""
    historical_requests = SM.GetFlatFile(os.environ['path_to_historical_requests'], historical_request_file_name + '.csv')
    all_calibrated_sites = SM.GetFlatFile(os.environ['path_to_local_sensors'], os.environ['param_file_name'])
    HAND_classes = np.loadtxt(os.path.join(os.environ['path_to_topos_and_textures'], HAND_name + '.txt'))
    texture_classes = HAND.GetJSON(os.path.join(os.environ['path_to_topos_and_textures'], texture_json_name + '.json'))
    return historical_requests, all_calibrated_sites, HAND_classes, texture_classes
                
def GetHistoricalTemporalScale(historical_requests, year_col, DOY_col, min_hours):
    """Given a data frame containing requests for estimates at various locations and times (historical_requests), along with the
    column in which the year and day-of-year are found (year_col), (DOY_col), return, as datetimes, the first and last required
    times, along with the number of hours between them + a (min_hours) column...meaning precip is needed BEFORE the first time_stamp"""
    max_year, min_year = np.max(historical_requests[year_col]), np.min(historical_requests[year_col])
    last_year_df, first_year_df = historical_requests[historical_requests[year_col] == max_year], historical_requests[historical_requests[year_col] == min_year]
    first_DOY, last_DOY = np.min(first_year_df[DOY_col]), np.max(last_year_df[DOY_col])  
    first_time, last_time = Local.ConvertYearDOY_to_Datetime(min_year, first_DOY), Local.ConvertYearDOY_to_Datetime(max_year, last_DOY)               
    num_hours = int((last_time - first_time).total_seconds()/3600) + min_hours
    return last_time, num_hours
    
def ProduceSiteData(local_precip_df, current_time, min_hours):
    """Given a data from of local precip (local_precip_df) and the number of hours of precip (min_hours) 
    required back from (current_time)"""
    last_P_ind = Precip.GetLastPrecipInd(local_precip_df, current_time, 'Year', 'DOY')
    site_data = local_precip_df[(last_P_ind-min_hours):last_P_ind]
    L = len(site_data)
    site_data['SM'] = [0.1 for i in range(L)];  site_data['Invalid'] = [False for i in range(L)]
    return site_data
    
def GetLocalDict_w_LocalSensors(local_dict, similar_local_sensors):
    new_dict = {}
    for network, site_code in zip(similar_local_sensors.Network, similar_local_sensors.Site_Code):
        local_key = network + "_" + site_code
        new_dict[local_key] = local_dict[local_key]
    return new_dict
            
def GenerateEstimatesAtAllPoints(historical_requests, all_calibrated_sites, HAND_classes, texture_classes, output_file_name, depth, min_hours,
    last_time, num_hours, HAND_min_northing, HAND_min_easting, HAND_cellsize, HAND_nrows, HAND_ncols, hydro_class, zone, precip_dict):
    old_northing, old_easting, local_dict, precip_site_series, sm_estimates, topos, textures = -9999,-9999, {}, {}, [], [], []
    for i, (northing, easting, year, DOY) in enumerate(zip(historical_requests.Northing, historical_requests.Easting, historical_requests.Year, historical_requests.DOY)):
        print("Site#", i, "Northing:", northing, "Easting:", easting, "Year:", year, "DOY:", DOY)        
        if (northing != old_northing or easting != old_easting): #we need to re-determine texture/topo info
            lat, lon = SPF.UTMtoLL(23, northing, easting, zone) 
            topo_class = HAND.GetHANDClassAtPoint(HAND_classes, northing, easting, HAND_min_northing, HAND_min_easting, 
                                                  HAND_cellsize, HAND_nrows, HAND_ncols)
            texture_class = Local.GetTextureClass(texture_classes, northing, easting)
            similar_cal_sensors = SM.GetSimilarSensors(all_calibrated_sites, lat, lon, hydro_class, topo_class, texture_class,
                      similar_class_map, similar_texture_map, similar_topo_map) #we can use their parameters
            similar_local_sensors = Local.ChooseLocalSensors(similar_cal_sensors, northing, easting, zone) #we can actually consider their insitu values
            local_dict = Local.GenerateLocalSensorDictionary(similar_local_sensors, os.environ['path_to_local_sensors'], local_dict)
            print("Texture:", texture_class, "Topo:", topo_class, "Models:", [s for s in similar_cal_sensors.Site_Code], "Locals:", [s for s in similar_local_sensors.Site_Code])
        current_time = Local.ConvertYearDOY_to_Datetime(year, DOY); location_key = str(northing) + "_" + str(easting)
        if location_key not in precip_site_series.keys():
            precip_site_series[location_key] = Precip.ReturnPrecipitationSeries_InSitu(northing, easting, precip_dict, 'P', last_time, num_hours)
        site_data = ProduceSiteData(precip_site_series[location_key], current_time, min_hours)    
        sm_estimates.append(SM.ProduceSMEstimateUsingInSituAndModel(similar_cal_sensors, GetLocalDict_w_LocalSensors(local_dict, similar_local_sensors), 
                                                                    northing, easting, depth, site_data, current_time))
        topos.append(topo_class); textures.append(texture_class)        
        old_northing, old_easting = northing, easting #update the values
    historical_requests['SMest_' + str(depth)] = sm_estimates; historical_requests['Topo'] = topos; historical_requests['Texture'] = textures 
    historical_requests.to_csv(os.path.join(os.environ['path_to_historical_requests'], 'HistoricalResults', output_file_name + '.csv'), index = False)
    return historical_requests    

def AnalyzeHistoricalGroupings(file_names):
    path_to_historical_requests = os.path.join('test_data', 'HistoricalAnalyses', 'HistoricalResults')
    for file_name in file_names:
        historical_requests = pd.read_csv(os.path.join(path_to_historical_requests, file_name + '.csv'))
        last_time, num_hours = GetHistoricalTemporalScale(historical_requests, 'Year', 'DOY', _min_hours)
        sm_dict = Precip.GeneratePrecipSeriesDictionary(local_precip_sensors, 'SM', 'Year', 'DOY', num_hours + 48, last_time + dt.timedelta(hours = 24), os.environ['path_to_local_sensors'])
        needed_columns = _needed_columns['SMAPVEX04'] if 'SMAPVEX04' in file_name else _needed_columns['SMAPVEX15']
        site_id_col = [c for c in historical_requests.columns if 'Site' in c][0]
        site_features = ['Topo','Texture','ClosestSensor','Dist',site_id_col,'Northing','Easting','n']
        unique_site_ids = list(ReturnUniqueList(historical_requests[site_id_col]))
        comparison_tuples = GetComparisonTuples(needed_columns, ['Closest', 'InvDistEst'])
        results_dict = GenerateResultsDict(comparison_tuples, site_features, _metrics)
        for site_id in unique_site_ids:
            closest_sm, closest_sensor, dists, inv_dist_interp = [], [], [], []; print("Processing:", site_id)
            sub_requests = historical_requests[historical_requests[site_id_col] == site_id]
            first_ind = sub_requests.index[0]
            northing, easting = sub_requests.Northing[first_ind], sub_requests.Easting[first_ind]
            current_time, num_hours = GetHistoricalTemporalScale(sub_requests, 'Year', 'DOY', 1)
            sm_series = Precip.ReturnPrecipitationSeries_InSitu(northing, easting, sm_dict, 'SM', current_time + dt.timedelta(hours = 12), num_hours + 24) #buffer
            for year, DOY in zip(sub_requests.Year, sub_requests.DOY):
                current_time = Local.ConvertYearDOY_to_Datetime(year, DOY)
                active_sensors = SM.GetAllLocalSensors(sm_dict, northing, easting, current_time)
                sensor, value, dist = SM.GetClosestSensorValue(active_sensors, 'sm_val', 'dist', 'sm_sensor')
                closest_sm.append(value); closest_sensor.append(sensor); dists.append(dist)
                inv_dist_interp.append(sm_series.SM[Precip.GetLastPrecipInd(sm_series, current_time, 'Year', 'DOY')])
            sub_requests['Closest'] = closest_sm; sub_requests['ClosestSensor'] = closest_sensor; sub_requests['Dist'] = dists
            sub_requests['InvDistEst'] = inv_dist_interp; sub_requests = ScaleSM(sub_requests, needed_columns)
            valid_df = GatherValidValues(sub_requests, needed_columns, [0 for c in needed_columns], [1 for c in needed_columns])
            results_dict = EvaluateGroup(valid_df, site_features, results_dict, comparison_tuples, _metrics)
        results_df = pd.DataFrame(results_dict)
        results_df.to_csv(os.path.join(path_to_historical_requests, file_name + 'results.csv'), index = False)
    return None

def ScaleSM(df, columns):
    """Given a (df), if any of the listed (columns) are 0-100 scaled, rather than 0-1, scale them down"""
    for c in columns:
        if np.max(df[c]) > 1:
            df[c] = [v*1.0/100 for v in df[c]]
    return df
    
def GenerateResultsDict(comparison_tuples, other_columns, metrics):
    D = {}    
    for c in comparison_tuples:
        for m in metrics:
            D[c[0] + ':' + c[1] + ':' + m] = []
    for col in other_columns:
        D[col] = []  
    return D
    
def GetComparisonTuples(needed_columns, other_columns):
    """Given a list of (needed_columns), typically the names for the gravimetric, logger, and model estimates and a list of
    (other_columns) with which to compare them, return a list of tuples containing columns to be analyzed as pairs."""
    sm_products_to_compare, all_columns = [], needed_columns + other_columns
    for i, col_one in enumerate(all_columns):
        for j, col_two in enumerate(all_columns):
            if i>j:
                sm_products_to_compare.append((col_one, col_two))
    return sm_products_to_compare
    
def GatherValidValues(df, needed_columns, min_vals, max_vals):
    """Given a dataframe (df), a list of (needed_columns) that must be retained, and numerical criteria
    (min_vals) and (max_vals) into which those values must fall, return a sub_df that contains only the rows
    that meet all of those criteria."""
    for col, min_val, max_val in zip(needed_columns, min_vals, max_vals):
        df = df[np.logical_and(df[col] > min_val, df[col] < max_val)]
    return df

def EvaluateGroup(valid_df, site_features, results_dict, comparison_tuples, metrics):
    L = len(valid_df); first_ind = valid_df.index[0] if L > 0 else -99
    for feature in site_features:
        if feature == 'n':
            results_dict['n'].append(L)
        else:
            results_dict[feature].append(valid_df[feature][first_ind] if L > 0 else -99)
    if L < 3:
        for k in results_dict.keys():
            if k not in site_features:
                results_dict[k].append(-99)
    else:
        for c in comparison_tuples:
            (rho, RMSE_base, RMSE_off, RMSE_all, offset, w, heuristic) = SM.EvaluateModelPerformance(valid_df[c[0]], valid_df[c[1]], 1, 0)
            for m in metrics:
                k = ':'.join([c[0],c[1],m])
                if m == 'Corr':
                    results_dict[k].append(rho)
                elif m == 'RMSE':
                    results_dict[k].append(RMSE_base)
                elif m == 'RMSE_off':
                    results_dict[k].append(RMSE_off)
                elif m == 'offset':
                    results_dict[k].append(offset)                    
    return results_dict
    
    
def ReturnUniqueList(seq, keepstr=True):
    """Given a sequence (seq), return the unique elements therein in list form."""
    t = type(seq)
    if t == str:
        t = (list, ''.join)[bool(keepstr)]
    seen = []
    return t(c for c in seq if not (c in seen or seen.append(c)))

if __name__ == "__main__":
    null, historical_request_file_name, output_file_name, HAND_name, texture_json_name, zone, hydro_class = sys.argv
    Env.AddEnvironmentVariables() #gather environment variables
    HAND_min_northing, HAND_min_easting, HAND_cellsize, nrows, ncols = HAND.GetHandFileInfo(HAND_name)
    similar_class_map, similar_texture_map, similar_topo_map = SM.GetSimilarClassesAndTextures()    
    historical_requests, all_calibrated_sites, HAND_classes, texture_classes = GetBackgroundInformation(historical_request_file_name,
                                                                                    output_file_name, HAND_name, texture_json_name)
    historical_requests = historical_requests.sort(['Northing','Easting']) #avoids re-gathering information
    local_precip_sensors = all_calibrated_sites[all_calibrated_sites.Site_Code.isin(_precip_sensors)]                                                                                
    last_time, num_hours = GetHistoricalTemporalScale(historical_requests, 'Year', 'DOY', _min_hours)
    precip_dict = Precip.GeneratePrecipSeriesDictionary(local_precip_sensors, 'P', 'Year', 'DOY', num_hours, last_time, os.environ['path_to_local_sensors'])
    
    ###Generate the historical estimates
    historical_requests = GenerateEstimatesAtAllPoints(historical_requests, all_calibrated_sites, HAND_classes, texture_classes, output_file_name, _depth_cm, 
    _min_hours, last_time, num_hours, HAND_min_northing, HAND_min_easting, HAND_cellsize, nrows, ncols, hydro_class, zone, precip_dict)

    ###Add information regarding the closest value, the inv-dist^2-estimate, and compare to estimates and gravimetric samples
