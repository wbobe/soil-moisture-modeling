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
import datetime as dt
import SoilMoisture_CalVal as SM
import LocalSensors as Local
import Precip
import HAND
import SpatialProjectionFunctions as SPF
import pandas as pd

_hydro_class = 'IACJ' #for a demo, we'll assume hydroclimatic class does not vary within the test space
_zone = '12R' #assume a northing/easting zone that is unchanged within the context of the demo
_depth_cm = 5
_precip_sensors = ['WG0' + str(i) for i in range(1,10)] + ['WG' + str(i) for i in range(10,20)]
_min_hours = 1000 #how many hours before the start of the analysis must precipitation be available

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
            
def GenerateEstimatesAtAllPoints(historical_requests, all_calibrated_sites, HAND_classes, texture_classes, output_file_name, depth, min_hours,
    last_time, num_hours, HAND_min_northing, HAND_min_easting, HAND_cellsize, HAND_nrows, HAND_ncols, hydro_class, zone, precip_dict):
    old_northing, old_easting, local_dict, precip_site_series, sm_estimates = -9999,-9999, {}, {}, []
    for northing, easting, year, DOY in zip(historical_requests.Northing, historical_requests.Easting, historical_requests.Year, historical_requests.DOY):
        print("Northing:", northing, "Easting:", easting, "Year:", year, "DOY:", DOY)        
        if (northing != old_northing or easting != old_easting): #we need to re-determine texture/topo info
            lat, lon = SPF.UTMtoLL(23, northing, easting, zone) 
            topo_class = HAND.GetHANDClassAtPoint(HAND_classes, northing, easting, HAND_min_northing, HAND_min_easting, 
                                                  HAND_cellsize, HAND_nrows, HAND_ncols)
            texture_class = Local.GetTextureClass(texture_classes, northing, easting)
            similar_cal_sensors = SM.GetSimilarSensors(all_calibrated_sites, lat, lon, hydro_class, topo_class, texture_class,
                      similar_class_map, similar_texture_map, similar_topo_map) #we can use their parameters
            similar_local_sensors = Local.ChooseLocalSensors(similar_cal_sensors, northing, easting, zone) #we can actually consider their insitu values
        local_dict = Local.GenerateLocalSensorDictionary(similar_local_sensors, os.environ['path_to_local_sensors'], local_dict)
        current_time = Local.ConvertYearDOY_to_Datetime(year, DOY); location_key = str(northing) + "_" + str(easting)
        if location_key not in precip_site_series.keys():
            precip_site_series[location_key] = Precip.ReturnPrecipitationSeries_InSitu(northing, easting, precip_dict, 'P', last_time, num_hours)
        site_data = ProduceSiteData(precip_site_series[location_key], current_time, min_hours)    
        sm_estimates.append(SM.ProduceSMEstimateUsingInSituAndModel(similar_cal_sensors, local_dict, northing, easting, depth, site_data, current_time))
    return sm_estimates    

if __name__ == "__main__":
    null, historical_request_file_name, output_file_name, HAND_name, texture_json_name, zone, hydro_class = sys.argv
    Env.AddEnvironmentVariables() #gather environment variables
    HAND_min_northing, HAND_min_easting, HAND_cellsize, nrows, ncols = HAND.GetHandFileInfo(HAND_name)
    similar_class_map, similar_texture_map, similar_topo_map = SM.GetSimilarClassesAndTextures()    
    historical_requests, all_calibrated_sites, HAND_classes, texture_classes = GetBackgroundInformation(historical_request_file_name,
                                                                                    output_file_name, HAND_name, texture_json_name)
    local_precip_sensors = all_calibrated_sites[all_calibrated_sites.Site_Code.isin(_precip_sensors)]                                                                                
    last_time, num_hours = GetHistoricalTemporalScale(historical_requests, 'Year', 'DOY', _min_hours)
    precip_dict = Precip.GeneratePrecipSeriesDictionary(local_precip_sensors, 'P', 'Year', 'DOY', num_hours, last_time, os.environ['path_to_local_sensors'])
    
    sm_estimates = GenerateEstimatesAtAllPoints(historical_requests, all_calibrated_sites, HAND_classes, texture_classes, output_file_name, _depth_cm, 
    _min_hours, last_time, num_hours, HAND_min_northing, HAND_min_easting, HAND_cellsize, nrows, ncols, hydro_class, zone, precip_dict)
