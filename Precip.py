# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 12:02:11 2016
@author: Evan Coopersmith

Modules to deliver precipitation time series data at the locations where requested
"""

import os
import numpy as np
import pandas as pd
import LocalSensors as Local

def ReturnPrecipitationSeries_NLDAS(lat, lon, num_timestamps, current_timestamp, interval = 'hour'):
    """Given a (lat), (lon) coordinate, a number of timestamps to retrieve (num_timestamps),
    the (current_timestamp) from which we must begin the look back, and the (interval) with which 
    we choose to retrieve - the default is 'hour', and otherwise we can retrieve daily totals.
    For example, if interval = 'hour', the current time is 12:00AM on Jan3rd, and num_timestamps
    is, say 12, this function should return a list of the form [P_12pm, P_1pm, P_2pm, ..., P_11pm]
    all from Jan2nd...**In MM.***"""
    return precip_list, years, DOYs
    
def ReturnPrecipitationSeries_InSitu(northing, easting, precip_dict, p_col, current_timestamp, num_timestamps):
    """Using the existing (local_sensors) dictionary, the column name for precipitation in each site's data frame (p_col), 
    the (current_timestamp) in datetime format, and the number of timestamps to look back, return the precipitation time series, 
    at a (northing), (easting) along with years, and days-of-years in a dataframe."""
    p_northings, p_eastings = [],[]
    for key in precip_dict.keys():
        p_northings.append(precip_dict[key]['Northing']); p_eastings.append(precip_dict[key]['Easting'])
    p_weights = GetPrecipWeights(p_northings, p_eastings, northing, easting)
    last_P_ind = GetLastPrecipInd(precip_dict[key]['precip_df'], current_timestamp, 'Year', 'DOY')
    first_P_ind = max(0, last_P_ind - num_timestamps)
    precip_list, years, DOYs = [], [], []
    for p_ind in range(first_P_ind, last_P_ind):
        precip_list.append(np.sum([w*max(precip_dict[k]['precip_df'][p_col][p_ind],0) for w,k in zip(p_weights, precip_dict.keys())]))
        years.append(precip_dict[key]['precip_df']['Year'][p_ind]); DOYs.append(precip_dict[key]['precip_df']['DOY'][p_ind]);
    return pd.DataFrame({'P':precip_list, 'Year':years, 'DOY':DOYs})
    
def GetPrecipWeights(p_northings, p_eastings, northing, easting):
    p_distances = [max((n - northing)**2 + (e - easting)**2,1000) for n,e in zip (p_northings, p_eastings)]
    inv_dists = [1/p for p in p_distances]; inv_dist_sum = np.sum(inv_dists) 
    return [d/inv_dist_sum for d in inv_dists]
    
def GeneratePrecipSeriesDictionary(local_precip_sensors, p_col, year_col, DOY_col, num_timestamps, last_time, path_to_precip):
    """This function exists largely to function in a world without ubiquitous precipitation data.  In this case,
    a data frame of (local_precip_sensors) is provided - sites that contain precip time series data available in
    our database, (p_col) denotes the column holding precipitation data, (year_col) and (DOY_col) denote year, and day-of-year
    columns, respectively. The most recent datetime is denoted by (last_time) and the number of timestamps to look back (num_timestamps) is measured in hours."""    
    precip_dict = {}
    for network, site_code, northing, easting in zip(local_precip_sensors.Network, local_precip_sensors.Site_Code, 
                                                     local_precip_sensors.Northing, local_precip_sensors.Easting):
        name_tag = ((network + "_" + str(site_code)) if network in ['CRN', 'SCAN'] else site_code)
        sensor_df = pd.read_csv(os.path.join(path_to_precip, network, name_tag + '.csv'))
        last_P_ind = GetLastPrecipInd(sensor_df, last_time, year_col, DOY_col)
        if last_P_ind > 0:
            key = network + '_' + site_code; precip_dict[key] = {}
            p_fac = 25.4 if network == 'SCAN' else 1 #SCAN precip is in inches
            first_P_ind = max(last_P_ind - num_timestamps, 0)
            precip_dict[key]['Northing'] = northing; precip_dict[key]['Easting'] = easting
            Ps, DOYs = [p*p_fac for p in sensor_df[p_col][first_P_ind: last_P_ind + 1]], [DOY for DOY in sensor_df[DOY_col][first_P_ind: last_P_ind + 1]]            
            Years = [p*p_fac for p in sensor_df[year_col][first_P_ind: last_P_ind + 1]]   
            precip_df = pd.DataFrame({"P":Ps, "DOY":DOYs, "Year":Years}); precip_dict[key]['precip_df'] = precip_df   
    return precip_dict
    
def GetLastPrecipInd(sensor_df, last_time, year_col, DOY_col):
    year, DOY = Local.ConvertDatetime_to_year_and_DOY(last_time)
    sub_df = sensor_df[sensor_df[year_col] == year]
    if len(sub_df) == 0: return -1
    sub_df = sub_df[np.logical_and(sub_df[DOY_col] < DOY + 0.001, sub_df[DOY_col] > DOY - 0.001)]
    return (int(sub_df.index[0]) if len(sub_df > 0) else -1)    
    
def AdjustPrecipUnit(df, p_col, adjust_fac):
    """If a data frame (df) contains a column for a precipitation time series (p_col) whose unit needs to be adjusted,
    return that column altered by a multiplicative constant (adjust_fac).  A deep copy is required to avoid altering
    the original data frame."""
    df2 = Local.very_deep_copy(df); df2[p_col] = [p * adjust_fac for p in df2[p_col]]
    return df2
    