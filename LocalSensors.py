# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 09:04:58 2016
@author: Evan Coopersmith

Accesses, processes, and incorporates information from local in situ resources
"""
import Env
import os
import json
import datetime as dt
import dateutil
import math
import pandas as pd
import HydroclimaticClassification as HydroClass

_insitu_proximity_thresh = 10000 #euclidian northing/easting distance (in m) required for insitu estimate to be considered
_temporal_thresh = 3600*6 #in seconds...how recent must data be to be deployed?

def GetJSON(f_path, f_name):
	"""Read in a JSON of (f_name) at (f_path)."""
	json_data = open(os.path.join(f_path, f_name)).read()
	return json.loads(json_data)

def GetTextureClass(textures, northing, easting):
    """Given a json (textures), return the class associated with the (northing), (easting) given."""
    cellsize, nrows, ncols, min_northing, min_easting = textures['cellsize'], textures['nrows'], textures['ncols'], textures['yllcorner'], textures['xllcorner']
    max_northing, max_easting = min_northing + (nrows-1)*cellsize, min_easting + (ncols-1)*cellsize 
    if northing < min_northing:
        print("Northing of", northing, "is below the southern boundary of", min_northing); northing = min_northing
    if easting < min_easting:
        print("Easting of", easting, "is beyond the western boundary of", min_easting); easting = min_easting
    if northing > max_northing:
        print("Northing of", northing, "is above the northern boundary of", max_northing); northing = max_northing
    if easting > max_easting:
        print("Easting of", easting, "is beyond the eastern boundary of", max_easting); easting = max_easting
    texture_key, texture_ind = max_northing - int(round((max_northing - northing)/cellsize, 0))*cellsize, int(round((easting - min_easting)/cellsize, 0))
    return textures[str(texture_key)][texture_ind]

def very_deep_copy(self):
	"""Make a true 'deep copy' of a data frame object to enable debugging"""
	return pd.DataFrame(self.values.copy(), self.index.copy(), self.columns.copy())
 
def UpdateLocalSensors(path_to_sensor_file, sensor_file_name):
    """Updates local sensor data to the present, incorporating most recently available information,
    then dumps the updated information to the appropriate .json file."""
    #ANDREW, THIS IS ALMOST CERTAINLY GOING TO BE ON YOUR PLATE AT SOME POINT.
    return None

def CalculateHoursBetween(t1, t2):
    """Given (t1) and (t2) in ISO-8601 format, return the number of hours by which t2 trails t1."""        
    if type(t1) == str: t1 = dateutil.parser.parse(t1) #converts ISO-8601 to datetime
    if type(t2) == str: t2 = dateutil.parser.parse(t2) #otherwise, assumes datetime    
    gap = (t1 - t2).total_seconds()
    return int((gap+1800)/3600)
    
def GetMostRecentDatetime(df, year_col, month_col, day_col, hour_col):
    """Given a data_frame (df) containing a column for year, month, day, and hour
    (year_col), (month_col), (day_col), (hour_col), return the datetime of the final
    row of the dataframe."""
    end_row = df.tail(1); ind = int(end_row.index)
    return dt.datetime(int(end_row.loc[ind][year_col]), int(end_row.loc[ind][month_col]), 
                       int(end_row.loc[ind][day_col]), int(end_row.loc[ind][hour_col]))
    
def ChooseLocalSensors(similar_sensors, northing, easting, zone):
    """Given a data frame of (similar_sensors), the (northing), (easting), and (zone) of a given location,
    use a spatial threshold (dist_thresh) to determine if a sensor is close enough for inclusion""" 
    sub_sensors = similar_sensors[similar_sensors.Zone == zone]
    sub_sensors['LocalDist'] = [math.sqrt((n - northing)**2 + (e - easting)**2) for n,e in zip(sub_sensors.Northing, sub_sensors.Easting)]
    return sub_sensors[sub_sensors.LocalDist < _insitu_proximity_thresh]
    
def SensorIsFresh(local_sensor, temp_thresh, current_time):
    """Return a boolean reflecting whether a local sensor has sufficiently recent data and thus, can be used."""
    last_P_time_gap = (current_time - dateutil.parser.parse(local_sensor['last_P_time'])).total_seconds()
    last_SM_time_gap = (current_time - dateutil.parser.parse(local_sensor['last_SM_time'])).total_seconds() 
    return (True if (last_P_time_gap < temp_thresh and last_SM_time_gap < temp_thresh) else False)     
    
def GetLastSM_and_Time(sensor_df, sm_col, year_col, DOY_col):
    last_ind = sensor_df.index[-1]
    last_year, last_DOY = sensor_df[year_col][last_ind], sensor_df[DOY_col][last_ind]
    return sensor_df[sm_col][last_ind], ConvertYearDOY_to_Datetime(last_year, last_DOY) 
    
def ConvertYearDOY_to_Datetime(year, DOY):
    """Given a (year) and a day-of-year (DOY), return a datetime object, noting the DOY runs from 1-to-366."""
    return RoundTime(dt.datetime(year,1,1,0,0) + dt.timedelta(days = (DOY-1)))
    
def ConvertDatetime_to_year_and_DOY(date_time):
    "Given a (date_time) value, return a year and day-of-year, noting that DOY runs from 1-to-366."""
    date_time = RoundTime(date_time)
    return date_time.year, date_time.timetuple().tm_yday + date_time.hour/24 + date_time.minute/1440
    
def RoundTime(datet, roundTo=3600):
   """Round a datetime object to any time laps in seconds. 
   roundTo : Closest number of seconds to round to, default 1 minute.  """
   seconds = (datet - datet.min).seconds
   rounding = (seconds+roundTo/2) // roundTo * roundTo  # // is a floor division, not a comment on following line:
   return datet + dt.timedelta(0,rounding-seconds,-datet.microsecond)    
    
def GenerateLocalSensorDictionary(similar_local_sensors, path_to_local_sensors, local_dict):
    for network, site_code, northing, easting in zip(similar_local_sensors.Network, similar_local_sensors.Site_Code,
                                                     similar_local_sensors.Northing, similar_local_sensors.Easting):
        name_tag = ((network + "_" + str(site_code)) if network in ['CRN', 'SCAN'] else site_code)
        key = network + '_' + site_code 
        if key not in local_dict.keys(): #then add this new sensor
            local_dict[key] = {}
            sensor_df = pd.read_csv(os.path.join(path_to_local_sensors, network, name_tag + '.csv'))
            local_dict[key]['site_data'] = sensor_df; local_dict[key]['northing'] = northing; local_dict[key]['easting'] = easting
            last_SM, last_time = GetLastSM_and_Time(sensor_df, 'SM', 'Year', 'DOY')
            local_dict[key]['last_SM'] = last_SM; local_dict[key]['last_time'] = last_time        
    return local_dict    