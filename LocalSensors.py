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

_insitu_proximity_thresh = 0.5 #euclidian lat-lon distance required for insitu estimate to be considered
_temporal_thresh = 3600*6 #in seconds...how recent must data be to be deployed?

def GetJSON(f_path, f_name):
	"""Read in a JSON of (f_name) at (f_path)."""
	json_data = open(os.path.join(f_path, f_name)).read()
	return json.loads(json_data)

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
    
def ChooseSimilarSensors(local_sensors, lat, lon, texture_class, topo_class, dist_thresh, temp_thresh):
    """Given a dictionary of (local_sensors), the (lat), (lon), (texture_class), and (topo_class)
    of a given location...note the 'class' variables assume integral values, a (dist_thresh)
    required for a sensor to be close enough for inclusion, and a (temp_thresh) to determine if
    data have gone stale return a list of dictionaries containing the appropriate sensors."""
    similar_texture_map = HydroClass.GetSimilarTexturesMap()
    list_of_similar_sensors, current_time = [], dt.datetime.now()
    current_time = dt.datetime(2016,1,6,20) #for debugging...keeps the data 'fresh' 
    for sensor_code in local_sensors.keys():
        local_sensor = local_sensors[sensor_code] #the sensor to be examined
        if SensorIsSimilar(local_sensor, lat, lon, texture_class, topo_class, dist_thresh, temp_thresh, 
           current_time, similar_texture_map): list_of_similar_sensors.append(local_sensor)
    return list_of_similar_sensors
    
def SensorIsSimilar(local_sensor, lat, lon, texture_class, topo_class, dist_thresh, temp_thresh, 
                    current_time, similar_texture_map):
    """Return a boolean reflecting whether a local sensor is similar and thus, applicable."""
    distance_to_sensor = math.sqrt((lat - local_sensor['lat'])**2 +(lon - local_sensor['lon'])**2)    
    last_P_time_gap = (current_time - dateutil.parser.parse(local_sensor['last_P_time'])).total_seconds()
    last_SM_time_gap = (current_time - dateutil.parser.parse(local_sensor['last_SM_time'])).total_seconds() 
    return (True if (last_P_time_gap < temp_thresh and last_SM_time_gap < temp_thresh and
                    distance_to_sensor < dist_thresh and topo_class == local_sensor['topo_class'] and
                    (texture_class == local_sensor['texture_class'] or local_sensor['texture_class'] in similar_texture_map[texture_class])) else False)     
 
def GetSimilarLocalSensors(path_to_local_sensors, local_sensor_file_name, lat, lon, texture_class, topo_class): 
    local_sensors = GetJSON(path_to_local_sensors, local_sensor_file_name)
    return ChooseSimilarSensors(local_sensors, lat, lon, texture_class, topo_class, _insitu_proximity_thresh, _temporal_thresh) 

if __name__ == "__main__":   
    Env.AddEnvironmentVariables()
    lat, lon, texture_class, topo_class = -32, -110, 3, 1 #for testing only
    #example call to get similar local sensors:
    similar_sensors = GetSimilarLocalSensors(os.environ['path_to_local_sensors'], os.environ['local_sensor_file_name'], 
            lat, lon, texture_class, topo_class)
