# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 09:25:14 2016
@author: Evan Coopersmith

Environment variable .json containing information regarding local directories
"""
import os

def AddEnvironmentVariables():
    os.environ['param_file_name'] = "ValidSites_USCRN_SCAN_ARS.csv"
    os.environ['path_to_local_sensors'] = os.path.join('sensor_information')
    os.environ['path_to_topos_and_textures'] = os.path.join('sensor_information','Topos_and_Textures')
    os.environ['path_to_historical_requests'] = os.path.join('test_data', 'HistoricalAnalyses')
    os.environ['local_sensor_file_name'] = "local_sensors.json"
    os.environ['sample_soil_moisture_file'] = "TestSM.csv"
    return None
