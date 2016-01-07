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
    os.environ['local_sensor_file_name'] = "local_sensors.json"
    os.environ['sample_soil_moisture_file'] = "TestSM.csv"
    return None
