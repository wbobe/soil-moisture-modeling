# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 19:31:35 2015

These modules are designed to receive, as inputs, NLDAS precipitation, USGS
streamflow, and potentially, estimates of potential evapotranspiration,
and return relevant information for hydroclimatic classification.

@author: Evan Coopersmith
"""

from dateutil.parser import parse
import datetime as dt
import numpy as np
import pandas as pd
import os
import Data

_1st_quadrant = {'S1':{'Var': 'max_day_q', 'Cut': 95, True: 'S2', False: 'S3'},
                 'S2':{'Var': 'max_day_p', 'Cut': 91, True: 'LWC', False: 'LPC'},
                 'S3':{'Var': 'max_day_q', 'Cut': 135, True: 'LPM', False: 'S4'},
                 'S4':{'Var': 'max_day_p', 'Cut': 54, True: 'S5', False: 'LPQ'},
                 'S5':{'Var': 'arid_ind', 'Cut': 1, True: 'LBMH', False: 'LBMS'}}
_2nd_quadrant = {'S1': {'Var': 'seasonality', 'Cut': 9999, True: 'LJ', False: 'LJ'}}
_3rd_quadrant = {'S1':{'Var': 'arid_ind', 'Cut': 1.0029, True: 'S2', False: 'S3'},
                 'S2':{'Var': 'max_day_q', 'Cut': 160, True: 'S4', False: 'S5'},
                 'S3':{'Var': 'max_day_q', 'Cut': 120, True: 'S6', False: 'S7'},
                 'S4':{'Var': 'max_day_q', 'Cut': 58, True: 'XHD', False: 'S8'},
                 'S5':{'Var': 'arid_ind', 'Cut': 0.5, True: 'IVD', False: 'ITF'},
                 'S6':{'Var': 'arid_ind', 'Cut': 1.25, True: 'ITC', False: 'S9'},
                 'S7':{'Var': 'max_day_p', 'Cut': 45, True: 'XSMB', False: 'ISQJ'},
                 'S8':{'Var': 'seasonality', 'Cut': 0.3857, True: 'S10', False: 'S11'},
                 'S9':{'Var': 'seasonality', 'Cut': 0.5236, True: 'S12', False: 'XSC'},
                 'S10':{'Var': 'max_day_q', 'Cut': 105, True: 'ITC', False: 'IHM'},
                 'S11':{'Var': 'arid_ind', 'Cut': .4684, True: 'XVM', False: 'XTM'},
                 'S12':{'Var': 'max_day_p', 'Cut': 200, True: 'ISCJ', False: 'ISCB'}}
_4th_quadrant = {'S1':{'Var': 'max_day_q', 'Cut': 103, True: 'S2', False: 'S3'},
                 'S2':{'Var': 'max_day_p', 'Cut': 52, True: 'XADB', False: 'S4'},
                 'S3':{'Var': 'max_day_q', 'Cut': 219, True: 'IAQ', False: 'IAF'},
                 'S4':{'Var': 'seasonality', 'Cut': .4581, True: 'IACJ', False: 'XACJ'}}    

_quadrant_splits = {'S1': {'Var': 'seasonality', 'Cut': 0.2564, True: 'S2', False: 'S3'},
                    'S2': {'Var': 'max_day_p', 'Cut': 152, True: _1st_quadrant, False: _2nd_quadrant}, 
                    'S3': {'Var': 'arid_ind', 'Cut': 1.9171, True: _3rd_quadrant, False: _4th_quadrant}}

_NLDAS_columns = ['Date','Lat','Lon','Precip','ShortWaveFlux','LongWaveFlux']
_USGS_columns = ['USGS_code','Date','Q']

def GetSimilarSensorMap():
    """When calibrated models within the particular class are unavailable, certain
    classes are sufficiently similar to warrant inclusion.  This function returns a hard-coded
    dictionary of which classes are pseudo-substitutes.  The classes themselves are the keys."""
    return {'LWC': ['LPC','LJ'], 'LPC': ['LWC','LJ','LPM'], 'LPM':['LPC','LBMH','LBMS'],
            'LBMH': ['LBMS','LPM','LPQ','IHM'], 'LBMS':['LBMH','LPQ','LPM'], 'LPQ':['LBMH','LBMS','LPM'],
            'LJ':['LPC','LWC','LPM'], 'XHD':['IVD','XVM'], 'ITC':['ISQJ','ISCJ','ISCB'],
            'IHM':['XTM','XVM','LBMH'], 'XVM':['XTM','IHM','IVD'], 'XTM':['XVM','IHM'],
            'IVD':['XVM','XHD'], 'ITF':[], 'ISCJ':['ISCB','ISQJ','ITC'], 'ISCB':['ISCJ','ITC','XSC'],
            'XSC':['ISCB','ISCJ','XSMB'], 'XSMB':['XSC','XADB'], 'ISQJ':['ITC','ISCJ','IACJ'],
            'XADB':['XSMB'], 'IACJ':['XACJ','ISCJ','ISQJ','IAF'], 'XACJ':['IACJ','XADB'], 'IAQ':['IACJ','IAF'], 'IAF':['IAQ','IACJ']}
            
def GetSimilarTexturesMap():
    return {1:[2,3], 2:[1,3], 3:[1,2,6,7], 4:[5,6,8,9], 5:[4], 6:[3,4,7,8,9], 7:[3,6,9,10],
            8:[4,8,9], 9:[4,6,7,8,10,12], 10:[7,9,12], 11:[8,9,12], 12:[7,8,9], 13:[]}
            
def GetSimilarTopoMap():
    return {0:[1], 1:[0,2], 2:[1,3], 3:[2,4], 4:[3]}

def GetUSGSFlat():
    """Returns a flat USGS file containing streamflow time series data"""
    Data.get_bounded_data(os.environ.get('USGS_TABLE_NAME'), bounds)

def GetNLDASFlat():
    """Returns a flat NLDAS file containing precipitation and potential evap."""
    Data.get_bounded_data(os.environ.get('NLDAS_TABLE_NAME'), bounds)
    
def GenerateAverageByDayOfYear(df, date_col, variable_col):
    """Given a dataframe (df), a (date_col) containing date information in the format
    '2013-01-12 22:00:00-05', average (variable_col) by day of year, returning a 365 item
    list containing the average on each of those days of year."""
    days_of_year = [min(parse(d).timetuple().tm_yday, 365) for d in df[date_col]] #force leap years to 365
    daily_averages = [(0,0) for i in range(365)] #tuple contains (ave_val, num_count)
    for day, value in zip(days_of_year, df[variable_col]):
        day_index = day - 1 #indexing 0-364
        ave_val, num_elements = daily_averages[day_index]
        daily_averages[day-1] = ((ave_val * num_elements + value)/(num_elements+1),num_elements+1)
    return [daily_val[0] for daily_val in daily_averages]
    
def CalculateSeasonality(precip_time_series):
    """Calculate seasonality from a precipitation time series according to
    methodology from Walsh & Lawler (1981), modified by Coopersmith et al (2012)"""
    tot_precip, mean_precip = np.sum(precip_time_series), np.mean(precip_time_series)
    return sum([abs(p - mean_precip)/tot_precip for p in precip_time_series])
    
def GetMaxIndex(time_series):
    """Given a (time_series) of indeterminate length (probably 365 in most cases),
    return the index where the maximum value occurs"""
    max_val = np.max(time_series)
    return list(time_series).index(max_val)
    
def EstimatePE(radiation_time_series, temperature_time_series):
    """Estimate the total potential evapotranspiration as described by 
    Oudin (2010), Eq. 3.  Note, (radiation_time_series) is to be given in
    daily MJ/m^2, (temperature_time_series) is the mean 2m air-temp in deg.C"""  
    return np.sum([(0.408*re*(t+5)/100 if t > 5 else 0) for t, re in zip(temperature_time_series, radiation_time_series)])
    
def SmoothCircularTimeSeries(time_series, window):
    """Given a (time_series) of a circular nature (i.e. the average precipitation on day
    365 is adjacent to average precipitation on day 1), smooth to the average value within 
    the (window) of indicies specified."""
    smoothed_values, L, half_window = [], len(time_series), int(window/2)
    for index in range(L):
        if index < half_window: #we are too close to the beginning of the list
            smoothed_values.append(np.mean(time_series[:index] + 
            time_series[(L-(half_window-index)):] + time_series[index:(index+half_window)]))
        elif index > (L-half_window):
            smoothed_values.append(np.mean(time_series[(index-half_window):index] + 
            time_series[0:(index-L+half_window)] + time_series[index:L]))
        else:
            smoothed_values.append(np.mean(time_series[(index-half_window):(index+half_window)]))
    return smoothed_values        

def DetermineClassification(class_tree, access_key, seasonality, arid_ind, max_day_p, max_day_q):
    """Given the four relevant indicators, return the classifications according to the tree
    defined in (Coopersmith et al, 2012)."""
    split_var = class_tree[access_key]['Var']
    if split_var == 'seasonality': split_value = seasonality
    if split_var == 'arid_ind': split_value = arid_ind
    if split_var == 'max_day_p': split_value = max_day_p
    if split_var == 'max_day_q': split_value = max_day_q    
    boolean_split = True if split_value <= class_tree[access_key]['Cut'] else False
    next_split = class_tree[access_key][boolean_split]
    if type(next_split) == dict:
       return DetermineClassification(next_split, 'S1', seasonality, arid_ind, max_day_p, max_day_q) 
    elif next_split[0] == 'S':
       return DetermineClassification(class_tree, next_split, seasonality, arid_ind, max_day_p, max_day_q)         
    else:
       return next_split
       
def CalculateFourIndicators(nldas2_box, usgs_box):
    smoothed_daily_precip = SmoothCircularTimeSeries(GenerateAverageByDayOfYear(nldas2_box, 'Date', 'Precip'),30)
    seasonality, max_day_p = CalculateSeasonality(smoothed_daily_precip), GetMaxIndex(smoothed_daily_precip)         
    arid_ind = 1.2 #NOT IMPLEMENTED.
    smoothed_daily_runoff = SmoothCircularTimeSeries(GenerateAverageByDayOfYear(usgs_box, 'Date', 'Q'),30)     
    max_day_q = GetMaxIndex(smoothed_daily_runoff) 
    return seasonality, arid_ind, max_day_p, max_day_q
    
# Example call to getData
bounds = {'lat_min':'-81', 'lat_max':'-82', 'lon_min':'-81', 'lon_max':'-82'}
#nldas2_box = Data.get_bounded_data(os.environ.get('USGS_TABLE_NAME'), bounds)
nldas2_box = pd.read_csv(os.path.join('test_data','nldas2_hourly_2015_10_21_07_32_39.csv'), 
                        header = None)
nldas2_box.columns = _NLDAS_columns

usgs_box = pd.read_csv(os.path.join('test_data','usgs_gauge_2015_10_21_16_19_53.csv'),
                       header = None)
usgs_box.columns = _USGS_columns
seasonality, arid_ind, max_day_p, max_day_q = CalculateFourIndicators(nldas2_box, usgs_box)
site_class = DetermineClassification(_quadrant_splits, 'S1', seasonality, arid_ind, max_day_p, max_day_q)