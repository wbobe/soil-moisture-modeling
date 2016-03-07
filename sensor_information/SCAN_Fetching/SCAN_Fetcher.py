# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 12:57:27 2016
@author: Evan Coopersmith
"""

"""Fetch files from SCAN, generate .csv files for soil moisture calibration."""

import os
import pandas as pd
import numpy as np
import requests
import datetime as dt
import time

_site_codes = pd.read_csv('SCAN_codes.csv')
_SCAN_url = "http://wcc.sc.egov.usda.gov/nwcc/view"
_current_year = dt.datetime.now().year
_SCAN_years = [i for i in range(2001,_current_year + 1)] 
_max_hourly_P = 2 #in inches...anything more in an hour is considered invalid
#### ****HourlyPrecipitation is given in INCHES.  Convert to mm before writing to file."""	
_sm_conversion = 0.01 #convert from 0-100 to 0-1
_p_conversion = 25.4 #convert from inches to mm
	
def Write_Request_To_File(request_data, file_path, site_name, year, month):
	all_text = request_data.text
	output_csv_file = open(os.path.join(file_path, "_".join([site_name, year, month]) + '.csv'), 'w')
	output_csv_file.write(all_text)
	output_csv_file.close()
	return None	
	
def request_SCAN_month_data_at_site(NOAA_url, site_code, year, month):	
	print('Making request to SCAN for site: ', site_code, "Year: ", year, "Month:", month)
	payload = {'intervalType' : 'View Historic', 'report' : 'ALL', 'timeseries' : 'Hourly', 'format' : 'copy', 
                 'sitenum' : site_code, 'interval' : 'MONTH', 'year' : year, 'month' : month, 'day' : '', 'userEmail' : ''}
	return requests.post(NOAA_url, data=payload)
	
def DownloadAllSCAN(site_codes, SCAN_url, SCAN_years):
    for site_code in site_codes.Code:
        output_folder = str(site_code) + '_monthly_files'
        if not os.path.exists(output_folder): os.makedirs(output_folder)
        for year in SCAN_years:
            for month in range(1,13):
                month_str = '0' + str(month) if month < 10 else str(month)
                if not os.path.exists(os.path.join(str(site_code) + '_monthly_files', "_".join([str(site_code), str(year), str(month)]) + '.csv')):
                    time.sleep(10)  # Delay for 10 seconds
                    r = request_SCAN_month_data_at_site(SCAN_url, str(site_code), str(year), month_str)
                    Write_Request_To_File(r, output_folder, str(site_code), str(year), str(month))
    return None

def ProcessAllSites(site_codes, SCAN_years, sm_conversion, p_conversion):
	for site_code in site_codes.Code:
		ProcessSCAN_files_by_site(site_code, SCAN_years, sm_conversion, p_conversion)
	return None
	
def ProcessSCAN_files_by_site(site_code, SCAN_years, sm_conversion, p_conversion):
    SCAN_df = {'SM5':[], 'SM10':[], 'SM20':[], 'P':[], 'Year':[], 'DOY':[]} #depth in cm
    folder_name =  str(site_code) + '_monthly_files'   
    for year in SCAN_years:	
        for month in range(1,13):
            file_path = os.path.join(folder_name, "_".join([str(site_code), str(year), str(month)]) + '.csv')
            if os.path.exists(file_path):
                print("Processing: ", site_code, "month", month, "in year", year)
                SCAN_df = ProcessSCAN_csv_file(SCAN_df, file_path)
    SCAN_df = pd.DataFrame(SCAN_df)
    SCAN_df = CleanSCANData(SCAN_df, sm_conversion, p_conversion)            
    SCAN_df.to_csv(os.path.join("_".join([str(site_code), 'complete_record']) + '.csv'), index = False)
    return None

def CleanSCANData(SCAN_df, sm_conversion, p_conversion):
    for col in SCAN_df.columns:
        if col == 'P':
            SCAN_df['P'] = [p*p_conversion for p in SCAN_df.P]
        elif 'SM' in col:
            SCAN_df[col] = [sm*sm_conversion for sm in SCAN_df[col]]
    return SCAN_df
	
def GetSCANColumns(SCAN_df, SCAN_file_columns):
	matched_columns = [] #the columns from the SCAN file that correspond to the columns of our new dataframe
	for col in SCAN_df.keys():
         if 'SM' in col: #this is a soil moisture column:
             depth = str(int(int(col.split('SM')[-1])/2.5))
             sm_columns = [c for c in SCAN_file_columns if ('SMS' in c and depth in c)]
             matched_columns.append(sm_columns[0] if len(sm_columns) > 0 else 'DoesNotExist')
         elif col == 'P': #the precip column  
             p_interval_columns = [c for c in SCAN_file_columns if 'PRCP' in c]
             if len(p_interval_columns) > 0: 
                 matched_columns.append(p_interval_columns)
             else:
                 p_cumulative_columns = [c for c in SCAN_file_columns if 'PREC' in c]
                 matched_columns.append(p_cumulative_columns if len(p_cumulative_columns) > 0 else 'DoesNotExist')
         elif col == 'DOY':
             matched_columns.append(['Date','Time']) #for calculating the year and DOY...
         else:
             matched_columns.append('DoesNotExist')
	return matched_columns		
 
def ProcessSCAN_csv_file(SCAN_df, file_path):
	month_f = pd.read_csv(file_path, skiprows = 3); previous_P = 0
	if 'Date' not in month_f.columns: month_f = pd.read_csv(file_path, skiprows = 2) #cases of missing whitespace  
	if type(month_f.Date[0]) != str: #empty file, no data...
         return SCAN_df    
	else:
         SCAN_columns = GetSCANColumns(SCAN_df, month_f.columns)
         for i in range(len(month_f)):
             if type(month_f.Date[i]) == str:
                 for df_col, file_col in zip(SCAN_df.keys(), SCAN_columns):
                     if 'Date' in file_col: #fill date and time stamp
                         year, DOY = GenerateYear_and_DOY(month_f['Date'][i], month_f['Time'][i])
                         SCAN_df['Year'].append(year); SCAN_df['DOY'].append(DOY)
                     elif 'SM' in file_col:
                         SCAN_df[df_col].append(month_f[file_col][i])
                     elif 'PR' in file_col or (type(file_col) == list and 'Date' not in file_col):
                         previous_P, P = GetSCANPrecip(month_f, file_col, i, previous_P)
                         SCAN_df[df_col].append(round(P,2) if P < _max_hourly_P else 0)
                     elif 'DoesNotExist' in file_col and df_col not in ['DOY','Year']: #MISSING
                         SCAN_df[df_col].append(-99)
	return SCAN_df
 
def GetSCANPrecip(month_f, p_columns, index, previous_P = 0):
     if 'PRCP' in p_columns[0]:
        return 0, np.mean([month_f[c][index] for c in p_columns if month_f[c][index] > 0]) #average of valid P readings
     else:
        P = np.mean([max(month_f[c][index],0) for c in p_columns if max(month_f[c][index],0) < _max_hourly_P])
        return P, P - previous_P

def GenerateYear_and_DOY(SCAN_date, SCAN_time):
    if '-' in SCAN_date:
        SCAN_datetime = dt.datetime.strptime(SCAN_date + ' ' + SCAN_time, '%Y-%m-%d %H:%M')
    else:
        SCAN_datetime = dt.datetime.strptime(SCAN_date + ' ' + SCAN_time, '%m/%d/%Y %H:%M')        
    return SCAN_datetime.year, SCAN_datetime.timetuple().tm_yday + SCAN_datetime.hour/24 + SCAN_datetime.minute/1440
 
if __name__ == "__main__":
	DownloadAllSCAN(_site_codes, _SCAN_url, _SCAN_years)
	ProcessAllSites(_site_codes, _SCAN_years, _sm_conversion, _p_conversion)
