# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 09:47:08 2015
@author: Evan Coopersmith

Functions herein allow calibration and validation of soil moisture models
"""
import os
import pandas as pd
import numpy as np
import math
from scipy.stats.stats import pearsonr
import GeneticAlgorithms_for_SM as GA

_param_file_name = "ValidSites_USCRN_SCAN_ARS.csv"
"""These values are default assumptions required to calibrate system memory"""
_default_triad = [0.015, 0, 0] #for estimation of system memory, a flat evapotranspiration rate
_long_beta_time = 2000 #maximum possible number of hours to consider precip history
_interval = 50 #incremental number of hours to consider when calibrating a sys memory
_thresh = 0.97 #how similar must the shorter beta series be to the longer one to be acceptable?
"""These values are for calibration of the sinusoidal loss function"""
_grow_start, _grow_end = 100, 300 #days of year, between which we'll attempt to estimate SM


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
    print("Generating 'beta series'"); BetaSeries = GetBeta(SiteData[precip_col], EtaSeries, z, parameters[-1])		
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

def GetSimilarSensors(parameter_file, lat, lon, hydroclimatic_class, topo_class, texture_class):
    """Given a (parameter_file), the relevant (lat) and (lon) where estimates are required, the (hydroclimatic_class),
    (topo_class), and (texture_class), return the indices of parameter_file for site whose parameters can be used to
    generate estimates.  Additionally, return a weight associated with each set of parameters."""        
    #THIS FUNCTION HAS NOT YET BEEN WRITTEN
    return sensor_indices, sensor_weights  

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
    
if __name__ == "__main__":    
    ###Example soil moisture estimation procedure
    parameter_file = GetFlatFile("sensor_information", _param_file_name)
    parameter_list = GetParameters(parameter_file, 23) #random row, number 23
    SiteData = GetFlatFile('sensor_information', 'TestSM.csv')[3000:10000]
    SiteData = CalculateSMEstimates(SiteData, parameter_list, 5, 'SM', 'Precip', 'DOY')
    
    ###Example calibration of soil moisture parameters
    SiteData['Invalid'] = GenerateInvalidSeq(SiteData, 'SM', 'Precip', 'DOY')
    SiteData['BetaSeries'] = [0 for i in range(len(SiteData))]
    parameter_list = CalibrateSMParameters(SiteData, 5, 'SM', 'Precip', 'DOY')