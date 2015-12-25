# -*- coding: utf-8 -*-
"""
Created on Thu Dec 24 15:59:39 2015
@author: Evan Coopersmith

These modules implement various machine learning tools that are broadly applicable to
soil moisture modelling.
"""
import functools as ft
import numpy as np
import math

_proportion_of_similar_examples = 0.01 #proportion of examples within training set to be used to set k in KNN

def KNN(examples_to_improve, historical_examples, ind_vars, dep_var):
    """This is a standard KNN process, beginning with a data frame of (examples_to_improve) with the KNN 
    process, a data frame of (historical examples) to be used as a training set, a list of independent
    variables (ind_vars) that should correspond to columns in both data frames and a dependent variable,
    (dep_var), which also should correspond to columns within the two dataframes."""
    k = int(len(historical_examples) * _proportion_of_similar_examples) #size of our 'similar set'
    stdevs = [np.std(historical_examples[i]) for i in ind_vars]
    ##This is a fast KNN, for which everything is held constant, apart from the vector we are testing
    pfnx_KNN = ft.partial(RunKNN_on_one_example, historical_examples = historical_examples, ind_vars = ind_vars, dep_var = dep_var, k = k, stdevs = stdevs)
    print("Generating", len(examples_to_improve), "predictions") 
    predictions = examples_to_improve.apply(pfnx_KNN, axis = 1)
    return [p for p in predictions]
    
def RunKNN_on_one_example(single_example, historical_examples, ind_vars, dep_var, k, stdevs):
	"""Given a (single_example) to predict, essentially, one row of the data-frame,
     run the KNN algorithm with (k) examples	and return a prediction"""
	distances = GetDistances(historical_examples[ind_vars], single_example[ind_vars], stdevs)
	historical_examples['Distances'] = distances
	historical_examples = historical_examples.sort('Distances')
	return np.mean(historical_examples[dep_var][0:int(k)])
	
def GetDistances(dataset, single_example, stdevs):
	"""Given an example (x_vec, list) calculate the euclidian distances between this example and all examples within the dataset.	
	The argument, circular (a list of tuples) defines which variables are circular and about which value they reset to zero.
	For instance, times of day are circular between 0 and 1 (time .999 is very similar to time .001, e.g.) while times of year are 
	circular between 0 and 365 (day 363 is very similar to day 4, e.g.)."""
	
	pfnx_euc_dist = ft.partial(FasterEuclidianDistance, single_example = np.array(single_example), stds = stdevs)
	return dataset.apply(pfnx_euc_dist, axis = 1)
 
def FasterEuclidianDistance(old_example, single_example, stds):
	"""Given a single historical example (old_example), and the current example (single_example), 
     determine the Euclidian distance between the two..."""
	return math.sqrt(np.sum(((np.array(old_example) - single_example)/stds)**2))