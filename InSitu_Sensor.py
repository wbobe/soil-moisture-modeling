# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 16:39:35 2016
@author: Evan Coopersmith

Defines the class 'Sensor', used for local in-situ sensors rather than USCRN gauges or 
SCAN gauges...these would be the sensors installed in a given parcel for personal use.
"""

class Sensor(object):
	"""Creates a sensor object containing the relevant features and the appropriate data streams
	attached to that particular sensor."""    

	def __init__(self, lat, lon, texture_class, topo_class, precip_hist, most_recent_SM, 
              last_valid_time):
