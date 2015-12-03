import math

def utmToLatLng(zone, easting, northing, northernHemisphere=True):
	"""Given a zone, easting, and northing, return a latitude and longitude pair.  In the absence of a
	fourth argument, assume the point in question is located in the northern hemisphere.""" 
	
	if not northernHemisphere:
		northing = 10000000 - northing

	a = 6378137
	e = 0.081819191
	e1sq = 0.006739497
	k0 = 0.9996

	arc = northing / k0
	mu = arc / (a * (1 - math.pow(e, 2) / 4.0 - 3 * math.pow(e, 4) / 64.0 - 5 * math.pow(e, 6) / 256.0))

	ei = (1 - math.pow((1 - e * e), (1 / 2.0))) / (1 + math.pow((1 - e * e), (1 / 2.0)))

	ca = 3 * ei / 2 - 27 * math.pow(ei, 3) / 32.0

	cb = 21 * math.pow(ei, 2) / 16 - 55 * math.pow(ei, 4) / 32
	cc = 151 * math.pow(ei, 3) / 96
	cd = 1097 * math.pow(ei, 4) / 512
	phi1 = mu + ca * math.sin(2 * mu) + cb * math.sin(4 * mu) + cc * math.sin(6 * mu) + cd * math.sin(8 * mu)

	n0 = a / math.pow((1 - math.pow((e * math.sin(phi1)), 2)), (1 / 2.0))

	r0 = a * (1 - e * e) / math.pow((1 - math.pow((e * math.sin(phi1)), 2)), (3 / 2.0))
	fact1 = n0 * math.tan(phi1) / r0

	_a1 = 500000 - easting
	dd0 = _a1 / (n0 * k0)
	fact2 = dd0 * dd0 / 2

	t0 = math.pow(math.tan(phi1), 2)
	Q0 = e1sq * math.pow(math.cos(phi1), 2)
	fact3 = (5 + 3 * t0 + 10 * Q0 - 4 * Q0 * Q0 - 9 * e1sq) * math.pow(dd0, 4) / 24

	fact4 = (61 + 90 * t0 + 298 * Q0 + 45 * t0 * t0 - 252 * e1sq - 3 * Q0 * Q0) * math.pow(dd0, 6) / 720

	lof1 = _a1 / (n0 * k0)
	lof2 = (1 + 2 * t0 + Q0) * math.pow(dd0, 3) / 6.0
	lof3 = (5 - 2 * Q0 + 28 * t0 - 3 * math.pow(Q0, 2) + 8 * e1sq + 24 * math.pow(t0, 2)) * math.pow(dd0, 5) / 120
	_a2 = (lof1 - lof2 + lof3) / math.cos(phi1)
	_a3 = _a2 * 180 / math.pi

	latitude = 180 * (phi1 - fact1 * (fact2 + fact3 + fact4)) / math.pi

	if not northernHemisphere:
		latitude = -latitude

	longitude = ((zone > 0) and (6 * zone - 183.0) or 3.0) - _a3

	return (latitude, longitude)
	

#!/usr/bin/env python

# Lat Long - UTM, UTM - Lat Long conversions

from math import pi, sin, cos, tan, sqrt

#LatLong- UTM conversion..h
#definitions for lat/long to UTM and UTM to lat/lng conversions
#include <string.h>

_deg2rad = pi / 180.0
_rad2deg = 180.0 / pi

_EquatorialRadius = 2
_eccentricitySquared = 3

_ellipsoid = [
#  id, Ellipsoid name, Equatorial Radius, square of eccentricity	
# first once is a placeholder only, To allow array indices to match id numbers
	[ -1, "Placeholder", 0, 0],
	[ 1, "Airy", 6377563, 0.00667054],
	[ 2, "Australian National", 6378160, 0.006694542],
	[ 3, "Bessel 1841", 6377397, 0.006674372],
	[ 4, "Bessel 1841 (Nambia] ", 6377484, 0.006674372],
	[ 5, "Clarke 1866", 6378206, 0.006768658],
	[ 6, "Clarke 1880", 6378249, 0.006803511],
	[ 7, "Everest", 6377276, 0.006637847],
	[ 8, "Fischer 1960 (Mercury] ", 6378166, 0.006693422],
	[ 9, "Fischer 1968", 6378150, 0.006693422],
	[ 10, "GRS 1967", 6378160, 0.006694605],
	[ 11, "GRS 1980", 6378137, 0.00669438],
	[ 12, "Helmert 1906", 6378200, 0.006693422],
	[ 13, "Hough", 6378270, 0.00672267],
	[ 14, "International", 6378388, 0.00672267],
	[ 15, "Krassovsky", 6378245, 0.006693422],
	[ 16, "Modified Airy", 6377340, 0.00667054],
	[ 17, "Modified Everest", 6377304, 0.006637847],
	[ 18, "Modified Fischer 1960", 6378155, 0.006693422],
	[ 19, "South American 1969", 6378160, 0.006694542],
	[ 20, "WGS 60", 6378165, 0.006693422],
	[ 21, "WGS 66", 6378145, 0.006694542],
	[ 22, "WGS-72", 6378135, 0.006694318],
	[ 23, "WGS-84", 6378137, 0.00669438]
]

#Reference ellipsoids derived from Peter H. Dana's website- 
#http://www.utexas.edu/depts/grg/gcraft/notes/datum/elist.html
#Department of Geography, University of Texas at Austin
#Internet: pdana@mail.utexas.edu
#3/22/95

#Source
#Defense Mapping Agency. 1987b. DMA Technical Report: Supplement to Department of Defense World Geodetic System
#1984 Technical Report. Part I and II. Washington, DC: Defense Mapping Agency

#def LLtoUTM(int ReferenceEllipsoid, const double Lat, const double Long, 
#			 double &UTMNorthing, double &UTMEasting, char* UTMZone)

def LLtoUTM(ReferenceEllipsoid, Lat, Long):
	"""converts lat/long to UTM coords.  Equations from USGS Bulletin 1532 
	East Longitudes are positive, West longitudes are negative. 
	North latitudes are positive, South latitudes are negative
	Lat and Long are in decimal degrees
	Written by Chuck Gantz- chuck.gantz@globalstar.com
	
	Set the first parameter to 23 for WGS-84... ~EJC"""

	a = _ellipsoid[ReferenceEllipsoid][_EquatorialRadius]
	eccSquared = _ellipsoid[ReferenceEllipsoid][_eccentricitySquared]
	k0 = 0.9996

	#Make sure the longitude is between -180.00 .. 179.9
	LongTemp = (Long+180)-int((Long+180)/360)*360-180 # -180.00 .. 179.9

	LatRad = Lat*_deg2rad
	LongRad = LongTemp*_deg2rad

	ZoneNumber = int((LongTemp + 180)/6) + 1
  
	if Lat >= 56.0 and Lat < 64.0 and LongTemp >= 3.0 and LongTemp < 12.0:
		ZoneNumber = 32

	# Special zones for Svalbard
	if Lat >= 72.0 and Lat < 84.0:
		if  LongTemp >= 0.0  and LongTemp <  9.0:ZoneNumber = 31
		elif LongTemp >= 9.0  and LongTemp < 21.0: ZoneNumber = 33
		elif LongTemp >= 21.0 and LongTemp < 33.0: ZoneNumber = 35
		elif LongTemp >= 33.0 and LongTemp < 42.0: ZoneNumber = 37

	LongOrigin = (ZoneNumber - 1)*6 - 180 + 3 #+3 puts origin in middle of zone
	LongOriginRad = LongOrigin * _deg2rad

    #compute the UTM Zone from the latitude and longitude
	UTMZone = "%d%c" % (ZoneNumber, _UTMLetterDesignator(Lat))

	eccPrimeSquared = (eccSquared)/(1-eccSquared)
	N = a/sqrt(1-eccSquared*sin(LatRad)*sin(LatRad))
	T = tan(LatRad)*tan(LatRad)
	C = eccPrimeSquared*cos(LatRad)*cos(LatRad)
	A = cos(LatRad)*(LongRad-LongOriginRad)

	M = a*((1
            - eccSquared/4
            - 3*eccSquared*eccSquared/64
            - 5*eccSquared*eccSquared*eccSquared/256)*LatRad 
           - (3*eccSquared/8
              + 3*eccSquared*eccSquared/32
              + 45*eccSquared*eccSquared*eccSquared/1024)*sin(2*LatRad)
           + (15*eccSquared*eccSquared/256 + 45*eccSquared*eccSquared*eccSquared/1024)*sin(4*LatRad) 
           - (35*eccSquared*eccSquared*eccSquared/3072)*sin(6*LatRad))
    
	UTMEasting = (k0*N*(A+(1-T+C)*A*A*A/6
                        + (5-18*T+T*T+72*C-58*eccPrimeSquared)*A*A*A*A*A/120)
                  + 500000.0)

	UTMNorthing = (k0*(M+N*tan(LatRad)*(A*A/2+(5-T+9*C+4*C*C)*A*A*A*A/24
                                        + (61
                                           -58*T
                                           +T*T
                                           +600*C
                                           -330*eccPrimeSquared)*A*A*A*A*A*A/720)))

	if Lat < 0:
		UTMNorthing = UTMNorthing + 10000000.0; #10000000 meter offset for southern hemisphere
	return (UTMZone, UTMEasting, UTMNorthing)


def _UTMLetterDesignator(Lat):
	"""This routine determines the correct UTM letter designator for the given latitude
	returns 'Z' if latitude is outside the UTM limits of 84N to 80S
	Written by Chuck Gantz- chuck.gantz@globalstar.com"""

	if 84 >= Lat >= 72: return 'X'
	elif 72 > Lat >= 64: return 'W'
	elif 64 > Lat >= 56: return 'V'
	elif 56 > Lat >= 48: return 'U'
	elif 48 > Lat >= 40: return 'T'
	elif 40 > Lat >= 32: return 'S'
	elif 32 > Lat >= 24: return 'R'
	elif 24 > Lat >= 16: return 'Q'
	elif 16 > Lat >= 8: return 'P'
	elif  8 > Lat >= 0: return 'N'
	elif  0 > Lat >= -8: return 'M'
	elif -8> Lat >= -16: return 'L'
	elif -16 > Lat >= -24: return 'K'
	elif -24 > Lat >= -32: return 'J'
	elif -32 > Lat >= -40: return 'H'
	elif -40 > Lat >= -48: return 'G'
	elif -48 > Lat >= -56: return 'F'
	elif -56 > Lat >= -64: return 'E'
	elif -64 > Lat >= -72: return 'D'
	elif -72 > Lat >= -80: return 'C'
	else: return 'Z'	# if the Latitude is outside the UTM limits

#void UTMtoLL(int ReferenceEllipsoid, const double UTMNorthing, const double UTMEasting, const char* UTMZone,
#			  double& Lat,  double& Long )

def UTMtoLL(ReferenceEllipsoid, northing, easting, zone):
	"""converts UTM coords to lat/long.  Equations from USGS Bulletin 1532 
	East Longitudes are positive, West longitudes are negative. 
	North latitudes are positive, South latitudes are negative
	Lat and Long are in decimal degrees. 
	Written by Chuck Gantz- chuck.gantz@globalstar.com
	Converted to Python by Russ Nelson <nelson@crynwr.com>
	
	Set first parameter to 23 for WGS84 ~EJC"""

	k0 = 0.9996
	a = _ellipsoid[ReferenceEllipsoid][_EquatorialRadius]
	eccSquared = _ellipsoid[ReferenceEllipsoid][_eccentricitySquared]
	e1 = (1-sqrt(1-eccSquared))/(1+sqrt(1-eccSquared))
	#NorthernHemisphere; //1 for northern hemispher, 0 for southern

	x = easting - 500000.0 #remove 500,000 meter offset for longitude
	y = northing

	ZoneLetter = zone[-1]
	ZoneNumber = int(zone[:-1])
	if ZoneLetter >= 'N':
		NorthernHemisphere = 1  # point is in northern hemisphere
	else:
		NorthernHemisphere = 0  # point is in southern hemisphere
		y -= 10000000.0         # remove 10,000,000 meter offset used for southern hemisphere

	LongOrigin = (ZoneNumber - 1)*6 - 180 + 3  # +3 puts origin in middle of zone

	eccPrimeSquared = (eccSquared)/(1-eccSquared)

	M = y / k0
	mu = M/(a*(1-eccSquared/4-3*eccSquared*eccSquared/64-5*eccSquared*eccSquared*eccSquared/256))

	phi1Rad = (mu + (3*e1/2-27*e1*e1*e1/32)*sin(2*mu) 
               + (21*e1*e1/16-55*e1*e1*e1*e1/32)*sin(4*mu)
               +(151*e1*e1*e1/96)*sin(6*mu))
	phi1 = phi1Rad*_rad2deg;

	N1 = a/sqrt(1-eccSquared*sin(phi1Rad)*sin(phi1Rad))
	T1 = tan(phi1Rad)*tan(phi1Rad)
	C1 = eccPrimeSquared*cos(phi1Rad)*cos(phi1Rad)
	R1 = a*(1-eccSquared)/pow(1-eccSquared*sin(phi1Rad)*sin(phi1Rad), 1.5)
	D = x/(N1*k0)

	Lat = phi1Rad - (N1*tan(phi1Rad)/R1)*(D*D/2-(5+3*T1+10*C1-4*C1*C1-9*eccPrimeSquared)*D*D*D*D/24
                                          +(61+90*T1+298*C1+45*T1*T1-252*eccPrimeSquared-3*C1*C1)*D*D*D*D*D*D/720)
	Lat = Lat * _rad2deg

	Long = (D-(1+2*T1+C1)*D*D*D/6+(5-2*C1+28*T1-3*C1*C1+8*eccPrimeSquared+24*T1*T1)
            *D*D*D*D*D/120)/cos(phi1Rad)
	Long = LongOrigin + Long * _rad2deg
	return (Lat, Long)

#if __name__ == '__main__':
#    (z, e, n) = LLtoUTM(23, 40 + (6 + 18.3591/60)/60, -(88 + (13 + 9.52349/60)/60))
#    print z, e, n
#    print UTMtoLL(23, e, n, z)
