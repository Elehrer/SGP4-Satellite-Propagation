# -*- coding: utf-8 -*-
"""
Spyder Editor
"""
#everything between here and line 169 is for calculations, just go to main function!
from __future__ import division
from sgp4.earth_gravity import wgs72
from sgp4.io import twoline2rv
from time import gmtime, strftime, mktime, localtime
import time

from numpy import (sin, cos, tan, sqrt, radians, arctan2, hypot, degrees, mod)

class EarthEllipsoid:

    def __init__(self):
        self.a = 6378137.  # semi-major axis [m]
        self.f = 1 / 298.2572235630  # flattening
        self.b = self.a * (1 - self.f)  # semi-minor axis

def ecef2geodetic(x, y=None, z=None, ell=EarthEllipsoid(), deg=True):
    if y is None:
        x, y, z = _depack(x)
    """Algorithm is based on
    http://www.astro.uni.torun.pl/~kb/Papers/geod/Geod-BG.htm
    This algorithm provides a converging solution to the latitude
equation
    in terms of the parametric or reduced latitude form (v)
    This algorithm provides a uniform solution over all latitudes as it
does
    not involve division by cos(phi) or sin(phi)
    """
    ea = ell.a
    eb = ell.b
    rad = hypot(x, y)
# Constant required for Latitude equation
    rho = arctan2(eb * z, ea * rad)
# Constant required for latitude equation
    c = (ea**2 - eb**2) / hypot(ea * rad, eb * z)
# Starter for the Newtons Iteration Method
    vnew = arctan2(ea * z, eb * rad)
# Initializing the parametric latitude
    v = 0
    count = 0
    while (v != vnew).any() and count < 5:
        v = vnew.copy()
#%% Newtons Method for computing iterations
        vnew = v - ((2 * sin(v - rho) - c * sin(2 * v)) /
                    (2 * (cos(v - rho) - c * cos(2 * v))))
        count += 1

#%% Computing latitude from the root of the latitude equation
    lat = arctan2(ea * tan(vnew), eb)
    # by inspection
    lon = arctan2(y, x)

    alt = ((rad - ea * cos(vnew)) * cos(lat)) + \
        ((z - eb * sin(vnew)) * sin(lat))

    if deg:
        return degrees(lat), degrees(lon), alt
    else:
        return lat, lon, alt  # radians

#this is from PySatel and gives same result to EIGHT decimal places
def cbrt(x):
    if x >= 0:
        return pow(x, 1.0/3.0)
    else:
        return -pow(abs(x), 1.0/3.0)

#use this to convert -------------------------------------------------
def ecef2ned(x, y, z, lat0, lon0, h0, ell=EarthEllipsoid(), deg=True):
    e, n, u = ecef2enu(x, y, z, lat0, lon0, h0, ell, deg=deg)
    print("ENU are %f,%f,%f" % (e,n,u))
    return n, e, -u

def ecef2enu(x, y, z, lat0, lon0, h0, ell=EarthEllipsoid(), deg=True):
    x0, y0, z0 = geodetic2ecef(lat0, lon0, h0, ell, deg=deg)
    print("results are %f,%f,%f" % (x0, y0, z0))
    print("results2 are %f,%f,%f" % (x-x0, y-y0, z-z0))
    return _uvw2enu(x - x0, y - y0, z - z0, lat0, lon0, deg=deg)
    
    
def geodetic2ecef(lat, lon, alt, ell=EarthEllipsoid(), deg=True):
    if deg:
        lat = radians(lat)
        lon = radians(lon)
    # radius of curvature of the prime vertical section
    N = get_radius_normal(lat, ell)
    # Compute cartesian (geocentric) coordinates given  (curvilinear) geodetic
    # coordinates.
    x = (N + alt) * cos(lat) * cos(lon)
    y = (N + alt) * cos(lat) * sin(lon)
    z = (N * (ell.b / ell.a)**2 + alt) * sin(lat)
    return x, y, z

def _uvw2enu(u, v, w, lat0, lon0, deg):
    if deg:
        lat0 = radians(lat0)
        lon0 = radians(lon0)
    t = cos(lon0) * u + sin(lon0) * v
    East = -sin(lon0) * u + cos(lon0) * v
    Up = cos(lat0) * t + sin(lat0) * w
    North = -sin(lat0) * t + cos(lat0) * w
    return East, North, Up

def geodetic2enu(lat, lon, h, lat0, lon0, h0, ell=EarthEllipsoid(), deg=True):
    x1, y1, z1 = geodetic2ecef(lat, lon, h, ell, deg=deg)
    x2, y2, z2 = geodetic2ecef(lat0, lon0, h0, ell, deg=deg)
    dx = x1 - x2
    dy = y1 - y2
    dz = z1 - z2
    return _uvw2enu(dx, dy, dz, lat0, lon0, deg=deg)


def geodetic2ned(lat, lon, h, lat0, lon0, h0, ell=EarthEllipsoid(), deg=True):
    e, n, u = geodetic2enu(lat, lon, h, lat0, lon0, h0, ell, deg=deg)
    return n, e, -u


def ned2aer(n, e, d, deg=True):
    return enu2aer(e, n, -d, deg=deg)

"""
def ned2ecef(n, e, d, lat0, lon0, h0, ell=EarthEllipsoid(), deg=True):
    return enu2ecef(e, n, -d, lat0, lon0, h0, ell, deg=deg)


def ned2geodetic(n, e, d, lat0, lon0, h0, ell=EarthEllipsoid(), deg=True):
    x, y, z = enu2ecef(e, n, -d, lat0, lon0, h0,  ell, deg=deg)
    return ecef2geodetic(x, y, z, ell, deg=deg)
"""

def get_radius_normal(lat_radians, ell):
    a = ell.a
    b = ell.b
    return a**2 / sqrt(
        a**2 * (cos(lat_radians))**2 + b**2 *
        (sin(lat_radians))**2)


def _depack(x0):
    if x0.ndim > 2:
        raise TypeError('I expect Nx3 or 3XN triplets')
    m, n = x0.shape
    if m == 3:  # 3xN triplets
        x = x0[0, :]
        y = x0[1, :]
        z = x0[2, :]
    elif n == 3:  # Nx3 triplets
        x = x0[:, 0]
        y = x0[:, 1]
        z = x0[:, 2]
    else:
        raise TypeError('I expect an Nx3 or 3xN input of x,y,z')
    return x, y, z
    
def enu2aer(e, n, u, deg=True):
    r = hypot(e, n)
    slantRange = hypot(r, u)
    elev = arctan2(u, r)
    az = mod(arctan2(e, n), 2 * arctan2(0, -1))
    if deg:
        return degrees(az), degrees(elev), slantRange
    else:
        return az, elev, slantRange  # radians    

#--------------------------------------------------------

def getTime(timeType = gmtime()):
    currentTime = strftime("%Y %m %d %H %M %S", timeType)
    year = int(currentTime[0:4])
    month = int(currentTime[5:7])
    day = int(currentTime[8:10])
    hour = int(currentTime[11:13])
    minute = int(currentTime[14:16])
    second = int(currentTime[17:19])
    curEpoch = int(mktime(localtime()))
    return year, month, day, hour, minute, second, curEpoch

def getInputs():
    header = input('Please enter the satellite title: ') #header is whatever text the user wants at the top of their file
    line1 = input('Please enter line one of the TLE: ')
    line2 = input('Please enter line two of the TLE : ')
    
    decimationRate = int(input('Please enter decimation rate (in seconds): '))

    startTime = int(input('Please enter start time (in terms of Epoch): '))
    stopTime = int(input('Please enter stop time (in terms of Epoch): '))

    return header, line1,line2,decimationRate,startTime,stopTime

def writeToFiles(header, startTime, stopTime, decimationRate, satellite, lat=33.0754490, lon = -117.2492680, alt = 68.25): #these lattitude, longitude, and altitude are just test coordinates
    
    epochToCheck = startTime - 1    
    
    year, month, day, hour, minute, second, epoch = getTime()    
    filename = ("Results For " + header)
    #toOpenFileAsAppend = 'a'g
    file = open("XYZ " + filename, 'w')  #opens XYZ file
    file2 = open("ENU " + filename, 'w')  #opens ENU file
    
    #print(int(mktime(localtime())) >= startTime)
    #print(int(mktime(localtime())) <= stopTime)

    #all header lines written to files start with "%" so they are easier to import into matlab
    file.write("% XYZ Results For " + header + "\n")  #writes title of document
    file.write("%% Run on %s/%s/%s \n" % (month, day, year))  #writes date to document
    file.write("% Epoch \t X-Coord \t Y-Coord \t Z-Coord \n")  #labelling columns: Epoch, X, Y, Z

    file2.write("% ENU Results For " + header + "\n")  #writes title of document, looks like "%ENU Results For GPS BIIA-10 (PRN 32)"
    file2.write("%% Latitude, Longitude, and Altitude of User are %f, %f, and %f, respectively\n" % (lat, lon, alt))    
    file2.write("%% Run on %s/%s/%s \n" % (month, day, year))  #writes date to document, looks like "%Run on 7/20/16"    
    file2.write("% Epoch \tEast \t North \tUp \tAzimuth \tElevation \n")  #labelling columns: Epoch, E, N, U, Az, Ele. Looks like "% Epoch 		 North 		 East 		 Up 		 Azimuth 		 Elevation"
    
    while(epochToCheck <= stopTime):
        epochToCheck+=1
        year, month, day, hour, minute, second = getTime(time.localTime(epochToCheck))
        
        
        position, velocity = satellite.propagate(year,month,day,hour,minute,second)
        
        if (epochToCheck % decimationRate == 0):
            print("Recording Data!")
            file.write(str(epochToCheck) + "\t")

            file.write(str(position[0]) + "\t")

            file.write(str(position[1]) + "\t")

            file.write(str(position[2]) + "\n")

            n, e, u = ecef2ned(position[0], position[1], position[2], lat, lon, alt)
            length = sqrt(n**2+e**2+u**2)
            file2.write(str(epochToCheck) + "\t") #writes epoch
            file2.write(str(e/length) + "\t")  #writes east unit vector
            file2.write(str(n/length) + "\t")  #writes north unit vector
            file2.write(str(u/length) + "\t")  #writes up unit vector
            
            az, elev, slantRange = enu2aer(e, n, u)  #azimuth and elevation both in degrees          
            
            file2.write(str(az) + "\t") #writes azimuth
            file2.write(str(elev) + "\n") #writes elevation
            
            
        epochToCheck+=1
    file.close()
    file2.close()
    print("Finished")

def main():
    fileHeader, TLEline1, TLEline2, decimationRate, startTime, stopTime = getInputs() #first, get all input values from user
    satellite = twoline2rv(TLEline1,TLEline2,wgs72) #next, create an sgp4 satellite based on inputted TLEs and the wgs72 gravitational model 
    writeToFiles(fileHeader, startTime, stopTime, decimationRate, satellite) #last, writes two files, one for XYZ coordinates during epoch range and one for ENU (plus azimuth and elevation) during the time range

main()  
