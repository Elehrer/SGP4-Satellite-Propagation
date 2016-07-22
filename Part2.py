# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 22:50:48 2016

@author: Elliott
"""

import os
import numpy
from numpy import linalg

constellationName = input("Please enter constellation name: ")
file = open("/Users/Elliott/.spyder2-py3/" + "Final %s Results" % (constellationName), 'w')
list = os.listdir("/Users/Elliott/.spyder2-py3")
lineToCopy = 5

array1 = []
satelliteNumber = -1 #AKA "N", set to -1 so on first satellite will be 0
#need to loop through this process

for i in list:
    if (("ENU" in i) and (constellationName in i)): #these are the satellites we want
        satelliteNumber+=1
        file1 = open("/Users/Elliott/.spyder2-py3/" + i)
        print(file1.readline()) #first four lines are header
        print(file1.readline())
        print(file1.readline())
        print(file1.readline())
        
        
        #array1.append([float(x) for x in file1.readline().split()])
        
        firstLine = file1.readline()
        array1.append([float(x) for x in firstLine.split()]) #add file line to array1
        if(array1[satelliteNumber][5]<5): #if the satellite is less than 5 degrees, its info gets written over
            del array1[satelliteNumber]            
            
            file1.close()
            continue
            
                
        print(array1[satelliteNumber][0]) 
        print(array1[satelliteNumber][1])
        print(array1[satelliteNumber][2])
        print(array1[satelliteNumber][3])
        print(array1[satelliteNumber][4])
        print(array1[satelliteNumber][5])
        
        print(array1.__len__())
        file1.close()
        


G = numpy.zeros(shape=(satelliteNumber+1,4)) #tall, wide
while(satelliteNumber>0):
    G[satelliteNumber][0] = array1[satelliteNumber][1] #0 for first satellite, 1 for east
    G[satelliteNumber][1] = array1[satelliteNumber][2] #0 for first satellite, 2 for north
    G[satelliteNumber][2] = array1[satelliteNumber][3] #0 for first satellite, 3 for up
    G[satelliteNumber][3] = -1
    satelliteNumber-=1
    
transposeG = G.transpose()
print(G)
print(transposeG)
dotResult = numpy.dot(transposeG , G)
print(dotResult)
print(linalg.det(dotResult))
print(dotResult[0][0]*dotResult[0][1]*dotResult[0][2]*dotResult[0][3]*dotResult[0][0]**10)
matrixD = linalg.inv(dotResult)
print(matrixD)


file.write("%s\t%s\t%s\t%s\t%s" % (array1[0][0], str(matrixD[0][0]), str(matrixD[1][1]), str(matrixD[2][2]), str(matrixD[3][3])))
file.close()