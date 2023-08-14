# -*- coding: utf-8 -*-
"""
Code for the Cavendish Experiment Lab
"""
import numpy as np
import scipy.optimize as opt
import scipy.stats as sta
import matplotlib.pyplot as plt

"""
Unpacking Data
"""
noweights = np.loadtxt("noweights.txt", dtype = str, delimiter = "\t",\
                       skiprows = 2, unpack = True)
nwtime = np.array([])
for t in noweights[0]:
    tcorrec = float(noweights[0][0].split(":")[0])*3600 +\
    float(noweights[0][0].split(":")[1])*60 + float(noweights[0][0].split(":")[2])
    nwtime = np.append(nwtime, float(t.split(":")[0])*3600 +\
                       float(t.split(":")[1])*60 + float(t.split(":")[2]) - tcorrec)
nwposition = np.array([])
for p in noweights[1]:
    nwposition = np.append(nwposition, float(p))
    
ccweights = np.loadtxt("ccweights.txt", dtype = str, delimiter = "\t",\
                       skiprows = 2, unpack = True)
cctime = np.array([])
for t in ccweights[0]:
    tcorrec = float(ccweights[0][0].split(":")[0])*3600 +\
    float(ccweights[0][0].split(":")[1])*60 + float(ccweights[0][0].split(":")[2])
    cctime = np.append(cctime, float(t.split(":")[0])*3600 +\
                       float(t.split(":")[1])*60 + float(t.split(":")[2]) - tcorrec)
ccposition = np.array([])
for p in ccweights[1]:
    ccposition = np.append(ccposition, float(p))
    
noweights2 = np.loadtxt("noweights2.txt", dtype = str, delimiter = "\t",\
                        skiprows = 2, unpack = True)
nw2time = np.array([])
for t in noweights2[0]:
    tcorrec = float(noweights2[0][0].split(":")[0])*3600 +\
    float(noweights2[0][0].split(":")[1])*60 + float(noweights2[0][0].split(":")[2])
    nw2time = np.append(nw2time, float(t.split(":")[0])*3600 +\
                        float(t.split(":")[1])*60 + float(t.split(":")[2]) - tcorrec)
nw2position = np.array([])
for p in noweights2[1]:
    nw2position = np.append(nw2position, float(p))
    
cweights2 = np.loadtxt("cweights2.txt", dtype = str, delimiter = "\t",\
                       skiprows = 2, unpack = True)
c2time = np.array([])
for t in cweights2[0]:
    tcorrec = float(cweights2[0][0].split(":")[0])*3600 +\
    float(cweights2[0][0].split(":")[1])*60 + float(cweights2[0][0].split(":")[2])
    c2time = np.append(c2time, float(t.split(":")[0])*3600 +\
                       float(t.split(":")[1])*60 + float(t.split(":")[2]) - tcorrec)
c2position = np.array([])
for p in cweights2[1]:
    c2position = np.append(c2position, float(p))
    
ccweights2 = np.loadtxt("ccweights2.txt", dtype = str, delimiter = "\t",\
                        skiprows = 2, unpack = True)
cc2time = np.array([])
for t in ccweights2[0]:
    tcorrec = float(ccweights2[0][0].split(":")[0])*3600 +\
    float(ccweights2[0][0].split(":")[1])*60 + float(ccweights2[0][0].split(":")[2])
    cc2time = np.append(cc2time, float(t.split(":")[0])*3600 +\
                        float(t.split(":")[1])*60 + float(t.split(":")[2]) - tcorrec)
cc2position = np.array([])
for p in ccweights2[1]:
    cc2position = np.append(cc2position, float(p))
    
"""
Calculations of G
"""
#Period calculations
nwT = (nwtime[np.where(nwposition == 0.248)[0][-1]] -\
              nwtime[np.where(nwposition == 0.244)[0][0]])/11
ccT = (cctime[np.where(ccposition == 0.224)[0][-1]] -\
              cctime[np.where(ccposition == 0.216)[0][0]])/10
nw2T = (nw2time[np.where(nw2position == 0.243)[0][-1]] -\
                nw2time[np.where(nw2position == 0.238)[0][0]])/10
c2T = (c2time[np.where(c2position == 0.253)[0][-1]] -\
              c2time[np.where(c2position == 0.304)[0][0]])/11
cc2T = (cc2time[np.where(cc2position == 0.219)[0][-1]] -\
                cc2time[np.where(cc2position == 0.241)[0][0]])/8

#Calculation of G
b = 0.06
d = 0.05
S = np.abs(np.mean(nwposition) - np.mean(cc2position))
m = 1.5
T = cc2T
L = 4.586
G = (np.pi**2*b**2*d*S)/(m*T**2*L)

#Error calculations
be = 0.0005
Se = np.sqrt(2*0.0000005**2)
Te = 0.0005
Le = 0.0005
error = ((S*b**2)/(L*T**2))*np.sqrt(S*b**2*((Se/S)**2 + (2*be/b)**2) +\
         L*T**2*((Le/L)**2 + (2*Te/T)**2))

"""
Calculations for Questions
"""
#Symmetry of E
cS = np.abs(np.mean(nwposition) - np.mean(c2position))
ccS = np.abs(np.mean(nwposition) - np.mean(cc2position))

#Factor of other weight
beta = b**3/(b**2 + 4*d**2)**1.5
betae = beta*np.sqrt((9*be**2/b**2) + ((9*b**2*be**2)/(b**2 + 4*d**2)**2))

"""
Plots
"""
#No Weights - Day 1
plt.figure(0)
plt.plot(nwtime, nwposition)
plt.xlabel("Time (s)")
plt.ylabel("Position (m)")
plt.title("No Large Masses")
plt.savefig("noweights.pdf")

#Counterclockwise Weights - Day 1
plt.figure(1)
plt.plot(cctime, ccposition)
plt.xlabel("Time (s)")
plt.ylabel("Position (m)")
plt.title("Large Masses in Counterclockwise Configuration")

#No Weights - Day 2
plt.figure(2)
plt.plot(nw2time, nw2position)
plt.xlabel("Time (s)")
plt.ylabel("Position (m)")
plt.title("No Large Masses")

#Clockwise Weights - Day 2
plt.figure(3)
plt.plot(c2time, c2position)
plt.xlabel("Time (s)")
plt.ylabel("Position (m)")
plt.title("Large Masses in Clockwise Configuration")
plt.savefig("cweights.pdf")

#Counterclockwise Weights - Day 2
plt.figure(4)
plt.plot(cc2time, cc2position)
plt.xlabel("Time (s)")
plt.ylabel("Position (m)")
plt.title("Large Masses in Counterclockwise Configuration")
plt.savefig("ccweights.pdf")