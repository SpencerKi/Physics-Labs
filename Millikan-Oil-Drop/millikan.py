# -*- coding: utf-8 -*-
"""
Code for Analysis of the Millikan Oil-Drop Experiment
"""
import pandas as pd
import numpy as np
import scipy.signal as sig
import matplotlib.pyplot as plt

"""
The following block of code was used to analse the individual excel files 
containing the paths of the oil droplets. ".xlsx" was replaced by the filename 
in question, with X,Y being the observed boundaries on terminal velocity 
movement and A,B being the observed boundaries on upwards movement. All data 
was saved as part of "millkian.txt".

data = pd.read_excel(".xlsx", na_values = "#NV").fillna(0).to_numpy(np.int64).flatten()
time = np.linspace(0, (len(data)-1), len(data))

plt.figure(0)
plt.scatter(time, data)
plt.xlabel("Frame")
plt.ylabel("Pixels from Top")

print(sig.argrelextrema(data, np.less))
print(sig.argrelextrema(data, np.greater))

vt = 10*(data[X]-data[Y])/(520*(time[X]-time[Y]))
v2 = 10*(data[A]-data[B])/(520*(time[A]-time[B]))
print("vt: " + str(vt))
print("v2: " + str(v2))
"""
"""
Unpacking the Data
"""
data = np.loadtxt("millikan.txt", delimiter = ",", skiprows = 1, unpack = True)
stop_volt = data[1]#in volts
up_volt = data[2]#in volts
term_velo = data[3]*10**-3#in m/s
up_velo = data[4]*10**-3#in m/s
C = 2.0232e-10#Calculated constant. No uncertainty as calculated from known values.
volt_err = 5#Reading error of voltage measurements
#Due to discrepancy in size compared to above, velocity reading error is negligble.

"""
Calculation of Q
"""
#Calculations of Q via Method 1
method1 = np.array([])
for i in range(len(term_velo)):
    method1 = np.append(method1, C*term_velo[i]**1.5/stop_volt[i])
method1err = (5/stop_volt)*method1#Propagation of error for method 1.
#Calculations of Q via Method 2
method2 = np.array([])
for j in range(len(term_velo)):
    method2 = np.append(method2, C*(term_velo[j]+up_velo[j])*term_velo[j]**0.5/up_volt[j])
method2err = (5/up_volt)*method2#Propagation of error for method 2.

"""
Calculation of Radius
"""
#All values given in the lab manual.
n = 1.827e-5
g = 9.8
p = 875.3 - 1.204

r = np.sqrt(9*n*term_velo/(2*g*p))#As derived in the report.

np.corrcoef(r, method1)
np.corrcoef(r, method2)
np.corrcoef(r, method1err)
np.corrcoef(r, method2err)
"""
Plots
"""

plt.figure(0)
plt.hist(method1)
plt.xlabel("Q")
plt.ylabel("Occurences")
plt.title("Calculations with Method 1")
plt.savefig("1.jpg")

plt.figure(1)
plt.hist(method2)
plt.xlabel("Q")
plt.ylabel("Occurences")
plt.title("Calculations with Method 2")
plt.savefig("2.jpg")

plt.figure(2)
plt.hist(method1)
plt.hist(method2)
plt.xlabel("Q")
plt.ylabel("Occurences")
plt.title("Comparions Between Both Methods")
plt.savefig("3.jpg")