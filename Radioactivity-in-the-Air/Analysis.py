# -*- coding: utf-8 -*-
"""
General Code for the Entire Lab
"""
import numpy as np
import scipy.optimize as opt
import scipy.stats as sp
import matplotlib.pyplot as plt

#defining the average and standard deviation functions to be used later
def average(container):
    avg = np.sum(container)/len(container)
    return(avg)

def stdev(container):
    mean_x = average(container)
    #an array of the differences between values of the array and the mean of the array
    avdiff_x = container - mean_x
    stdev_x = np.sqrt((sum(avdiff_x**2))/(len(container)-1))
    #tells the function to return the standard deviation of the array
    return(stdev_x)

#propagation of error functions
def adsubprop(err_x, err_y): #err_x is the error in the argument x
    err_adsub=np.sqrt(err_x**2+err_y**2)
    return(err_adsub)

def muldivprop(x, y, z, err_x, err_y):
    err_muldiv=z*(np.sqrt((err_x/x)**2+(err_y/y)**2))
    return(err_muldiv)

def exprop(exponent, x, err_x):
    err_exp=exponent*(x**exponent)*(err_x/x)
    return(err_exp)

"""
Code for the Cesium-137 Sample
"""
#model function for a logistic curve with plateau
def p(v, k, L):
    return L/(1+np.exp(-k*v))

#importing the data collected from the cesium sample in the Geiger-Muller Tube
cesium_data = np.loadtxt("650-820 volt data.txt", skiprows=2, unpack=True)
cesium_data_av = average(cesium_data[1])
cesium_data_er = np.sqrt(cesium_data_av)#counting error

#calculate the average counts detected for each voltage level
v1 = average(cesium_data[1][0:12])
v2 = average(cesium_data[1][12:24])
v3 = average(cesium_data[1][24:36])
v4 = average(cesium_data[1][36:48])
v5 = average(cesium_data[1][48:60])
v6 = average(cesium_data[1][60:72])
v7 = average(cesium_data[1][72:84])
v8 = average(cesium_data[1][84:93])
v9 = average(cesium_data[1][93:102])
v10 = average(cesium_data[1][102:110])
v11 = average(cesium_data[1][110:118])
v12 = average(cesium_data[1][118:127])

voltages = np.array([650, 660, 670, 680, 700, 710, 720, 740, 760, 780, 800, 820])
cesium_counts = np.array([v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12])

#Python refuses to execute this properly; see discussion.
p_opt_cesium, p_cov_cesium = opt.curve_fit(p, voltages, cesium_counts)

plt.figure(0)
plt.scatter(voltages, cesium_counts)
plt.errorbar(voltages, cesium_counts, yerr= cesium_data_er, fmt='-', ls="none")
#plt.plot(voltages, p(voltages, p_opt_cesium[0], p_opt_cesium[1]), 'm--', label="Linear Regression Fit")
plt.title("Geiger Counter Operating Voltage")
plt.xlabel("Voltage (V)")
plt.ylabel("Mean detected radioactive incidences in 12 second interval")
plt.savefig("CesiumPlateau.pdf")

"""
Code for Unpacking the Air Sample
"""
#importing the data collected from the Air Sample in the Geiger-Muller Tube
air_data = np.loadtxt("AirRadioactivityNov12-1208PM.txt", skiprows=2, unpack=True)
air_data_av = average(air_data[1])
air_data_er = np.sqrt(air_data_av)#counting error

#importing the data collected on the background radiation
background = np.loadtxt("BackgroundNov12-1017AM.txt", skiprows=2, unpack=True)
background_av = average(background[1])
background_er = np.sqrt(background_av)#counting error

"""
Code for the Descendents of Thorium-232
"""
#linear model function
def f(x, a, b):
    return a*x+b

#cleaning the data solely for the long-lived descendents
sample_long = air_data[0][90:]
counts_long = air_data[1][90:]#-background_av #adjusted for background radiation
#the below code accounts for samples where correction for background leads to negative results
counts_long = np.round(counts_long)
for i in range (0, len(counts_long)):
    if counts_long[i] < 1:
        counts_long[i] = 1
counts_long_av = average(counts_long)
counts_long_er = adsubprop(air_data_er, background_er)
lin_counts_long = np.log(counts_long)#counts_long linearised logarithmically
lin_counts_long_er = np.abs(counts_long_er/counts_long_av)#propagation of error

#curve fit parameters for linear regression
p_opt_long , p_cov_long = opt.curve_fit(f, sample_long, lin_counts_long)

#plotting logarithmic linear regression
plt.figure(1)
plt.scatter(sample_long, lin_counts_long, label="Experimental Decay")
plt.plot(sample_long, f(sample_long, p_opt_long[0], p_opt_long[1]),\
         'm--', label="Linear Regression Fit")
plt.errorbar(sample_long, f(sample_long, p_opt_long[0], p_opt_long[1]),\
             yerr=lin_counts_long_er, fmt='r-', ls="none")
plt.title("Linearized Plot of Thorium-232 Descendent Decay")
plt.xlabel("Sample Number")
plt.ylabel("Natural Logarithm of Counts Per 20 Seconds")
plt.legend()
plt.savefig("ThoriumDecay.pdf")

print("The reduced chi-sqaured value for the long-lived isotope model is ")
print(sp.chisquare(f(sample_long, p_opt_long[0], p_opt_long[1]), lin_counts_long)[0]\
      /(len(lin_counts_long)-2))

"""
Code for the Descendents of Uranium-238
"""
#Whyte and Taylor's relation between activity rate and detector efficiency
def wat(t, a, r):
    eb = 0.80#detector efficiency for Lead-214
    ec = 0.95#detector efficiency for Bismuth-214
    lb = np.log(2)/(10.64*60/3)#decay constant for Lead-214 in 20 second intervals.
    lc = np.log(2)/(1.00*60/3)#decay constant for Bismuth-214 in 20 second intervals.
    return a*(((1+(ec/eb)*(lc/(lc-lb)))*np.exp(-lb*t))+((ec/eb)*(r-(lc/(lc-lb)))*np.exp(-lc*t)))

#unpacking the imported Geiger-Muller Tube data
sample_short = air_data[0][:90]
counts_short = air_data[1][:90]-f(sample_short, p_opt_long[0], p_opt_long[1])-background_av
counts_short_av = average(counts_short)
counts_short_er = adsubprop(air_data_er, background_er)
p_opt_short , p_cov_short = opt.curve_fit(wat, sample_short, counts_short)

#plotting logarithmic linear regression
plt.figure(2)
plt.scatter(sample_short, counts_short, label="Experimental Decay")
plt.plot(sample_short, wat(sample_short, p_opt_short[0], p_opt_short[1]),\
         'm--', label="Whyte and Taylor Model Fit")
plt.errorbar(sample_short, wat(sample_short, p_opt_short[0], p_opt_short[1]),\
             yerr=counts_short_er, fmt='r-', ls="none")
plt.title("Whyte and Taylor Model of Uranium-238 Descendent Decay")
plt.xlabel("Sample Number")
plt.ylabel("Counts Per 20 Seconds")
plt.legend()
plt.savefig("UraniumDecay.pdf")

print("The reduced chi-sqaured value for the short-lived isotope model is ")
print(sp.chisquare(wat(sample_short, p_opt_short[0], p_opt_short[1]),\
                   counts_short)[0]/(len(counts_short)-2))