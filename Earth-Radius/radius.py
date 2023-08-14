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

"""
Code for Unpacking the Data
"""
data1 = np.loadtxt("grav.txt", skiprows=1, unpack=True)  #Day1 measurements
data2 = np.loadtxt("grav2.txt", skiprows=1, unpack=True) #Day2 measurements

step_count = 22 #number of stair steps per floor.
step_height = 0.176 #each step is 0.176m high.
meter_constant = 0.10023 #gravimeter constant

#Corrections due to Earth rotation from reference measurements
corr1 = (data1[2][0] + data1[2][7]+ data1[2][11])/3
corr2 = (data1[2][11] + data1[2][18]+ data1[2][23])/3
corr3 = (data2[2][0] + data2[2][3]+ data2[2][7])/3

#Day 1, 10AM
elevation1 = (data1[1][:13] - 3) * step_count * step_height
acceleration1 = (data1[2][:13] - corr1) * meter_constant #in mgals
si_acceleration1 = acceleration1 * 10**-5 #in m/s^2

#Day 1, 11AM
elevation2 = (data1[1][12:-2] - 3) * step_count * step_height
acceleration2 = (data1[2][12:-2] - corr2) * meter_constant #in mgals
si_acceleration2 = acceleration2 * 10**-5 #in m/s^2

#Day 2, 11AM
elevation3 = (data2[1][:-3] - 3) * step_count * step_height
acceleration3 = (data2[2][:-3] - corr3) * meter_constant #in mgals
si_acceleration3 = acceleration3 * 10**-5 #in m/s^2

#Total Measurements
elevation = np.concatenate((elevation1, elevation2, elevation3))
acceleration = np.concatenate((acceleration1, acceleration2, acceleration3))
si_acceleration = np.concatenate((si_acceleration1, si_acceleration2, si_acceleration3))

"""
Code for Curve Fitting and Plotting
"""
#Returns gravitational acceleration as a function of elevation.
#Constant gr is reference gravitational acceleration divided by reference Earth radius.
def acc(elev, gr):
    return -2 * gr * elev

#mgal curve fits
p_opt1, p_cov1 = opt.curve_fit(acc, elevation1, acceleration1)
p_opt2, p_cov2 = opt.curve_fit(acc, elevation2, acceleration2)
p_opt3, p_cov3 = opt.curve_fit(acc, elevation3, acceleration3)
p_opt, p_cov = opt.curve_fit(acc, elevation, acceleration)

#m/s^2 curve fits
p_opt1_si, p_cov1_si = opt.curve_fit(acc, elevation1, si_acceleration1)
p_opt2_si, p_cov2_si = opt.curve_fit(acc, elevation2, si_acceleration2)
p_opt3_si, p_cov3_si = opt.curve_fit(acc, elevation3, si_acceleration3)
p_opt_si, p_cov_si = opt.curve_fit(acc, elevation, si_acceleration)

plt.figure(0) #mgal plot
plt.scatter(elevation1, acceleration1, label="Nov. 28, 10AM")
plt.plot(elevation1, acc(elevation1, p_opt1[0]), label="Nov. 28, 10AM")
plt.scatter(elevation2, acceleration2, label="Nov. 28, 11AM")
plt.plot(elevation2, acc(elevation2, p_opt2[0]), label="Nov. 28, 11AM")
plt.scatter(elevation3, acceleration3, label="Dec. 3, 11AM")
plt.plot(elevation3, acc(elevation3, p_opt3[0]), label="Dec. 3, 11AM")

plt.plot(elevation, acc(elevation, p_opt[0]), label="Combined Fit")
plt.title("Grav. Accel. Change (mgal) vs. Elevation Change")
plt.xlabel("Change in Elevation from Reference (m)")
plt.ylabel("Change in Grav. Accel. from Reference (mgal)")
plt.legend()
plt.savefig("mgalplot.pdf")

plt.figure(1) #m/s^2 plot
plt.scatter(elevation1, si_acceleration1, label="Nov. 28, 10AM")
plt.plot(elevation1, acc(elevation1, p_opt1_si[0]), label="Nov. 28, 10AM")
plt.scatter(elevation2, si_acceleration2, label="Nov. 28, 11AM")
plt.plot(elevation2, acc(elevation2, p_opt2_si[0]), label="Nov. 28, 11AM")
plt.scatter(elevation3, si_acceleration3, label="Dec. 3, 11AM")
plt.plot(elevation3, acc(elevation3, p_opt3_si[0]), label="Dec. 3, 11AM")

plt.plot(elevation, acc(elevation, p_opt_si[0]), label="Combined Fit")
plt.title("Grav. Accel. Change (ms^-2) vs. Elevation Change")
plt.xlabel("Change in Elevation from Reference (m)")
plt.ylabel("Change in Grav. Accel. from Reference (ms^-2)")
plt.yticks(rotation=45)
plt.legend()
plt.savefig("ms2plot.pdf")

"""
Code for Calculation and Statistics
"""
ref_grav = 9.804253 #reference value of gravity

#Reduced Chi-Squared tests for every SI plot
print("Chi-Square Test Statistics:")
print(sp.chisquare(acc(elevation1, p_opt1_si[0]), si_acceleration1)[0]/(len(si_acceleration1)-1))
print(sp.chisquare(acc(elevation2, p_opt2_si[0]), si_acceleration2)[0]/(len(si_acceleration2)-1))
print(sp.chisquare(acc(elevation3, p_opt3_si[0]), si_acceleration3)[0]/(len(si_acceleration3)-1))
print(sp.chisquare(acc(elevation, p_opt_si[0]), si_acceleration)[0]/(len(si_acceleration)-1))

#Array of all SI curve fit parameters
param_array = np.array([p_opt1_si[0], p_opt2_si[0], p_opt3_si[0], p_opt_si[0]])

#Calculatin Earth's radius from each SI fitted parameter
param_array = param_array ** -1 * ref_grav
print("Radius of the Earth Resutls:")
for param in param_array:
    print(param)

print("Radius Standard Deviation")
print(stdev(param_array))