import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sp
import scipy.stats as sta

#reusing my average function from before
def average(container):
    avg = np.sum(container)/len(container)
    return(avg)

#propagation of error modules
def adsubprop(err_x, err_y): #err_x is the error in the argument x
    err_adsub=np.sqrt(err_x**2+err_y**2)
    return(err_adsub)

def muldivprop(x, y, z, err_x, err_y):#z is the final result of multiplication/division of x and y
    err_muldiv=z*(np.sqrt((err_x/x)**2+(err_y/y)**2))
    return(err_muldiv)

def exprop(exponent, x, err_x):
    err_exp=exponent*(x**exponent)*(err_x/x)
    return(err_exp)

#defines a model function for the curve fit to use
def model_function(x, a, b):
    return a*x+b

#defining constants
perm =  4e-7*np.pi #universal constant
turns = 130 #given constant
R = 15.85*0.01 #radius of the coil
K = ((4/5)**(3/2))*perm*turns/R #constant; B-field/current (coil constant)
k = K/np.sqrt(2) #coil constant

#loading data
data = np.loadtxt('Fixed-Voltage-Charge-to-Mass.txt', skiprows = 1, unpack = True)
vdata = data[0] #
idata = data[1] #amperes
radius = (data[2]/2)*0.01 #radius of the electron beam path
inv_r = 1/radius #1/radius

#error for V is recorded as the standard deviation in V
#error for I is the reading error
#reading error for the ruler was recorded as +/- 0.5mm = 0.005m
#this value was used as the instrument error was given at 0.04%, an order of
#magnitude below reading error.
err_i = data[3]
err_v = np.std(vdata)
err_r = 0.0005

B = ((4/5)**(3/2))*((perm*turns*idata)/radius)

#curve fitting of 1/r as a linear function of current by Equation 7
p_opt, p_cov = sp.curve_fit(model_function, idata, inv_r,(1,0), err_i, True)

#Derived in discussion
Bext = p_opt[1]*K/p_opt[0]
var_invr = p_cov[0,1]

#plotting
plt.figure()
plt.plot(np.arange(-0.1,2.6,0.1), p_opt[0]*np.arange(-0.1,2.6,0.1)+p_opt[1],\
         label = 'Line of Best Fit')
plt.grid(True)
plt.scatter(idata, inv_r, label = "Experimental Data")
plt.errorbar(idata, inv_r, yerr = var_invr, ls = 'none', label = 'Variance')
plt.ylabel("Inverse of Radius (1/M)")
plt.xlabel("Current (A)")
plt.title("Plot of the inverse of the radius of the beam path as a function of " +\
          "current through the Helmholtz coil", fontsize = 8)
plt.legend()
plt.savefig("ChuKi I vs r-1.pdf")

#calculation of charge/mass by the method described in the lab report;
i0 = -Bext/k
idiff = idata-i0
cm = average((vdata/((radius*k)**2))*((1/(idiff)**2)))

#propagation of error
#k is constants multiplied, propagating no error
#treating i0 as constant
err_invert_r = exprop(-1, radius, err_r)
err_invert_r_sq = exprop(2, 1/radius, err_invert_r)
err_diff_i = adsubprop(0, err_i)
err_invert_diff_i = exprop(-1, idiff, err_diff_i)
err_invdiff_i_sq = exprop(2, 1/(idiff**2), err_invert_diff_i)
z = inv_r**2*1/k**2
err_a = muldivprop(inv_r**2, 1/k**2, z, err_invert_r_sq, 0)
err_b = muldivprop(z, vdata, z*vdata, err_a, err_v)
err_tot = muldivprop(z*vdata, 1/idiff**2, z*vdata*1/idiff**2, err_b, err_invdiff_i_sq)

#reduced chi-square calculation
fit = sta.chisquare(inv_r, p_opt[0]*idata+p_opt[1])[0]/2

#print-out
print("The final, calculated charge-to-mass ratio of an electron is")
print(cm)
print("The observed errors for each datapoint are:")
print(err_tot)
print("The resultant reduced chi-square statistic is")
print(fit)