# -*- coding: utf-8 -*-
"""
Code for Thermoelectricity Lab
"""
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt

"""
Passive Heat Flow Experiment
"""
phf = np.loadtxt("passiveheatflow.txt", delimiter = ",", skiprows = 1, unpack = True)
phfi = phf[0]
phfierr = phf[1]
phftin = phf[2]
phftinerr = phf[3]
phftout = phf[4]
phftouterr = phf[5]
t0 = 21.4

phfpin = phfi * 5
phfpinerr = phfierr * 5
phftinphftoutdiff = phftin - phftout
phftinphftoutdifferr = np.sqrt(phftinerr**2 + phftout**2)
phftoutt0diff = phftout - t0

plt.figure(0)
plt.scatter(phfpin, phftinphftoutdiff, label="Input-Output Temperature Difference")
plt.scatter(phfpin, phftoutt0diff, label="Output-Ambient Temperature Difference")
plt.xlabel("Heat Flow into TEC (W)")
plt.ylabel("Temperature Difference (K)")
plt.title("Temperature Difference vs TEC Heat Flow")
plt.legend()

"""
Cooling Experiment A
"""
cooling1 = np.loadtxt("cooling1day1.csv", delimiter = ",", skiprows = 1, unpack = True)

cl1time = cooling1[0]
cl1i = cooling1[1]
cl1v = cooling1[2]
cl1tin = cooling1[4]
cl1tout = cooling1[3]
Pd = cl1i * cl1v

plt.figure(1)
plt.scatter(cl1time, cl1tin, label="Input Temperature")
plt.scatter(cl1time, cl1tout, label="Output Temperature")
plt.xlabel("Time (s)")
plt.ylabel("Temperature (C)")
plt.title("Temperature Change Over Time")
plt.legend()

plt.figure(2)
plt.scatter(cl1time, Pd)
plt.xlabel("Time (s)")
plt.ylabel("Power (W)")
plt.title("Electrical Power Input Over Time")

"""
Cooling Experiment B
"""
cooling2 = np.loadtxt("cooling2.csv", delimiter = ",", skiprows = 1, unpack = True)

cl2dci = cooling2[0]
cl2aci = cooling2[1]
cl2dcv = cooling2[2]
cl2sv = cooling2[3]
cl2tin = cooling2[5]
cl2tout = cooling2[4]

cl2pin = cl2aci * 5
cl2tempdiff = cl2tout - cl2tin
cl2pd = cl2dci * cl2dcv

plt.figure(3)
plt.scatter(cl2pin, cl2tempdiff)
plt.xlabel("Input Power (W)")
plt.ylabel("Input-Output Temperature Difference (K)")
plt.title("Temperature Difference vs. Input Power")

plt.figure(4)
plt.scatter(cl2pin, cl2pd)
plt.xlabel("Input Power (W)")
plt.ylabel("Total Electrical Power (W)")
plt.title("Total vs. Input Power")

plt.figure(5)
plt.scatter(cl2pin, cl2dcv)
plt.xlabel("Input Power (W)")
plt.ylabel("Total Voltage (V)")
plt.title("Total Voltage vs. Input Power")

"""
Calculation of Thermal Conductance of the System,
Seebeck Coefficient, and TEC Resistance
"""
def P(T, K):
    return K*T

kd_opt, kd_cov = opt.curve_fit(P, phftinphftoutdiff, phfpin)
sd_opt, sd_cov = opt.curve_fit(P, cl2tempdiff, cl2sv)
rd_opt, rd_cov = opt.curve_fit(P, cl2dci, cl2dcv)

Kd = kd_opt
Sd = sd_opt
Rd = rd_opt

Kderr = kd_cov
Sderr = sd_cov
Rderr = rd_cov

"""
Calculation of Hypothetical Lowest Achievable Temperature of the Thermal Reservoir
"""
Id = np.mean(cl2dci)
Iderr = 0.0005
def Tin(tout, i):
    return -(tout*(Sd*i+Kd)+1.5*Rd*i**2)/(2*Sd*i+Kd)

min_temp = Tin(t0, Id)

"""
Calculation of Optimum Drive Current
"""
def Pout(tin, i):
    return (tin)*Sd*i - 0.5*Rd*i**2-Kd*(35-tin) + Sd*(35-tin)*i + i**2*Rd

i_opt, i_cov = opt.curve_fit(Pout, 5, cl1i*cl1v)

OI = i_opt
OIerr = i_cov

"""
Calculation of Hypothetical Lowest Achievable Temperature Using Optimum Input Current
"""
new_min_temp = Tin(35, OI)