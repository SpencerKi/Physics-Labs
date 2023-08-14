 # -*- coding: utf-8 -*-
"""
Code for Kater's Pendulum Lab
"""
import numpy as np
import scipy.optimize as opt
import scipy.stats as sta
import matplotlib.pyplot as plt

#mass position and period are related via a quadratic equation.
def quadratic(x, c, b):
    return c*x**2 + b

#coarse adjustment mass positioning data.
coarse = np.loadtxt("coarse.txt", delimiter = ",", skiprows = 1, unpack = True)
coarsed = coarse[0]
coarsetop = coarse[1]
coarsebottom = coarse[2]

#curve-fitting the upright and inverted plots.
coarsetopopt, coarsetopcov = opt.curve_fit(quadratic, coarsed, coarsetop)
coarsebottomopt, coarsebottomcov = opt.curve_fit(quadratic, coarsed, coarsebottom)

#function of upright pendulum relationship.
def coarsetopquad(x):
    return coarsetopopt[0]*x**2 + coarsetopopt[1]

#function of inverted pendulum relationship.
def coarsebottomquad(x):
    return coarsebottomopt[0]*x**2 + coarsebottomopt[1]

#for use with fsolve().
def coarsediffquad(x):
    return (coarsetopopt[0] - coarsebottomopt[0])*x**2 + coarsetopopt[1] - coarsebottomopt[1]

#calculate intersection.
coarseloc = opt.fsolve(coarsediffquad, 26)
coarseperiod = coarsetopquad(coarseloc)

plt.figure(0)
plt.scatter(coarsed, coarsetop, label="Upright")
plt.scatter(coarsed, coarsebottom, label="Inverted")
plt.plot(coarsed, coarsetopquad(coarsed))
plt.plot(coarsed, coarsebottomquad(coarsed))
plt.plot(coarseloc, coarseperiod, 'ro', label="Intersection")
plt.xlabel("Mass Position")
plt.ylabel("Period of 8 Oscillations")
plt.title("Coarse Adjustment Mass, Upright & Inverted")
plt.legend()
plt.savefig("coarse.png")

#fine adjustment mass positioning data.
fine1 = np.loadtxt("fine1.txt", delimiter = ",", skiprows = 1, unpack = True)
fine1d = fine1[0]
fine1top = fine1[1]
fine1bottom = fine1[2]

#even finer positioning data.
fine2 = np.loadtxt("fine2.txt", delimiter = ",", skiprows = 1, unpack = True)
fine2d = fine2[0]
fine2top = fine2[1]
fine2bottom = fine2[2]

#combining the datasets.
finecombd = np.append(fine1d, fine2d)
finecombtop = np.append(fine1top, fine2top)
finecombbottom = np.append(fine1bottom, fine2bottom)

#curve-fitting the upright and inverted plots.
finecombtopopt, finecombtopcov = opt.curve_fit(quadratic, finecombd, finecombtop)
finecombbottomopt, finecombbottomcov = opt.curve_fit(quadratic, finecombd, finecombbottom)

#function of upright pendulum relationship.
def finecombtopquad(x):
    return finecombtopopt[0]*x**2 + finecombtopopt[1]

#function of inverted pendulum relationship.
def finecombbottomquad(x):
    return finecombbottomopt[0]*x**2 + finecombbottomopt[1]

#for use with fsolve().
def finecombdiffquad(x):
    return (finecombtopopt[0] - finecombbottomopt[0])*x**2 + finecombtopopt[1] - finecombbottomopt[1]

#calculate intersection.
finecombloc = opt.fsolve(finecombdiffquad, 5)
finecombperiod = finecombtopquad(finecombloc)

plt.figure(1)
plt.scatter(finecombd, finecombtop, label="Upright")
plt.scatter(finecombd, finecombbottom, label="Inverted")
plt.plot(finecombd, finecombtopquad(finecombd))
plt.plot(finecombd, finecombbottomquad(finecombd))
plt.plot(finecombloc, finecombperiod, 'ro', label="Intersection")
plt.xlabel("Mass Position")
plt.ylabel("Period of 8 Oscillations")
plt.title("Fine Adjustment Mass, Upright & Inverted")
plt.legend()
plt.savefig("fine.png")

#measured pendulum length.
L = 1.002
#calculated period.
T = finecombperiod / 8
#calculated period from course data.
cT = coarseperiod / 8

#relationship between pendulum length and period.
def grav(L, T):
    return ((2*np.pi)**2)*(L/T**2)

#calculate acceleration due to gravity.
g = grav(L, T)
#calculate acceleration due to gravity with coarse data.
cg = grav(L, cT)
#measurement error.
merr = (np.sqrt(g/((2*np.pi)**2))*np.sqrt((5e-7/finecombperiod)**2 + \
                (5e-5/L)**2))*(2*g/((2*np.pi)**2))*(2*np.pi)**2/((L/T)**2)
#measurement error with coarse data.
cmerr = (np.sqrt(g/((2*np.pi)**2))*np.sqrt((5e-7/finecombperiod)**2 + \
                (5e-5/L)**2))*(2*g/((2*np.pi)**2))*(2*np.pi)**2/((L/cT)**2)