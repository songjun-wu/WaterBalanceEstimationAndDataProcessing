import pandas as pd
import numpy as np
import math

def deltaE_18O(airtemp, rh, deltaL, deltaP):
    """
    Parameters:
        - airtemp: (array of) daily average air temperatures [Celsius].
        - rh: (array of) daily average relative humidity [%].
        - deltaL, P: (array of) isotopic composition in liquid phase, precipitation
    """
    alphaPlus = np.exp((-7.685 + 6.7123*(1e3/airtemp)-1.6664*(1e6/airtemp**2)+ 0.35041*(1e9/airtemp**3))/1000)
    epsilonPlus = (alphaPlus-1)*1000  # per mill
    epsilonK = (1-rh)*0.5*0.0286       # decimal
    deltaA = (deltaP - epsilonPlus)/alphaPlus
    deltaE = ((deltaL-epsilonPlus)/alphaPlus - rh*deltaA -epsilonK)/(1-rh-epsilonK/1000)

    return deltaE


def deltaE_2H(airtemp, rh, deltaL, deltaP):
    """
    Parameters:
        - airtemp: (array of) daily average air temperatures [Celsius].
        - rh: (array of) daily average relative humidity [%].
        - deltaL, P: (array of) isotopic composition in liquid phase, precipitation
    """
    alphaPlus = np.exp(
        (1158.8*(airtemp**3/1e9) - 1620.1*(airtemp**2/1e6) + 794.84*(airtemp/1e3) - 161.04 + 2.9992*(1e9/airtemp**3))/1000)
    epsilonPlus = (alphaPlus - 1) * 1000  # per mill
    epsilonK = (1 - rh) * 0.5 * 0.025  # decimal
    deltaA = (deltaP - epsilonPlus) / alphaPlus
    deltaE = ((deltaL - epsilonPlus) / alphaPlus - rh * deltaA - epsilonK) / (1 - rh - epsilonK / 1000)

    print(alphaPlus)
    print(epsilonPlus)
    print(epsilonK)
    print(deltaA)
    print(deltaE)

    return deltaE

deltaE_18O(276.7733,0.8033811,-0.0073185*1000, -0.0072222*1000)
deltaE_2H(276.7733,0.8033811,-0.054414625*1000, -0.0510233*1000)