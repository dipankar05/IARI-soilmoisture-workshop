# -*- coding: utf-8 -*-
"""
Created on Sun Mar 14 16:57:18 2021

@author: Dipankar
References
----------
Oh (2004): Quantitative retrieval of soil moisture content and surface roughness from multipolarized radar observations of bare soil surface. IEEE TGRS 42(3). 596-601.
"""

  # ---------------------------------------------------------------------------------------
  # Copyright (C) 2021 by Microwave Remote Sensing Lab, IITBombay http://www.mrslab.in
 
  # This program is free software; you can redistribute it and/or modify it
  # under the terms of the GNU General Public License as published by the Free
  # Software Foundation; either version 3 of the License, or (at your option)
  # any later version.
  # This program is distributed in the hope that it will be useful, but WITHOUT
  # ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  # FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
  # more details.
 
  # You should have received a copy of the GNU General Public License along
  # with this program; if not, see http://www.gnu.org/licenses/
  # ---------------------------------------------------------------------------------------
  
import numpy as np
#import matplotlib.pyplot as plt


c0=299792458.  # speed of light [m/s]

def f2lam(f):
    """
    given the frequency in GHz,
    return the wavelength [m]
    """
    return c0/(f*1.E9)

  

##-----------------------------------------------------------------------------        
## Oh2004(mv, ks, theta):      
theta = np.deg2rad(35)
freq = 3.0 # GHz
s = 0.01  # m
ks = (2.*np.pi/f2lam(freq))*s
mv = 0.20

p = 1 - (2.*theta/np.pi)**(0.35*mv**(-0.65)) * np.exp(-0.4 * ks**1.4)
q = 0.095 * (0.13 + np.sin(1.5*theta))**1.4 * (1-np.exp(-1.3 * ks**0.9))

## VH calculation
a = 0.11 * mv**0.7 * np.cos(theta)**2.2
b = 1 - np.exp(-0.32 * ks**1.8)
hv0 = a*b

vv0 = hv0 / q
hh0 = hv0 / q * p

## Backscatter in dB scale
vv0dB = 10*np.log10(vv0)
hh0dB = 10*np.log10(hh0)
hv0dB = 10*np.log10(hv0) 