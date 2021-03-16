# -*- coding: utf-8 -*-
"""
Created on Sun Mar 14 11:09:57 2021

@author: Dipankar
References
----------
Oh et al. (1992): An empirical model and an inversion technique for radar scattering from bare soil surfaces. IEEE TGRS 30(2). 370-381.

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

def Fresnel0(e):
    return np.abs( (1.-np.sqrt(e))/(1.+np.sqrt(e))   )**2.

c0=299792458.  # speed of light [m/s]

def f2lam(f):
    """
    given the frequency in GHz,
    return the wavelength [m]
    """
    return c0/(f*1.E9)

  
def Reflectivity(eps, theta):
    co = np.cos(theta)
    si2 = np.sin(theta)**2.
    rho_v = (eps*co-np.sqrt(eps-si2))/(eps*co+np.sqrt(eps-si2))
    rho_h = (co-np.sqrt(eps-si2))/(co+np.sqrt(eps-si2))
    v = np.abs(rho_v)**2.
    h = np.abs(rho_h)**2.
    return v,h


##-----------------------------------------------------------------------------        
      
eps = 11.3-1.5j   # 
G0 = Fresnel0(eps)  # nadir fresnel reflectivity
theta = np.deg2rad(40)
freq = 3.0  # GH
s = 0.01  # m
ks = (2.*np.pi/f2lam(freq))*s

a = 1./(3.*G0)
p = (1. - (2.*theta/np.pi)**a * np.exp(-ks))**2.
q = 0.23*(G0)**0.5 * (1.-np.exp(-ks))

G = Reflectivity(eps, theta)
vva = 0.7*(1.-np.exp(-0.65*ks**1.8))
vvb = np.cos(theta)**3. * (G[0]+G[1]) / np.sqrt(p)
vv0 = vva*vvb
hh0 = p * vv0
hv0 = q * vv0

## Backscatter in dB scale
vv0dB = 10*np.log10(vv0)
hh0dB = 10*np.log10(hh0)
hv0dB = 10*np.log10(hv0)     