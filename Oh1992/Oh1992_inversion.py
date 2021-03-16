# -*- coding: utf-8 -*-
"""
Created on Sun Mar 14 12:31:27 2021

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

## Description: Given sigma_0_vv, sigma_0_hh, and sigma_0_hv, the inverse
## model computes eps, s, and mv (assuming loamy soil)

##Input Variables:
#    %sigma_0_vv (dB), sigma_0_hh (dB), and sigma_0_hv (dB)
#    %theta: Incidence angle (degrees)
#    %f: frequency (GHz)
    
#%Output Product:
#    %eps: absolute value of complex dielectric constant of soil medium
#    %s: rms height (m)
#    %mv: volumetric moisture content (g/cm3)
    
#%Book Reference: Section 10-5


sigma0vvdB = -14.1
sigma0hhdB = -16.0
sigma0hvdB = -26.5
theta = 35.  ##Incidence angle
f = 5.0 ##GHz


k = 2*np.pi*f/0.3; #calculate the wave number 

theta_rad = theta*np.pi/180; #represent angle in radians

sigma_0_vv = np.power(10,(sigma0vvdB/10)) #%represent data in linear scale
sigma_0_hh = np.power(10,(sigma0hhdB/10))
sigma_0_hv = np.power(10,(sigma0hvdB/10))


p = sigma_0_hh / sigma_0_vv; #calculate the p-ratio
q = sigma_0_hv / sigma_0_vv; #calculate the q-ratio

gamma0 = np.arange(0.001,1.0, 0.001) # set Gamma0 range of values (fine increments)

a1 = 2.0*theta_rad / np.pi;

## calculate the error function for different values of Gamma0 and find the minimum

err = 1- a1**(1./(3.*gamma0)) *(1- q/(0.23*np.sqrt(gamma0))) - np.sqrt(p)

abs_err = np.abs(err);
min_err = np.min(abs_err); #find the value of minimum error

gamma0_min = gamma0[np.where(abs_err == min_err)]


            
eps = ( (1 + np.sqrt(gamma0_min)) /(1- np.sqrt(gamma0_min)) )**2

s = -1/k * np.log(1 - q /(0.23*np.sqrt(gamma0_min)))  ## meter scale
s0 = s*100 ##cm


#-----------------------------------------------------------
#------ estimate mv based on dielectric model of soil (loamy)
## Dobson et al. 1985 dielectic mixture model
## Eq ref. Ulaby et al (2014)
t= 30.; # temperature in degree C
S = 30.6 / 100; # fraction of Sand
C= 13.5 /100; # fraction of Clay

rho_b = 1.36; # soil bulk density g/cm3
f_hz = f * 1.0e9; # transform from GHz to Hz

beta1 = 1.27 - 0.519 * S - 0.152* C;  #eq: 4.68b

alpha = 0.65; # eq: 4.68a

#Dielectric Constant of Pure Water
    
ew_inf = 4.9; # eq: E.15
ew_0 = 88.045 - 0.4147 * t + 6.295e-4 * t**2 + 1.075e-5 * t**3;   
tau_w = (1.1109e-10 - 3.824e-12*t +6.938e-14*t*t - 5.096e-16*t**3)/2/np.pi

epsrW = ew_inf +(ew_0-ew_inf)/(1 + (2*np.pi*f_hz*tau_w)**2)


##--volumetric soil moisture estimation---------------
  
mvv = np.arange(0.001,0.5, 0.01)  
#-calculate error function
err2 = abs(eps**(alpha) - 1- 0.66*rho_b - mvv**beta1 * epsrW**alpha + mvv);

min_err2 = np.min(err2) #find minimum error
mv = mvv[np.where(err2 == min_err2)] #find mv corresponding to minimum error

print('Estimated rms height s (cm): ', s0)
print('Estimated volumetric soil moisture: ', mv)


