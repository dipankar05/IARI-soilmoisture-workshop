# -*- coding: utf-8 -*-
"""
Created on Sun Mar 14 18:53:14 2021

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


## Description: Given sigma_0_vv, sigma_0_hh, and sigma_0_hv, the inverse
## model computes s, and mv 

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

mv0 = np.arange(0.05,0.5,0.01) # set Gamma0 range of values (fine increments)



## First estimates s1 and mv1
ks = ((-1)*3.125*np.log(1 - sigma_0_hv/(0.11 * mv0**0.7 * (np.cos(theta_rad))**2.2)))**0.556
err = (1 - (2.*theta_rad/np.pi)**(0.35*mv0**(-0.65)) * np.exp(-0.4 * ks**1.4))-p
abs_err = np.abs(err);
min_err = np.min(abs_err); #find the value of minimum error
mv1 = mv0[np.where(abs_err == min_err)]
ks1 = ((-1)*3.125*np.log(1 - sigma_0_hv/(0.11 * mv1**0.7 * (np.cos(theta_rad))**2.2)))**0.556
s1 = ks1/k


## Second estimate s2 and mv2
ks2 = (np.log(1-(q/(0.095 * (0.13 + np.sin(1.5*theta_rad))**1.4))) /(-1.3))**(10./9.)
s2 = ks2/k

xx = (1-p)/np.exp(-0.4 * ks2**1.4)
if xx<=0:
    mv2 =0
else:
    yy = np.log(xx)/(0.35*np.log(2*theta_rad/np.pi))
    mv2 = yy**(-100/65)

         
          
          
## Third estimate mv3
mv3 = ((sigma_0_hv/(1 - np.exp(-0.32 * ks2**1.8)))/(0.11 * np.cos(theta_rad)**2.2))**(1/0.7)

## weighted average s and mv-------------------------------------
sf = (s1 + 0.25*s2)/(1+0.25)
mvf = (mv1+mv2+mv3)/3

print('Estimated rms height s (cm): ', sf*100)
print('Estimated volumetric soil moisture: ', mvf)