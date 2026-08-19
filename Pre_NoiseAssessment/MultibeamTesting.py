# -*- coding: utf-8 -*-
"""
Created on Wed Aug  6 17:32:09 2025

@author: nikol
"""

import numpy as np
from scipy.fft import rfft,rfftfreq,irfft
import yaml
import matplotlib.pyplot as plt
from scipy.io import loadmat
from scipy.interpolate import interp1d
from scipy.signal.windows import tukey

dtheta = 107/180*np.pi
da     = 5/180*np.pi

wd = 20.0

#%%
#Simple Directivity Plot 1D
along_a = np.arange(-80.0,80.0,0.1) #Angles
along_a_rad = along_a/180*np.pi+1e-15 #Angles in radians

e = np.pi*along_a_rad/dtheta #Element Transducer Directivity Parameter
de = 20.0*np.log10(np.abs(np.sin(e)/e)) #Element Transducer Directivity

g = np.pi/da*(np.sin(along_a_rad)-np.sin(0)) #Beam Pattern Parameter
dp = 20*np.log10(np.abs(np.sin(g)/g)) #Beam Pattern
#Plotting
plt.figure()
plt.plot(along_a,de,'-k',label='Across Track = Element Transducer')
plt.plot(along_a,dp,'--b',label='Along Track Directivity')
plt.grid()
plt.legend()
plt.ylabel('Directivity Pattern [dB]')
plt.xlabel('Angle [deg]')
plt.tight_layout()


#%%

# Directivity Plot 2D
x = np.arange(-100.0,100.0,1.0) #X-Vector
y = np.arange(0,100.0,1.0) #Y-Vector
xgrid, ygrid = np.meshgrid(x,y) #Grid
dist = np.sqrt(ygrid**2+xgrid**2+(wd)**2) #Distance to point
dist_h = np.sqrt(ygrid**2+xgrid**2) #Distance

plt.figure()

# Directivity
psi = np.abs(np.arctan(dist_h/(wd))+1e-19) #Angle point
e = np.pi*psi/dtheta #Element Transducer Directivity Parameter
de_psi = 20.0*np.log10(np.abs(np.sin(e)/e)) #Element Transducer Directivity

plt.subplot(1,3,1)
plt.pcolormesh(x,y,de_psi,vmin=-30,vmax=0,shading='auto')


# g = np.pi/da*(np.sin(psi)-np.sin(0)) #Beam Pattern Parameter
# da_psi = 20*np.log10(np.abs(np.sin(g)/g)) #Beam Pattern



phi = np.abs(np.arcsin(xgrid/dist)+1e-19)
g = np.pi/da*(np.sin(phi)-np.sin(0)) #Beam Pattern Parameter
da_phi = 20*np.log10(np.abs(np.sin(g)/g)) #Beam Pattern
#Plotting

plt.subplot(1,3,2)
plt.pcolormesh(x,y,da_phi,vmin=-30,vmax=0,shading='auto')

plt.subplot(1,3,3)
plt.pcolormesh(x,y,de_psi+da_phi,vmin=-30,vmax=0,shading='auto')
plt.colorbar()

#%%
# Underwater Noise Measures
#Spatial Coordinates
x = np.arange(-profile_len/2,profile_len/2,profile_ds)
y = np.arange(hydro_seldistmin,hydro_seldistmax,hydro_seldistste)
xgrid, ygrid = np.meshgrid(x,y)
dist = np.sqrt(xgrid**2+ygrid**2+(depth-td)**2) #Distance to point
dist_h = np.sqrt(xgrid**2+ygrid**2) #Horizontal distance
#Directivity
theta = np.abs(np.arctan(dist_h/(depth-td))+1e-19) #Angle
e = np.pi*theta/dtheta #Element Transducer Directivity Parameter
dfy = 20.0*np.log10(np.abs(np.sin(e)/e)) #Element Transducer Directivity
beta = np.abs(np.arctan(xgrid/dist)+1e-19)
a = np.pi/phi*(np.sin(beta)-np.sin(0)) #Beam Pattern Parameter
dfx = 20*np.log10(np.abs(np.sin(a)/a)) #Beam Pattern 
#Corrections
ed = 10*np.log10(sigdur) #Event duration for dB calculation
fw,_ = WeightFunction(group,f) #Weight in dB for functional hearing group
#tl = GeometricalSpreading(geom,0.0,dist,dBred,depth,safetydb) #Transmission loss by geometric spreading
tl = GeometricalSpreading(yaml_input['spreading'],0.0,dist_h,\
    yaml_input['spreading_dBred'],\
    yaml_input['spreading_wd'],\
    yaml_input['spreading_safety'],\
    yaml_input['kraken_count'],\
    yaml_input['kraken_model'],\
    yaml_input['kraken_depth'],\
    yaml_input['kraken_selection'])  
absorp = -SimpleAbsorption(sea,f)*dist #Absorption
#Calculation
spl_peak = spl+tl+dfx+dfy+absorp
spl_rms  = spl-3.+tl+fw+dfx+dfy+absorp
sel_shot = spl-3.+tl+ed+fw+dfx+dfy+absorp
sel_shot = 10**(sel_shot/10)  
sel_cum = 10.*np.log10(np.sum(sel_shot,axis=1))
#Plotting
plt.figure(figsize=(10,8))
plt.subplot(3,1,1)
plt.pcolormesh(x,y,spl_peak,shading='auto',vmin=100,vmax=200)