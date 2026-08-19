# -*- coding: utf-8 -*-
"""
Created on Wed Mar  5 09:08:51 2025

@author: nikol
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
# #Input parameters
# gd = 0.8            #Gun Depth below water surface
# v = 1500.           #Water velocity [m/s]
# T_dom = 0.00227     #Dominant period, spacing of inflection points [s]
# z = 40            #Water Depth [m]
# d_max = 300.         #Maximum distance for calculation

# #Dipole separation = imaginary ghost source
# sep = 2*gd          #Dipole distance [m]

# #Calculation of wavenumber
# f_dom = 10000#1/T_dom     #Dominant frequency [Hz]
# lamb = v/f_dom      #Wavelength [m]
# k=2*np.pi*lamb      #Wavenumber [1/m]

# #Calculation of the emergence angle, measured from the vertical
# d = np.arange(0.0,d_max,1.0)    #Distance vector [m]
# tau = np.arctan(d/z)

# #Calculation of directivity [ ]
# def directivity(angle,wavenum,dipolsep):
#     pfunc_theta = np.sin(wavenum*dipolsep/2*np.cos(angle))
#     return pfunc_theta

# dir_f = directivity(tau,k,sep) 
# dir_f_dB = 20*np.log10(dir_f)

# #Polar Plot
# theta = np.arange(0.,2*np.pi,0.1)
# dir_t = directivity(theta,k,sep) 
# dir_t_dB = 20*np.log10(dir_t)

# #Plot
# fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
# ax.plot(theta, dir_t_dB)
# #ax.set_rmax(0)
# #ax.set_rticks([0.5, 1, 1.5, 2])  # Less radial ticks
# ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
# ax.grid(True)
# #ax.set_title("A line plot on a polar axis", va='bottom')
# #plt.show()


# fig,ax = plt.subplots(1,1)
# ax.plot(d,dir_f_dB)
# ax.set_xlabel('Distance [m]')
# ax.set_ylabel('Dipole directivity [dB]')
# ax.grid()

#%%
Ios = 1.0
zr = np.arange(0.0,20.0,0.01)
zs = 0.7
f_dom =440.0     #Dominant frequency [Hz]
v = 1500.           #Water velocity [m/s]
lamb = v/f_dom      #Wavelength [m]
k=2*np.pi/lamb      #Wavenumber [1/m]
mu = -1 

x = np.arange(-25.0,25.0,0.01)
xv, zrv = np.meshgrid(x,zr)
r = np.sqrt(xv**2+zrv**2)+1e-13


#I = Ios*(1/r**2)*(1+mu**2+2*mu*(1-2*np.sin(k*zs*zrv/r)**2))
#I = (1+mu**2+2*mu*(1-2*np.sin(k*zrv*zs/r**2)**2))
#I = Ios*(1/r**2)*(4*np.sin(k*zrv*zs/r)**2)
I=np.sin(k*zrv*zs/r)**2
I=20.0*np.log10(I)

fix,ax = plt.subplots(1,1,)
fix.set_size_inches(7, 3)
fix.set_dpi(300)
pcm=ax.pcolormesh(xv,zrv,I,vmin=-60,vmax=0)
cb=plt.colorbar(pcm)
cb.set_label('Relative Intensity [dB]') 
ax.invert_yaxis()
ax.set_xlabel('Horizontal Distance from Source [m]')
ax.set_ylabel('Water Depth [m]')
plt.tight_layout()
fix.savefig('test2png.png', dpi=300)

zr = 20.
x = np.arange(0,2000.0,0.1)
r=np.sqrt(x**2+zr**2)
I=np.sin(k*zr*zs/r)**2
I=20.0*np.log10(I)

fix,ax = plt.subplots(1,1)
fix.set_size_inches(7, 3)
fix.set_dpi(300)
ax.plot(x,I)
ax.set_ylabel('Relative Intensity [dB]') 
ax.set_xlabel('Horizontal Distance from Source [m]')
ax.grid()
plt.tight_layout()

fix.savefig('test3png.png', dpi=300)



#%%
x=np.arange(0,10,0.01)
y=x**(1/3)

plt.figure(),plt.plot(x,y)

#%%
f=440 #frequency [Hz]
v=1500 #water vel [m/s]
w = 10 #Wind speed [m/s] Beaufort 5: 10
tau = np.arange(0.,90.,0.1)
tau_rad = tau/180*np.pi
lamb = v/f     #Wavelength [m]

k=2*np.pi/lamb      #Wavenumber [1/m]
hlist = np.array([0.0,.1,.2,.6,1,2,3,4,5.5,7,9])

fig,ax = plt.subplots(1,1) 
fig.set_size_inches(7, 4)
fig.set_dpi(300)
n_lines = np.size(hlist)
cmap = mpl.colormaps['cividis']
colors = cmap(np.linspace(0, 1, n_lines))


for ind,h in enumerate(hlist): 
    theta = 2*k*h*np.sin(tau_rad)
    r=-np.exp(-0.5*theta**2)
    labelstr = 'Beaufort '+str(ind)+': h='+str(h) +' m'
    ax.plot(tau,r,label=labelstr,color=colors[ind])
ax.grid()
plt.legend(loc='upper right', bbox_to_anchor=(1.2, 1.2))
ax.set_xlabel('Incidence Angle [°]')
ax.set_ylabel('Sea Surface \n Reflection Coefficient [ ]')
plt.tight_layout()
fig.savefig('test3png.png', dpi=300)

#rl = 8.6*10**-9*f**2*w**4*np.sin(tau)**2
#rl = 3.0*10**-4*f**2*h**2*tau**2