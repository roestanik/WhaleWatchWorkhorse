# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 09:04:19 2023

@author: nikol
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.fft import rfft,rfftfreq,irfft


def ReadTXT(name,headerlines,colums2read):
    """
    Function to read in txt files

    Parameters
    ----------
    name : STRING
        Full path of the file to be read in.
    headerlines : INT
        Number of lines to skip.
    colums2read : List of INT
        Columns in the file to be read in.

    Returns
    -------
    header : List of strings
        Header words of the columns read in.
    data_app : FLOAT Numpy array
        Numpy array of the data read in.

    """  
    #Open the file
    file = open(name)
    #Initialisation
    counter  = 0  # Counter for lines
    data_app = [] # Data variable
    header   = [] # Header variable
    
    #Main loop
    for line in file: #Skip trough all lines
        counter +=1 #Counter goes up every cycle
        line=line.strip() #Read the next line
        columns = line.split('\t') #Extract the columns from the line seperated by tabs
        if counter ==headerlines: #Get the headers of the columns of interest, headers are in the last line before the data
            for step,c in enumerate(colums2read): #Read the columns of interest
                header.append(columns[int(c)])   #save for output
        if counter>headerlines: #Read all lines containing data
            data = np.zeros(np.shape(colums2read))*np.nan #initialisation of the data in the line
            try: #Error handling: Do not stop if one line is broken/Too short
                for step,c in enumerate(colums2read):
                    #print(step,c)
                    #get depth if there are enough columns
                    if len(columns[int(c)])>0:
                        data[step]=float(columns[int(c)])    
            except: #Error handling:
                print('error reading line (not all data points given?)', counter, 'skipping further entries in this line')
            data_app.append(data) #Add data to output variable
    data_app=np.array(data_app) #Convert output to numpy array
    return header, data_app

header, data_app= ReadTXT('Profile.txt',1,[0,1])

distance = data_app[:,0]
distance = np.max(distance)-distance
depth    = data_app[:,1]

sortind = np.argsort(distance)
distance = distance[sortind]
depth = depth[sortind]
depth = -1.0*depth

distante_int = np.arange(0,np.floor(np.max(distance)/1000)*1000,900)
depth_int = interp1d(distance, depth)
depth_int = depth_int(distante_int)

plt.figure()
plt.plot(distance,depth)
plt.plot(distante_int,depth_int)
plt.gca().invert_yaxis()
plt.xlabel('Range [m]')
plt.ylabel('Depth [m]')
plt.grid()

out = open('actup_profile.bty','w+')
out.write(str(np.size(depth_int))+'\n')
for depthind, depthstep in enumerate(depth_int):
    out.write(str(np.round(distante_int[depthind]/1000,1))+'\t'+str(np.round(depthstep,1))+'\n')
    
out.close()

#%%
import seawater as sw

#acoustic properties
def coppens(temp,sal,depth):
    t=temp/10
    D=depth/1000
    
    c_0 = 1449.05 + 45.7*t - 5.21*t**2 + 0.23*t**3 + (1.333 - 0.126*t + 0.009*t**2)*(sal - 35)
    c_d = c_0 + (16.23 + 0.253*t)*D + (0.213-0.1*t)*D**2 + (0.016 + 0.0002*(sal-35))*(sal - 35)*t*D
    
    return c_d

lat = 36.34
layer1_top_sal=38.8
layer1_top_t = 32. 
layer1_top_d = 0   
layer1_top_v = sw.eos80.svel(layer1_top_sal, layer1_top_t, sw.eos80.pres(layer1_top_d, lat))
layer1_top_rho = sw.eos80.dens(layer1_top_sal, layer1_top_t, sw.eos80.pres(layer1_top_d, lat))
print('Layer 1:')
print('Top: ',layer1_top_v ,' m/s', layer1_top_rho, ' kg/m^3')

layer1_bottom_sal = 38.8
layer1_bottom_t   = 27. 
layer1_bottom_d    = 50   
layer1_bottom_v =sw.eos80.svel(layer1_bottom_sal, layer1_bottom_t, sw.eos80.pres(layer1_bottom_d, lat))
layer1_bottom_rho  = sw.eos80.dens(layer1_bottom_sal, layer1_top_t, sw.eos80.pres(layer1_bottom_d, lat))
print('Btm: ', layer1_bottom_v,' m/s',layer1_bottom_rho , ' kg/m^3')

layer2_top_sal=40.5
layer2_top_t = 23. 
layer2_top_d = 100.   
layer2_top_v = sw.eos80.svel(layer2_top_sal, layer2_top_t, sw.eos80.pres(layer2_top_d, lat))
layer2_top_rho = sw.eos80.dens(layer2_top_sal, layer2_top_t, sw.eos80.pres(layer2_top_d, lat))
print('Layer 2:')
print('Top: ', layer2_top_v,' m/s', layer2_top_rho, ' kg/m^3')

layer2_bottom_sal = 41
layer2_bottom_t   = 22. 
layer2_bottom_d    = 1800.       
layer2_bottom_v = sw.eos80.svel(layer2_bottom_sal, layer2_bottom_t, sw.eos80.pres(layer2_bottom_d, lat))
layer2_bottom_rho = sw.eos80.dens(layer2_bottom_sal, layer2_top_t, sw.eos80.pres(layer2_bottom_d, lat))
print('Btm: ',layer2_bottom_v ,' m/s',layer2_bottom_rho, ' kg/m^3')

fig,[ax1,ax2,ax3,ax4]=plt.subplots(1,4,sharey=True)
ax1.plot([layer1_top_t,layer1_bottom_t,layer2_top_t,layer2_bottom_t],\
         [layer1_top_d,layer1_bottom_d,layer2_top_d,layer2_bottom_d],'or')
ax1.plot([layer1_top_t,layer1_bottom_t,layer2_top_t,layer2_bottom_t],\
         [layer1_top_d,layer1_bottom_d,layer2_top_d,layer2_bottom_d],'-k')
ax1.grid()
ax1.set_xlabel('Temperature [deg]')    
ax1.set_ylabel('Depth [m]')

ax2.plot([layer1_top_sal,layer1_bottom_sal,layer2_top_sal,layer2_bottom_sal],\
         [layer1_top_d,layer1_bottom_d,layer2_top_d,layer2_bottom_d],'or')
ax2.plot([layer1_top_sal,layer1_bottom_sal,layer2_top_sal,layer2_bottom_sal],\
         [layer1_top_d,layer1_bottom_d,layer2_top_d,layer2_bottom_d],'-k')
ax2.grid()
ax2.set_xlabel('Salinity [psu]')    

ax3.plot([layer1_top_v,layer1_bottom_v,layer2_top_v,layer2_bottom_v],\
         [layer1_top_d,layer1_bottom_d,layer2_top_d,layer2_bottom_d],'or')
ax3.plot([layer1_top_v,layer1_bottom_v,layer2_top_v,layer2_bottom_v],\
         [layer1_top_d,layer1_bottom_d,layer2_top_d,layer2_bottom_d],'-k')
ax3.grid()
ax3.set_xlabel('PWave velocity [m/s]')    

ax4.plot([layer1_top_rho,layer1_bottom_rho,layer2_top_rho,layer2_bottom_rho],\
         [layer1_top_d,layer1_bottom_d,layer2_top_d,layer2_bottom_d],'or')
ax4.plot([layer1_top_rho,layer1_bottom_rho,layer2_top_rho,layer2_bottom_rho],\
         [layer1_top_d,layer1_bottom_d,layer2_top_d,layer2_bottom_d],'-k')
ax4.grid()
ax4.set_xlabel('Density [kg/m$^3$]')  

ax4.invert_yaxis()
plt.tight_layout()

#%%
from scipy.signal import welch,spectrogram
from matplotlib.mlab import psd
from scipy.signal.windows import tukey

t = np.array([0,0.0151234068764529,0.0551572645091517,0.0867920836013358,0.105204466403464,0.115247735808032,0.125289338575051,0.148764188699555,0.173924861293962,0.230984496449729,0.284704152779228,0.356892857587711,0.439158164551863,0.506312734876401,0.553314100889606,0.588558042189903,0.61539620406644,0.637190231633698,0.66571088270503,0.685864087016369,0.701006489769524,0.71110975812607,0.724574671597095,0.739712074437585,0.751476165700714,0.773285193005964,0.791720908733859,0.803446667333228,0.82021322696676,0.830276496021991,0.848788877077399,0.857224655977143,0.865657101601772,0.870751235037695,0.877521191118187,0.887642792487835,0.899438549864503,0.91456761951722,0.943081604038333,0.973256411466028,1.01854528880845,1.06719581126535,1.13768702714105,1.20651075355995,1.27366865715959,1.33917073776532,1.37948381293821,1.42483602250772,1.48191232403903,1.53227366852916,1.59103079261771,1.6497979165316,1.69346763690631,1.74386564742264,1.78252123296405,1.8228359747745,1.87152649653272,1.91349206142442,1.95713344896069,2.01757139560649,2.07465103041291,2.1518905360809,2.20562519214839,2.25768902548397,2.30639121370506,2.38362571946038,2.43063708529892,2.50451494585241,2.58343360744008,2.66403642486009,2.71777608084025,2.83699948441303,2.93103388255298,3.02674243669994,3.1056577650125,3.31052145376435,3.96875224091864,4.92755027105908])
p = np.array([0,4.78553137167026,-4.09316603988693,-0.467843114736398,-0.359287510773191,-0.300357290050557,-0.238325873497928,-0.1762844570092,-0.126646573848637,-0.067681353349923,-0.042831786965706,-0.014867274839292,0.0100035414089343,0.031761911899121,0.0597076741454163,0.0938470781078142,0.146587407090315,0.211728769417093,0.258267956731551,0.252080564975579,0.196270289963058,0.143557460804426,0.0846447399696966,0.0381380524471573,0.0195396274109436,0.0567702272673709,0.121909089610123,0.174638168664715,0.221168606035154,0.242884476797456,0.165368330958598,0.0909458809977934,0.0227258226969771,-0.0827110855481963,-0.181944352141302,-0.268770335430137,-0.346291481237181,-0.377292189609479,-0.318348218974911,-0.225289844217652,-0.135321415362749,-0.0515528781838084,0.0105235380808648,0.0508903335432298,0.0664463123734294,0.0385842995910939,0.0138047327587669,-0.0140722799274484,0.0138809822707864,0.0387280486711128,0.0635813650312986,0.0698275064113858,0.0481516353932134,0.00477239353259229,-0.0603239690821402,-0.0882047317443253,-0.0788648944863191,-0.054024078045988,-0.0229796199535888,-0.00122624943145189,0.0205246211066599,0.0236833165688237,0.0206221204826891,0.00205369525450561,-0.0103148382975808,0.00214744465442962,0.0114860319206258,0.0208446190587521,0.0240045645128606,0.017862172468663,0.0054973888926817,-0.000616253335566697,-0.00364744961367203,0.00262619159055566,0.0119885287045167,0,0,0])

t=t/10
p_Pa = p*10**5
dt = 0.0005

t_int = np.arange(0,0.5,dt)
p_int = interp1d(t,p_Pa,kind='cubic',bounds_error=False,fill_value=0)
p_int = p_int(t_int)
window = tukey(p_int.size)
window[0:int(window.size/2)]=1
p_int = window*p_int
p_int = p_int/np.max(np.abs(p_int))*0.48*10**6

#Fourier transform
f        = rfftfreq(t_int.size*2, dt)
p0 = 10**(-6)
FFT      = rfft(p_int/p0,t_int.size*2)
FFT_I    = 20*np.log10(np.abs(FFT)/np.size(FFT))

plt.figure()
plt.subplot(1,2,1)
plt.plot(t,p_Pa,label='Gundalf Export')
plt.plot(t_int,p_int,'--r',label='Interpolation')
plt.grid()
plt.legend()
plt.xlabel('Time [s]')
plt.ylabel('Pressure [Pa]')
plt.subplot(1,2,2)
plt.plot(f,FFT_I,'-r')
plt.grid()
plt.xlabel('Frequency [Hz]')
plt.ylabel('Pressure [dB re 1 $\mu$ Pa/Hz]')

print(20*np.log10(np.max(np.abs(p_int))*10**6))