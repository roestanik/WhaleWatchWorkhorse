# -*- coding: utf-8 -*-
"""
ImpactAssessment_Functions.py

Functions to generate plots and data for the underwater noise impact 
assessment caused by seismics.

Author:             Nikolas Römer-Stange
Initial Draft:      2021.11.02
Last Update:        2020.11.04

Dependencies:       Defined in section 01 of the script.

Units:              SI unless stated differently.

Nomenclature:       Defined upon first call of variable

References:         See subrutines/functions
    
Comments:           Noise impact assessment as required by Danish
                    and German authorities.

"""

import numpy as np
from scipy.fft import rfft,rfftfreq,irfft
import yaml
import matplotlib.pyplot as plt
from scipy.io import loadmat

#%% Plotting defaults
#Open the wavelet 
yaml_file = open("PlottingDefaults.yaml",'r')
yaml_input = yaml.safe_load(yaml_file)
yaml_file.close()
    
fminplot   = yaml_input['fminplot']
fmaxplot   = yaml_input['fmaxplot']
fsteppllot = yaml_input['fsteppllot']

distmin    = yaml_input['distmin']
distmax    = yaml_input['distmax']
distste    = yaml_input['distste']

seldistmin = yaml_input['seldistmin']
seldistmax = yaml_input['seldistmax']
seldistste = yaml_input['seldistste']

hydro_seldistmin = yaml_input['hydro_seldistmin']
hydro_seldistmax = yaml_input['hydro_seldistmax']
hydro_seldistste = yaml_input['hydro_seldistste']

hydro_beamx   = yaml_input['hydro_beamx']
hydro_beamy   = yaml_input['hydro_beamy']
hydro_beamdx  = yaml_input['hydro_beamdx']


#%%
def GetSourceInfo(source_yaml_input,group,plotbool=True):
    #Open the wavelet 
    yaml_file = open("../Data/Sources/"+source_yaml_input+'_Measurements.yaml','r')
    yaml_input = yaml.safe_load(yaml_file)
    yaml_file.close()
    measurements  = yaml_input[source_yaml_input]
    wavelet_t = np.array(measurements['t'])
    wavelet_p = np.array(measurements['p'])
    wavelet_dt = np.unique(np.round(np.diff(wavelet_t),decimals=15))[0]
    
    #Fourier transform
    f        = rfftfreq(wavelet_t.size, wavelet_dt)
    FFT      = rfft(wavelet_p)
    FFT_A    = np.abs(FFT)
    FFT_I    = 20*np.log10(FFT_A)
    
    #Min+MaxFrequency
    f_temp = f[FFT_I>=np.max(FFT_I)-6]
    fmin=np.min(f_temp)
    fmax=np.max(f_temp)
    
    #Pulse Length
    wavesum = np.cumsum(wavelet_p**2)
    wavesum = wavesum/np.max(wavesum)*100
    pulse_length = np.max(wavelet_t[wavesum<97.5])-np.max(wavelet_t[wavesum<2.5])
    #SPL
    SPL=20.0*np.log10(np.max(np.abs(wavelet_p))/10.0**(-6.0))
    #SEL
    SEL = 10.0*np.log10(np.sum(wavelet_p**2.0)*wavelet_dt*10.0**12.0)
    
    if plotbool:
        #Also make Hearing Threshold plot
        f_whale = np.linspace(fminplot,np.log10(fmaxplot),fsteppllot)
        f_whale=10.0**f_whale
        plt.figure()
        plt.subplot(2,1,1)
        plt.plot(wavelet_t,wavelet_p)
        plt.grid()
        plt.xlabel('Time [s]')
        plt.ylabel('Pressure [Pa] \n at 1 m distance ')
        plt.axis('tight')
        plt.subplot(2,1,2)
        plt.plot(f,FFT_I,label=source_yaml_input+ ' (1 m Distance)')
        plt.plot(f_whale,threshold(group,f_whale),label='Hearing Threshold '+group)
        plt.xlabel('Frequency [Hz]')
        plt.ylabel('Sound Power [dB re 1$\mu Pa$]')
        plt.xscale('log')
        plt.axis('tight')
        plt.xlim(left=10)
        plt.ylim([0,220])
        plt.grid()
        plt.legend()
        plt.tight_layout()
        plt.savefig('Seismics_Source',dpi=300)

    return wavelet_t,wavelet_p,f,FFT_I,fmin,fmax,pulse_length,SPL,SEL

def threshold(group,f,estimated=True):
    f = f/1000
    #Open the wavelet 
    yaml_file = open("../Data/FilterFunctions/Southall-et-al-2009_Audiograms.yaml",'r')
    yaml_input = yaml.safe_load(yaml_file)
    yaml_file.close()
    if estimated:
        yaml_input = yaml_input['EstimatedAudiogram']
    else:
        yaml_input = yaml_input['NormalizedAudiogram']
    params  = yaml_input[group]
    
    T = params['T0']+params['A']*np.log10(\
        1+params['F1']/f)+\
        (f/params['F2'])**params['B']
        
    return T
    
def tts_onset(group,f):
    f = f/1000
    #Open the wavelet 
    yaml_file = open("../Data/FilterFunctions/Southall-et-al-2009_Weights_TTS.yaml",'r')
    yaml_input = yaml.safe_load(yaml_file)
    yaml_file.close()
    params  = yaml_input[group]
    
    E = params['K']-10.0*np.log10(\
        (f/params['f1'])**(2*params['a'])/\
        (1+(f/params['f1'])**params['a']*(f/params['f2'])**2)**params['B'])
        
    return E

def WeightFunction(group,f):
    """
    Weight function for generic band pass filter for marine mammal hearing 
    sensitivities

    Parameters
    ----------
    C : float
        Fitting parameter vor vertical position.
    f : float
        Frequency of interest [Hz].
    f1 : float
        LF transition value [Hz].
    f2 : float
        HF transition value [Hz].
    a : float
        LF decline of function [ ].
    b : float
        HF decline of function [ ].

    Returns
    -------
    w : float
        Weighting function amplitude [dB re 1 mu Pa].

    """
    #Open the wavelet 
    yaml_file = open("../Data/FilterFunctions/Southall-et-al-2009_Weights_TTS.yaml",'r')
    yaml_input = yaml.safe_load(yaml_file)
    yaml_file.close()
    params  = yaml_input[group]
    
    f=f/1000.0
    
    num = (f/params['f1'])**(2*params['a'])
    den1= (1+(f/params['f1'])**2)**params['a']
    den2=(1+(f/params['f2'])**2)**params['B']
    frac = num/(den1*den2)
    w=params['C']+10*np.log10(frac)
    
    f=10.0**(w/20.0)
    
    return w,f

def FilterFunctionPlot(group):
    f_whale = np.linspace(fminplot,np.log10(fmaxplot),fsteppllot)
    f_whale=10.0**f_whale
    
    weight,filterfunc = WeightFunction(group,f_whale)
    
    plt.figure()
    plt.plot(f_whale,filterfunc)
    plt.xscale('log')
    plt.grid()
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Filter Function [ ]')
    plt.tight_layout()
    plt.savefig('Seismics_Filter',dpi=300)
    

#%% Geometrical Spreading

def GeometricalSpreading(assumption,I_ref,x,dBred=14,waterdepth=20, safety=4.75):
    if assumption == 'Spherical':
        I = Distance2Intensity_Spherical(I_ref,x)
    elif assumption == 'Elmer':
        I = Distance2Intensity_Elmer(I_ref,x)
    elif assumption == 'Custom':
        I = Distance2Intensity_Custom(I_ref,x,dBred)
    elif assumption == 'ShallowWater':
        I = Distance2Intensity_ShallowWater(I_ref,x,waterdepth, safety)
    elif assumption == 'Cylindrical':
        I = Distance2Intensity_Cylindrical(I_ref,x)
    return I

def Distance2Intensity_Spherical(I_ref,x):
    """
    Helper function to calculate the sound intensity depending on the distance.
    Assuming spherical spreading.
    Inverse function of Intensity2Distance

    Parameters
    ----------
    I_ref : Float
        Reference sound intensity level [dB (re 1$\,\mu$Pa peak)].
    x : Float
        Distance to sound source [m].

    Returns
    -------
    I : Float
        Intensity [dB (re 1$\,\mu$Pa peak)]] at certain distance.

    """    
    I=np.ones_like(x)*I_ref
    x=np.abs(x)
    if np.shape(x)==():
        I=I_ref+20*np.log10(1/x)
    else:       
        I[x>0]=I_ref+20*np.log10(1/x[x>0])
    return I

def Distance2Intensity_Elmer(I_ref,x):
    """
    Helper function to calculate the sound intensity depending on the distance.
    For BfN Formula Elmer et al 2007.

    Parameters
    ----------
    I_ref : Float
        Reference sound intensity level [dB (re 1$\,\mu$Pa peak)].
    x : Float
        Distance to sound source [m].
    dBred : Float
        Custom value for geometrical spreading e.g. 14 dB

    Returns
    -------
    I : Float
        Intensity [dB (re 1$\,\mu$Pa peak)]] at certain distance.

    """    
    I=np.ones_like(x)*I_ref
    x=np.abs(x)
    if np.shape(x)==():
        I=I_ref-(14.0+x*0.0002)*np.log10(x)
    else: 
        I[x>0]=I_ref-(14.0+x*0.0002)*np.log10(x[x>0])
    return I


def Distance2Intensity_Custom(I_ref,x,dBred):
    """
    Helper function to calculate the sound intensity depending on the distance.
    A custom value of geometrical spreading can be chosen (20: Spherical,
    10: Cylindrical, 14: Modelling shallow water).

    Parameters
    ----------
    I_ref : Float
        Reference sound intensity level [dB (re 1$\,\mu$Pa peak)].
    x : Float
        Distance to sound source [m].
    dBred : Float
        Custom value for geometrical spreading e.g. 14 dB

    Returns
    -------
    I : Float
        Intensity [dB (re 1$\,\mu$Pa peak)]] at certain distance.

    """    
    I=np.ones_like(x)*I_ref
    x=np.abs(x)
    if np.shape(x)==():
        I=I_ref+dBred*np.log10(1/x)
    else: 
        I[x>0]=I_ref+dBred*np.log10(1/x[x>0])
    return I

def Distance2Intensity_ShallowWater(I_ref,x,waterdepth, safety):
    """
    Helper function to calculate the sound intensity depending on the distance.
    The geometrical spreading can be adjusted for shallow water scenarios.
    Formula taken from:
    Duncan, A. J. & Parsons, M. J. G. How Wrong Can You Be? Can a Simple 
    Spreading Formula Be Used to Predict Worst-Case Underwater Sound Levels? 
    8 (2011).


    Parameters
    ----------
    I_ref : Float
        Reference sound intensity level [dB (re 1$\,\mu$Pa peak)].
    x : Float
        Distance to sound source [m].
    waterdepth : Float
        Custom value for geometrical spreading e.g. 14 dB
    safety : float
        Safety threshold (st) for mean intensity depending on probabilty of 
        exceedence (p):
            st(p=1)    =-15 dB
            st(p=0.357)= 0 dB 
            st(p=0.05) = 4.75 => significant
            st(p=0.01) = 6.6 dB => very significant
            st(p=0.001)= 8.35 dB => highly significant

    Returns
    -------
    I : Float
        Intensity [dB (re 1$\,\mu$Pa peak)]] at certain distance.

    """    
    I=np.ones_like(x)*I_ref
    x=np.abs(x)
    if np.shape(x)==():
        I=I_ref+10*np.log10(1/x)+10*np.log10(2/waterdepth)+safety
    else:
        I[x>0]=I_ref+10*np.log10(1/x[x>0])+10*np.log10(2/waterdepth)+safety
        
    return I



def Distance2Intensity_Cylindrical(I_ref,x):
    """
    Helper function to calculate the sound intensity depending on the distance.
    Assuming cylinidrical spreading.
    Inverse function of Intensity2Distance

    Parameters
    ----------
    I_ref : Float
        Reference sound intensity level [dB (re 1$\,\mu$Pa peak)].
    x : Float
        Distance to sound source [m].

    Returns
    -------
    I : Float
        Intensity [dB (re 1$\,\mu$Pa peak)]] at certain distance.

    """  
    I=np.ones_like(x)*I_ref
    x=np.abs(x)
    if np.shape(x)==():
        I=I_ref+10*np.log10(1/x)
    else:
        I[x>0]=I_ref+10*np.log10(1/x[x>0])
    return I


def LoadKraken(model,refdepth):
    data = loadmat('../Data/SpreadingModels/'+model+'.mat')
    r_d = data['Pos']['r'][0][0][0]['depth'][0][:]
    r_r = data['Pos']['r'][0][0][0]['range'][0][:]
    geom = 20*np.log10(np.abs(data['tlt'][0,:,:]))
    mask = [r_d==refdepth][0]
    spreading = geom[mask[:,0],:].T
    
    return np.append(r_r,spreading,axis=1)

def SpreadingPlot(spreading,spreading_dBred,spreading_wd, spreading_safety,kraken_count,kraken_model,kraken_depth):
    dist = np.arange(1.0,2000.0,1.0)
    
    plt.figure()
    plt.plot(dist,GeometricalSpreading(spreading,0,dist,spreading_dBred,spreading_wd, spreading_safety),label='Used approximation: '+spreading)
    if spreading != 'Spherical':
        plt.plot(dist,GeometricalSpreading('Spherical',0,dist,spreading_dBred,spreading_wd, spreading_safety),label='Spherical')
    if spreading != 'ShallowWater':
        plt.plot(dist,GeometricalSpreading('ShallowWater',0,dist,spreading_dBred,spreading_wd, spreading_safety),label='Shallow Water '+str(spreading_wd) + ' m')
    if spreading != 'Elmer':
        plt.plot(dist,GeometricalSpreading('Elmer',0,dist,spreading_dBred,spreading_wd, spreading_safety),label='Elmer')
    if spreading != 'Cylindrical':
        plt.plot(dist,GeometricalSpreading('Cylindrical',0,dist,spreading_dBred,spreading_wd, spreading_safety),label='Cylindrical')
    if kraken_count >0:
        for kraken_step in range(0,kraken_count):
            kraken = LoadKraken(kraken_model[kraken_step],kraken_depth)
            plt.plot(kraken[:,0],kraken[:,1],'--',label='KRAKEN-model '+ kraken_model[kraken_step])
    plt.legend()
    plt.axis('tight')
    plt.grid()
    plt.xlabel('Distance [m]')
    plt.ylabel('Geometrical Spreading [dB]')
    plt.tight_layout()
    plt.savefig('Seismics_Spreading',dpi=300)
    
#%% Absorption
def SimpleAbsorption(sea,f):
    """
    Calculation of the absorption in different seas based on: 
    Ainslie, M.A., McColm, J.G., 1998. 
    A simplified formula for viscous and chemical absorption in sea water. 
    The Journal of the Acoustical Society of America 103, 
    1671–1672. https://doi.org/10.1121/1.421258
.

    Parameters
    ----------
    sea : String
        Selection of the sea: baltic, pacific, red, arctic, custom (User input).
    f : Float
        Frequency [Hz].

    Returns
    -------
    abs_sea : Float
        Absorption as a function of the frequency [dB/m].

    """
    if sea =='Baltic':
        S       = 8.0 #Salinity [psu]
        T       = 4.0 #Temperature [Deg C]
        pH      = 7.9  #Acidity [pH]
        D       = 0.0 #Depth [km]
    elif sea =='North': #For HeXXX very shallow :)
        S       = 33 #Salinity [psu]
        T       = 4 #Temperature [Deg C]
        pH      = 8.0  #Acidity [pH]
        D       = 0.0 #Depth [km]
    elif sea =='Pacific':
        S       = 34 #Salinity [psu]
        T       = 4 #Temperature [Deg C]
        pH      = 7.7  #Acidity [pH]
        D       = 1.0 #Depth [km]
    elif sea =='Red':
        S       = 40 #Salinity [psu]
        T       = 22.0 #Temperature [Deg C]
        pH      = 8.2  #Acidity [pH]
        D       = 0.2 #Depth [km]
    elif sea == 'Arctic':
        S       = 30 #Salinity [psu]
        T       = -1.5 #Temperature [Deg C]
        pH      = 8.2  #Acidity [pH]
        D       = 0 #Depth [km]
    elif sea == 'Custom':
        S       = float(input("Salinity in psu:\n")) #Salinity [psu]
        T       = float(input("Temperature in deg Celsius:\n")) #Temperature [Deg C]
        pH      = float(input("Acidity in pH:\n"))  #Acidity [pH]
        D       = float(input("(Mean)Depth in km:\n")) #Depth [km]
    
    f_boric = 0.78*np.sqrt(S/35)*np.exp(T/26) #Relaxation Frequency Boric acid [Hz]
    f_mag   = 42*np.exp(T/17) #Relaxation Frequency Magnesium sulphate [Hz]
    f_calc  = f/1000
    
    #Calcuation of the absorption [dB/km]
    abs_sea = 0.106*(f_boric*f_calc**2)/(f_boric**2+f_calc**2)*np.exp((pH-8)/(0.56))+0.52*(1+T/43)*(S/35)*(f_mag*f_calc**2)/(f_mag**2+f_calc**2)*np.exp(-D/6)+0.00049*f_calc**2*np.exp(-(T/27+D/17))
    
    abs_sea = abs_sea/1000 #Conversion to [dB/m]
    return abs_sea


def AbsorptionPlot(sea):
    f_whale = np.linspace(fminplot,np.log10(fmaxplot),fsteppllot)
    f_whale=10.0**f_whale
    
    abs_sea = SimpleAbsorption(sea,f_whale)
    
    plt.figure()
    plt.plot(f_whale,abs_sea*1000)
    plt.xscale('log')
    plt.grid()
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Absorption [dB/km]')
    plt.tight_layout()
    plt.savefig('Seismics_Absorption',dpi=300)
    
#%%

def peakSPL(t,p):
    """
    Calculation peak Sound Pressure Level [dB re 1 mu Pa] as given in 
    Southall, B. L. et al. Marine Mammal Noise Exposure Criteria: Updated 
    Scientific Recommendations for Residual Hearing Effects. Aquat Mamm 45, 
    125–232 (2019).

    Parameters
    ----------
    t : float
        time vector [s].
    p : float
        pressure vector [Pa].

    Returns
    -------
    SPL : float
        Peak Sound Pressure Level [dB re 1 mu Pa].
    T : float
        Dominant period = time difference global maximum to global minimum [s].

    """
    #Determination Sound pressure level peak to peak
    SPL=20.0*np.log10(np.max(np.abs(p))/10.0**(-6.0))
    return SPL


def SignalSpreading(general_input,wavelet_t,FFT,FFT_p,f,dist):
    #Calculate the absorption
    abs_f = -SimpleAbsorption(general_input['abs_model'],f)
    abs_f_d = abs_f*dist
    filter_abs = 10**(abs_f_d/20)
    
    #Calculate the geometrical spreading
    I_ref = np.float(0)
    GeomSpread_dB = GeometricalSpreading(general_input['spreading'],I_ref,dist,\
        general_input['spreading_dBred'],\
        general_input['spreading_wd'],\
        general_input['spreading_safety'])       
    GeomSpread_factor = 10**(GeomSpread_dB/20)
    

    t_filtered      = wavelet_t[wavelet_t<0.3]
    FFT_filter      = GeomSpread_factor*filter_abs*np.abs(FFT)
    p_filtered      = irfft(FFT_filter*np.exp(1j*FFT_p))
    p_filtered      = p_filtered[wavelet_t<0.3]

    SPL = peakSPL(t_filtered,p_filtered)
  
    return SPL

def SPL_Loop(general_input):
    
    lim_file = open("../Data/FilterFunctions/SoundLimits.yaml",'r')
    lim_input = yaml.safe_load(lim_file)
    lim_file.close()
    
    #Import the sound limits
    pts_limit = lim_input['PTS']['SPL'][general_input['group']]
    tts_limit = lim_input['TTS']['SPL'][general_input['group']]
    
    #Import the wavelet
    wavelet_t,wavelet_p,f,_,_,_,_,_,_ = \
        GetSourceInfo(general_input['source'],general_input['group'],plotbool=False)
    np.seterr(divide = 'ignore') 

    #Fourier transformation of the wavelets
    dt       = np.unique(np.round(np.diff(wavelet_t),decimals=15))
    f        = rfftfreq(wavelet_t.size, dt)
    FFT      = rfft(wavelet_p)
    FFT_p    = np.angle(FFT) #Phase
    
         
    #Make a distance vector
    dist = np.arange(distmin,distmax,distste)
    spl = np.zeros_like(dist)
    
    for dind,dstep in enumerate(dist):
        spl[dind] = SignalSpreading(general_input,wavelet_t,FFT,FFT_p,f,dstep)
        
    pts_limit_dist = np.min(dist[spl<=pts_limit])
    tts_limit_dist = np.min(dist[spl<=tts_limit])
    
    plt.figure()
    plt.plot(dist,spl,'-k',label='SPL')
    plt.plot([distmin,distmax],[pts_limit,pts_limit],'-r',label='PTS Limit')
    plt.plot([distmin,distmax],[tts_limit,tts_limit],'--r',label='TTS Limit')
    plt.legend()
    plt.xlabel('Distance [m]')
    plt.ylabel('SPL$_{uw,0-peak}$ [dB re 1 $\mu$ Pa')
    plt.xlim([distmin,distmax])
    plt.tight_layout()
    plt.grid()
    plt.savefig('Seismics_SPLdist',dpi=300)
    
    return pts_limit_dist,tts_limit_dist

def SignalWeighting(general_input,wavelet_t,FFT,FFT_p,f,dist):
    #Calculate the absorption
    abs_f = -SimpleAbsorption(general_input['abs_model'],f)
    abs_f_d = abs_f*dist
    filter_abs = 10**(abs_f_d/20)
    
    #Calculate the geometrical spreading
    I_ref = np.float(0)
    GeomSpread_dB = GeometricalSpreading(general_input['spreading'],I_ref,dist,\
        general_input['spreading_dBred'],\
        general_input['spreading_wd'],\
        general_input['spreading_safety'])       
    GeomSpread_factor = 10.0**(GeomSpread_dB/20.0)
    
    #Filter function
    _,filterfunc = WeightFunction(general_input['group'],f)

    t_filtered      = wavelet_t[wavelet_t<0.3]
    FFT_filter      = GeomSpread_factor*filter_abs*filterfunc*np.abs(FFT)
    p_filtered      = irfft(FFT_filter*np.exp(1j*FFT_p))
    p_filtered      = p_filtered[wavelet_t<0.3]
  
    return p_filtered,t_filtered

def SEL_Loop(general_input):
    
    lim_file = open("../Data/FilterFunctions/SoundLimits.yaml",'r')
    lim_input = yaml.safe_load(lim_file)
    lim_file.close()
    
    #Import the sound limits
    pts_limit = lim_input['PTS']['SEL'][general_input['group']]
    tts_limit = lim_input['TTS']['SEL'][general_input['group']]
    
    #Import the wavelet
    wavelet_t,wavelet_p,f,_,_,_,_,_,_ = \
        GetSourceInfo(general_input['source'],general_input['group'],plotbool=False)
    np.seterr(divide = 'ignore') 

    #Fourier transformation of the wavelets
    dt       = np.unique(np.round(np.diff(wavelet_t),decimals=15))
    f        = rfftfreq(wavelet_t.size, dt)
    FFT      = rfft(wavelet_p)
    FFT_p    = np.angle(FFT) #Phase

    #Line parameters
    whale_y        = np.arange(seldistmin,seldistmax,seldistste)
    profile_dur    = general_input['profile_dur']*60.0*60.0
    profile_spd    = general_input['profile_spd']*1852.0/60.0/60.0
    profile_len    = profile_dur*profile_spd
    profile_sr     = general_input['profile_sr']
    profile_ds     = profile_spd*profile_sr
    profile_x      = np.arange(-profile_len/2,profile_len/2,profile_ds)

    
    sel         = np.zeros_like(whale_y)    
    
    for whale_y_ind, whale_y_step in enumerate(whale_y):
        print(whale_y_step)
        for profile_x_ind, profile_x_step in enumerate(profile_x):
            dist = np.sqrt(whale_y_step**2+profile_x_step**2)
            p_filtered,_ = SignalWeighting(general_input,wavelet_t,FFT,FFT_p,f,dist)
            sel[whale_y_ind] = sel[whale_y_ind]+np.sum(p_filtered**2.)
     
    sel         = 10.0*np.log10(sel*dt*10.0**12.0)
 
    plt.figure()
    plt.plot(whale_y,sel,'-k',label='SEL')
    plt.plot([distmin,distmax],[pts_limit,pts_limit],'-r',label='PTS Limit')
    plt.plot([distmin,distmax],[tts_limit,tts_limit],'--r',label='TTS Limit')
    plt.legend()
    plt.xlabel('Distance [m]')
    plt.ylabel('SEL$_{w}$ [dB re 1 $\mu$ Pa$^2$ s]')
    plt.xlim([np.min(whale_y),np.max(whale_y)])
    plt.tight_layout()
    plt.grid()
    plt.savefig('Seismics_SELdist',dpi=300)
    
    try:
        pts_limit_dist = np.min(whale_y[sel<=pts_limit])
        tts_limit_dist = np.min(whale_y[sel<=tts_limit])
    except:
        pts_limit_dist = 0.0
        tts_limit_dist = 0.0
    
    return pts_limit_dist,tts_limit_dist    

def splrms125ms(p,t):
    p0 = 10.0*10**(-6)
    p_c = p[t<0.125]
    qmw = np.sqrt(np.mean(p_c**2)/p0**2)
    Leq = 20.0*np.log10(qmw)
    return Leq

def SPLrms_Loop(general_input):
    
    lim_file = open("../Data/FilterFunctions/SoundLimits.yaml",'r')
    lim_input = yaml.safe_load(lim_file)
    lim_file.close()

    #Import the sound limits
    b_limit = lim_input['Behaviour']['SPLrms'][general_input['group']]

    
    #Import the wavelet
    wavelet_t,wavelet_p,f,_,_,_,_,_,_ = \
        GetSourceInfo(general_input['source'],general_input['group'],plotbool=False)
    np.seterr(divide = 'ignore') 

    #Fourier transformation of the wavelets
    dt       = np.unique(np.round(np.diff(wavelet_t),decimals=15))
    f        = rfftfreq(wavelet_t.size, dt)
    FFT      = rfft(wavelet_p)
    FFT_p    = np.angle(FFT) #Phase
    
         
    #Make a distance vector
    dist = np.arange(distmin,distmax,distste)
    spl = np.zeros_like(dist)
    
    for dind,dstep in enumerate(dist):
        p_filtered,t_filtered = SignalWeighting(general_input,wavelet_t,FFT,FFT_p,f,dstep)
        spl[dind] = splrms125ms(p_filtered,t_filtered)
        
    b_limit_dist = np.min(dist[spl<=b_limit])
    
    plt.figure()
    plt.plot(dist,spl,'-k',label='SPL')
    plt.plot([distmin,distmax],[b_limit,b_limit],'--r',label='Limit Behavioral Effect')
    plt.legend()
    plt.xlabel('Distance [m]')
    plt.ylabel('SPL$_{W,RMS}$ [dB re 1 $\mu$ Pa]')
    plt.xlim([distmin,distmax])
    plt.tight_layout()
    plt.grid()
    plt.savefig('Seismics_SplRmsdist',dpi=300)
    
    return b_limit_dist

#%% Hydroacoustics

def SingleBeam_Plotting(yaml_file):
    #Open the parameter file and determine all parameters of interest
    yaml_file = open(yaml_file)
    yaml_input = yaml.safe_load(yaml_file)
    yaml_file.close()
    spl = yaml_input['SPL']
    #v=yaml_input['v']
    f=yaml_input['f']
    phi = yaml_input['phi']/180*np.pi
    dtheta = yaml_input['dtheta']/180*np.pi
    sigdur = yaml_input['sigdur']
    group = yaml_input['group']
    td = yaml_input['towingdepth']
    sea = yaml_input['abs_model']
    dBred = yaml_input['spreading_dBred']
    geom= yaml_input['spreading']
    depth = yaml_input['spreading_wd']
    safetydb = yaml_input['spreading_safety']
    profile_dur = yaml_input['profile_dur']*60*60
    profile_spd = yaml_input['profile_dur']*1852/60/60
    profile_len = profile_dur*profile_spd
    profile_sr  = 1/yaml_input['profile_sf']
    profile_ds  = profile_spd*profile_sr
    
    #Import the noise limits
    lim_file = open("../Data/FilterFunctions/SoundLimits.yaml",'r')
    lim_input = yaml.safe_load(lim_file)
    lim_file.close()
    pts_spl = lim_input['PTS']['SPL'][yaml_input['group']]
    tts_spl = lim_input['TTS']['SPL'][yaml_input['group']]
    pts_sel = lim_input['PTS']['SEL'][yaml_input['group']]
    tts_sel = lim_input['TTS']['SEL'][yaml_input['group']]
    b_spl   = lim_input['Behaviour']['SPLrms'][yaml_input['group']]
    
    
    #Simple Directivity Plot 1D
    along_a = np.arange(-80.0,80.0,0.1) #Angles
    along_a_rad = along_a/180*np.pi #Angles in radians
    e = np.pi*along_a_rad/dtheta #Element Transducer Directivity Parameter
    dfy = 20.0*np.log10(np.abs(np.sin(e)/e)) #Element Transducer Directivity
    a = np.pi/phi*(np.sin(along_a_rad)-np.sin(0)) #Beam Pattern Parameter
    dfx = 20*np.log10(np.abs(np.sin(a)/a)) #Beam Pattern
    #Plotting
    plt.figure()
    plt.plot(along_a,dfy,'-k',label='Element Transducer Pattern')
    plt.plot(along_a,dfx,'--b',label='Beam Pattern')
    plt.plot(along_a,dfx+dfy,'-r',label='Combination')
    plt.grid()
    plt.legend()
    plt.ylabel('Directivity Pattern [dB]')
    plt.xlabel('Along-track angle [deg]')
    plt.tight_layout()
    plt.savefig('SingleBeam_Directivity1D', dpi=300)
    
    # Underwater Noise Measures
    #Spatial Coordinates
    x = np.arange(-profile_len/2,profile_len/2,profile_ds)
    y = np.arange(hydro_seldistmin,hydro_seldistmax,hydro_seldistste)
    xgrid, ygrid = np.meshgrid(x,y)
    dist = np.sqrt(xgrid**2+ygrid**2+(depth-td)**2) #Distance to point
    dist_h = np.sqrt(xgrid**2+ygrid**2) #Horizontal distance
    theta = np.abs(np.arctan(dist_h/(depth-td))+1e-19) #Angle
    #Directivity
    e = np.pi*theta/dtheta #Element Transducer Directivity Parameter
    dfy = 20.0*np.log10(np.abs(np.sin(e)/e)) #Element Transducer Directivity
    a = np.pi/phi*(np.sin(theta)-np.sin(0)) #Beam Pattern Parameter
    dfx = 20*np.log10(np.abs(np.sin(a)/a)) #Beam Pattern 
    #Corrections
    ed = 10*np.log10(sigdur) #Event duration for dB calculation
    fw,_ = WeightFunction(group,f) #Weight in dB for functional hearing group
    tl = GeometricalSpreading(geom,0.0,dist,dBred,depth,safetydb) #Transmission loss by geometric spreading
    absorp = -SimpleAbsorption(sea,f)*dist #Absorption
    #Calculation
    spl_peak = spl+tl+dfx+dfy+absorp
    spl_rms  = spl+tl-3.+fw+dfx+dfy+absorp
    sel_shot = spl-3.+tl+ed+fw+dfx+dfy+absorp
    sel_shot = 10**(sel_shot/10)  
    sel_cum = 10.*np.log10(np.sum(sel_shot,axis=1))
    #Plotting
    plt.figure(figsize=(10,8))
    plt.subplot(3,1,1)
    plt.pcolormesh(x,y,spl_peak,shading='auto')
    c=plt.colorbar()
    c.set_label('SPL [dB re 1$\mu$Pa]')
    plt.contour(x,y,spl_peak,[pts_spl,tts_spl],colors=['red','red'],linestyles=['--','-'])
    plt.xlabel('Along Track distance [m]')
    plt.ylabel('Across Track distance [m]')
    plt.subplot(3,1,2)
    plt.pcolormesh(x,y,spl_rms,shading='auto',vmin=80)
    c=plt.colorbar()
    c.set_label('RMS SPL [dB re 1$\mu$Pa]')
    plt.contour(x,y,spl_rms,[b_spl],colors=['red'],linestyles=['-'])
    plt.xlabel('Along Track distance [m]')
    plt.ylabel('Across Track distance [m]')
    plt.subplot(3,1,3)
    plt.plot(y,np.ones_like(y)*pts_sel,'-r',label = 'PTS limit')
    plt.plot(y,np.ones_like(y)*tts_sel,'--r',label = 'TTS limit')
    plt.plot(y,sel_cum,'-k')
    plt.grid()
    plt.xlim(0,1400)
    plt.legend()
    plt.xlabel('Across Track distance [m]')
    plt.ylabel('SEL [dB re 1$\mu$Pa$^2$]')
    plt.tight_layout()
    plt.savefig('SingleBeam_NoiseLimits', dpi=300)
    
    try:
        pts_spl_d = np.nanmax(dist[spl_peak>pts_spl])
    except:
        pts_spl_d = 0.0
    try:
        tts_spl_d = np.nanmax(dist[spl_peak>tts_spl])
    except:
        tts_spl_d = 0.0
    try:
        pts_sel_d = np.nanmax(y[sel_cum>pts_sel])
    except:
        pts_sel_d = 0.0
    try:
        tts_sel_d = np.nanmax(y[sel_cum>tts_sel])
    except:
        tts_sel_d = 0.0
    try:
        b_d = np.nanmax(dist[spl_rms>b_spl])
    except:
        b_d       = 0.0
        
    return pts_spl_d ,tts_spl_d, pts_sel_d, tts_sel_d, b_d
    
#%%
def MultiBeam_Plotting(yaml_file):
    #Open the parameter file and determine all parameters of interest
    yaml_file = open(yaml_file)
    yaml_input = yaml.safe_load(yaml_file)
    yaml_file.close()
    spl = yaml_input['SPL']
    f=yaml_input['f']
    phi = yaml_input['phi']/180*np.pi
    dtheta = yaml_input['dtheta']/180*np.pi
    sigdur = yaml_input['sigdur']
    td = yaml_input['towingdepth']
    group = yaml_input['group']
    sea = yaml_input['abs_model']
    dBred = yaml_input['spreading_dBred']
    geom= yaml_input['spreading']
    depth = yaml_input['spreading_wd']
    safetydb = yaml_input['spreading_safety']
    profile_dur = yaml_input['profile_dur']*60*60
    profile_spd = yaml_input['profile_dur']*1852/60/60
    profile_len = profile_dur*profile_spd
    profile_sr  = 1/yaml_input['profile_sf']
    profile_ds  = profile_spd*profile_sr
    
    #Import the noise limits
    lim_file = open("../Data/FilterFunctions/SoundLimits.yaml",'r')
    lim_input = yaml.safe_load(lim_file)
    lim_file.close()
    pts_spl = lim_input['PTS']['SPL'][yaml_input['group']]
    tts_spl = lim_input['TTS']['SPL'][yaml_input['group']]
    pts_sel = lim_input['PTS']['SEL'][yaml_input['group']]
    tts_sel = lim_input['TTS']['SEL'][yaml_input['group']]
    b_spl   = lim_input['Behaviour']['SPLrms'][yaml_input['group']]
    
    
    #Simple Directivity Plot 1D
    along_a = np.arange(-80.0,80.0,0.1) #Angles
    along_a_rad = along_a/180*np.pi+1e-15 #Angles in radians
    e = np.pi*along_a_rad/dtheta #Element Transducer Directivity Parameter
    dfy = 20.0*np.log10(np.abs(np.sin(e)/e)) #Element Transducer Directivity
    a = np.pi/phi*(np.sin(along_a_rad)-np.sin(0)) #Beam Pattern Parameter
    dfx = 20*np.log10(np.abs(np.sin(a)/a)) #Beam Pattern
    #Plotting
    plt.figure()
    plt.plot(along_a,dfy,'-k',label='Element Transducer Pattern')
    plt.plot(along_a,dfx,'--b',label='Beam Pattern')
    plt.plot(along_a,dfx+dfy,'-r',label='Combination')
    plt.grid()
    plt.legend()
    plt.ylabel('Directivity Pattern [dB]')
    plt.xlabel('Along-track angle [deg]')
    plt.tight_layout()
    plt.savefig('MultiBeam_Directivity1D', dpi=300)
    
    
    # # Directivity Plot 2D
    # x = np.arange(-hydro_beamx,hydro_beamx,hydro_beamdx) #X-Vector
    # y = np.arange(0,hydro_beamy,hydro_beamdx) #Y-Vector
    # xgrid, ygrid = np.meshgrid(x,y) #Grid
    # dist = np.sqrt(ygrid**2+xgrid**2+(depth-td)**2) #Distance to point
    # dist_h = np.sqrt(ygrid**2+xgrid**2) #Distance
    # # Directivity
    # theta = np.abs(np.arctan(dist_h/(depth-td))+1e-19) #Angle point
    # e = np.pi*theta/dtheta #Element Transducer Directivity Parameter
    # dfy = 20.0*np.log10(np.abs(np.sin(e)/e)) #Element Transducer Directivity
    # beta = np.abs(np.arctan(xgrid/dist)+1e-19)
    # a = np.pi/phi*(np.sin(beta)-np.sin(0)) #Beam Pattern Parameter
    # dfx = 20*np.log10(np.abs(np.sin(a)/a)) #Beam Pattern 
    # #Plotting
    # plt.figure()
    # plt.subplot(1,3,1)
    # plt.pcolormesh(x,y,dfx,vmin=-30,vmax=0,shading='auto')
    # plt.subplot(1,3,2)
    # plt.pcolormesh(x,y,dfy,vmin=-30,vmax=0,shading='auto')
    # plt.subplot(1,3,3)
    # plt.pcolormesh(x,y,dfx+dfy,vmin=-30,vmax=0,shading='auto')
    # plt.colorbar()
    # plt.savefig('MultiBeam_Directivity2D', dpi=300)
    
    
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
    tl = GeometricalSpreading(geom,0.0,dist,dBred,depth,safetydb) #Transmission loss by geometric spreading
    absorp = -SimpleAbsorption(sea,f)*dist #Absorption
    #Calculation
    spl_peak = spl+tl+dfx+dfy+absorp
    spl_rms  = spl+tl-3.+fw+dfx+dfy+absorp
    sel_shot = spl-3.+tl+ed+fw+dfx+dfy+absorp
    sel_shot = 10**(sel_shot/10)  
    sel_cum = 10.*np.log10(np.sum(sel_shot,axis=1))
    #Plotting
    plt.figure(figsize=(10,8))
    plt.subplot(3,1,1)
    plt.pcolormesh(x,y,spl_peak,shading='auto')
    c=plt.colorbar()
    c.set_label('SPL [dB re 1$\mu$Pa]')
    plt.contour(x,y,spl_peak,[pts_spl,tts_spl],colors=['red','red'],linestyles=['--','-'])
    plt.xlabel('Along Track distance [m]')
    plt.ylabel('Across Track distance [m]')
    plt.subplot(3,1,2)
    plt.pcolormesh(x,y,spl_rms,shading='auto',vmin=80)
    c=plt.colorbar()
    c.set_label('RMS SPL [dB re 1$\mu$Pa]')
    plt.contour(x,y,spl_rms,[b_spl],colors=['red'],linestyles=['-'])
    plt.xlabel('Along Track distance [m]')
    plt.ylabel('Across Track distance [m]')
    plt.subplot(3,1,3)
    plt.plot(y,np.ones_like(y)*pts_sel,'-r',label = 'PTS limit')
    plt.plot(y,np.ones_like(y)*tts_sel,'--r',label = 'TTS limit')
    plt.plot(y,sel_cum,'-k')
    plt.grid()
    plt.xlim(0,1400)
    plt.legend()
    plt.xlabel('Across Track distance [m]')
    plt.ylabel('SEL [dB re 1$\mu$Pa$^2$]')
    plt.tight_layout()
    plt.savefig('MultiBeam_NoiseLimits', dpi=300)
    
    try:
        pts_spl_d = np.nanmax(dist[spl_peak>pts_spl])
    except:
        pts_spl_d = 0.0
    try:
        tts_spl_d = np.nanmax(dist[spl_peak>tts_spl])
    except:
        tts_spl_d = 0.0
    try:
        pts_sel_d = np.nanmax(y[sel_cum>pts_sel])
    except:
        pts_sel_d = 0.0
    try:
        tts_sel_d = np.nanmax(y[sel_cum>tts_sel])
    except:
        tts_sel_d = 0.0
    try:
        b_d = np.nanmax(dist[spl_rms>b_spl])
    except:
        b_d       = 0.0
        
    return pts_spl_d ,tts_spl_d, pts_sel_d, tts_sel_d, b_d