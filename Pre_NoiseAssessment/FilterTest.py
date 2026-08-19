# -*- coding: utf-8 -*-
"""
Created on Tue Mar 11 16:04:24 2025

@author: nikol
"""

#%% 01) Import dependencies
#General dependencies
import numpy as np
import datetime
from pylatex import Document, Section, Subsection, \
     Figure, Alignat, PageStyle, Head, LineBreak, \
     simple_page_number, Foot, MiniPage, LargeText, MediumText, NewPage
from pylatex.utils import bold
import yaml

from scipy.fft import rfft,rfftfreq,irfft
import matplotlib.pyplot as plt
from scipy.io import loadmat
from scipy.interpolate import interp1d
from scipy.signal.windows import tukey

#Custom functions
import ImpactAssessment_Functions as iafunc

from matplotlib import rc
rc("pdf", fonttype=42)

#%%
yaml_file = "Seismics_Paremeters_SercelMicroGI_VHF.yaml"


general_file = open(yaml_file)
general_input = yaml.safe_load(general_file)
general_file.close()

#Make the wavelet plot
wavelet_t,wavelet_p,f,wavelet_fft_I,wavelet_fmin,wavelet_fmax,\
    wavelet_pulse_length,wavelet_SPL,wavelet_SEL \
    = iafunc.GetSourceInfo(general_input,general_input['group'])
    
drange = np.arange(1.0,5000.0,10.0)
maxp  = np.zeros_like(drange) 
for ind,dstep in enumerate(drange):
    print(dstep)
    p_filtered,t_filtered = iafunc.SignalWeighting(general_input,wavelet_t,wavelet_p,dstep)
    maxp[ind] = np.max(np.abs(p_filtered))
    
p0 = 10.0*10**(-6)
maxp = 20*np.log10(maxp/p0)
plt.figure()
plt.plot(drange,maxp)

#%%
group = general_input['group']

f_whale = np.linspace(0.1,np.log10(170000.),500)
f_whale=10.0**f_whale

t_whale=iafunc.threshold(group,f_whale)

plt.figure()
plt.plot(f_whale,t_whale,label='Hearing Threshold '+group)
plt.xlabel('Frequency [Hz]')
plt.ylabel(r"Sound Power" "\n" r"[dB re 1$\mu Pa$]")
plt.xscale('log')
plt.axis('tight')
#plt.ylim([0,220])

df = np.diff(f_whale)  # in Hz
df = np.append(df, df[-1])  # pad last value

E_threshold = np.sum(t_whale * df)  # µPa²·s
SEL_threshold = 10 * np.log10(E_threshold)



print(f"SEL_threshold = {SEL_threshold:.2f} dB re 1 µPa²·s")

#%%

# Example: Harbor porpoise audiogram (from Southall et al., 2007)
# Frequency (Hz), Threshold (dB re 1 µPa)
frequencies = np.array([1000, 2000, 5000, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000])
thresholds = np.array([50, 45, 40, 35, 30, 35, 40, 50, 60, 70, 80, 90, 100])  # dB re 1 µPa

# Convert dB to µPa² (energy)
threshold_energy = 10**(thresholds / 10)  # µPa²

# Frequency bandwidth (1/3-octave bands)
# Approximate Δf for each frequency
df = np.diff(frequencies)  # in Hz
df = np.append(df, df[-1])  # pad last value

# Compute total threshold energy (in µPa²·s)
E_threshold = np.sum(threshold_energy * df)  # µPa²·s

# Convert to SEL_threshold (dB re 1 µPa²·s)
SEL_threshold = 10 * np.log10(E_threshold)

plt.figure()
plt.plot(frequencies,thresholds)
print(f"SEL_threshold = {SEL_threshold:.2f} dB re 1 µPa²·s")


#%%
f_whale = np.linspace(0.1,np.log10(700000.),500)
f_whale=10.0**f_whale



#
dist =250.
zr=20.
#Calculate the geometrical spreading
abs_f = -iafunc.SimpleAbsorption(general_input['abs_model'],f)
abs_f_d = abs_f*dist
filter_abs = 10**(abs_f_d/20)
    
    
I_ref = 0.0
GeomSpread_dB = iafunc.GeometricalSpreading(general_input['spreading'],I_ref,dist,\
    general_input['spreading_dBred'],\
    general_input['spreading_wd'],\
    general_input['spreading_safety'],\
    general_input['kraken_count'],\
    general_input['kraken_model'],\
    general_input['kraken_depth'],\
    general_input['kraken_selection'])       
GeomSpread_factor = 10**(GeomSpread_dB/20)

#Calculate the directivity
#Get Water Depth [m]
zr = general_input['spreading_wd']  #z=> Depth; r => Receiver [m]
v = general_input['watervel']  #water velocity [m/s]
#Get Tow Depth [m]
source_yaml_input=general_input['source']
yaml_file = open("../Data/Sources/"+source_yaml_input+'_Specs.yaml','r')
yaml_input_spec = yaml.safe_load(yaml_file)
yaml_file.close()
zs  = yaml_input_spec[general_input['source']]['towdepth'] #z=> Depth; s => Source [m]
fd  = yaml_input_spec[general_input['source']]['fd'] #f=> frequency; d=> dominant [Hz]
#Calc directivity
d = np.sqrt(dist**2+zr**2)
direct_dB = iafunc.lloyds_mirror(v,fd,zr,zs,d)
direct = 10**(direct_dB/20)

fft_filter = GeomSpread_factor*filter_abs*direct

# #Fourier transform
wavelet_dt = np.unique(np.round(np.diff(wavelet_t),decimals=15))[0]
f        = rfftfreq(wavelet_t.size, wavelet_dt)
FFT      = rfft(wavelet_p*10**6)
FFT_A    = tukey(np.size(FFT),0.1)*np.abs(FFT)

FFT_I    = 20*np.log10(FFT_A/np.size(FFT)*fft_filter)



#
fig,ax = plt.subplots(1,1)
ax.plot(f_whale,iafunc.threshold(general_input['group'],f_whale),label='Hearing Threshold ')
ax.plot(f,FFT_I,label='Micro-GI')
ax.set_xlabel('Frequency [Hz]')
ax.set_ylabel(r'Sound Power [dB re 1$\mu Pa$/Hz]')
ax.set_xscale('log')
ax.axis('tight')
ax.set_xlim([1,100000])
ax.set_ylim([0,200])
ax.grid()
ax.legend()



#%%
yaml_file = "Seismics_Paremeters_SercelMicroGI_VHF.yaml"


general_file = open(yaml_file)
general_input = yaml.safe_load(general_file)
general_file.close()


#Import animal information
yaml_file = open("../Data/FilterFunctions/SoundLimits.yaml",'r')
yaml_input = yaml.safe_load(yaml_file)
yaml_file.close()
groupinfo = yaml_input['Description'][general_input['group']]


#Open the source config 
yaml_file = open("../Data/Sources/"+general_input['source']+'_Specs.yaml','r')
yaml_source = yaml.safe_load(yaml_file)
yaml_file.close()
sspecs  = yaml_source[general_input['source']]


#Make the wavelet plot
wavelet_t,wavelet_p,f,wavelet_fft_I,wavelet_fmin,wavelet_fmax,\
    wavelet_pulse_length,wavelet_SPL,wavelet_SEL \
    = iafunc.GetSourceInfo(general_input,general_input['group'])
      
lim_file = open("../Data/FilterFunctions/SoundLimits.yaml",'r')
lim_input = yaml.safe_load(lim_file)
lim_file.close()

#Import the sound limits
pts_limit = lim_input['PTS']['SPL'][general_input['group']]
tts_limit = lim_input['TTS']['SPL'][general_input['group']]


#Fourier transformation of the wavelets
dt       = np.unique(np.round(np.diff(wavelet_t),decimals=15))
f        = rfftfreq(wavelet_t.size, dt)
FFT      = rfft(wavelet_p)
FFT_p    = np.angle(FFT) #Phase

#%%   
#Make a distance vector
distmin=5.0
distmax=6.0
diste = 1.0

dstep = 20

#Calculate the absorption
abs_f = -iafunc.SimpleAbsorption(general_input['abs_model'],f)
abs_f_d = abs_f*dstep
filter_abs = 10**(abs_f_d/20)

#Calculate the geometrical spreading
I_ref = 0.0
GeomSpread_dB = iafunc.GeometricalSpreading(general_input['spreading'],I_ref,dstep,\
    general_input['spreading_dBred'],\
    general_input['spreading_wd'],\
    general_input['spreading_safety'],\
    general_input['kraken_count'],\
    general_input['kraken_model'],\
    general_input['kraken_depth'],\
    general_input['kraken_selection'])       
GeomSpread_factor = 10**(GeomSpread_dB/20)

#Calculate the directivity
#Get Water Depth [m]
zr = general_input['spreading_wd']  #z=> Depth; r => Receiver [m]
v = general_input['watervel']  #water velocity [m/s]
#Get Tow Depth [m]
source_yaml_input=general_input['source']
yaml_file = open("../Data/Sources/"+source_yaml_input+'_Specs.yaml','r')
yaml_input_spec = yaml.safe_load(yaml_file)
yaml_file.close()
zs  = yaml_input_spec[general_input['source']]['towdepth'] #z=> Depth; s => Source [m]
fd  = yaml_input_spec[general_input['source']]['fd'] #f=> frequency; d=> dominant [Hz]
#Calc directivity
d = np.sqrt(dstep**2+zr**2)
direct_dB = iafunc.lloyds_mirror(v,fd,zr,zs,d)
direct = 10**(direct_dB/20)

    #Filter function
if ((general_input['group']=='Fish') or\
    (general_input['group']=='ST') or \
    (general_input['group']=='Human')):
    filterfunc = np.ones_like(f)
else:
    _,filterfunc = iafunc.WeightFunction(general_input['group'],f)
        
# t_filtered      = wavelet_t[wavelet_t<0.45]
# FFT_filter      = GeomSpread_factor*filter_abs*np.abs(FFT)*direct
# p_filtered      = irfft(FFT_filter*np.exp(1j*FFT_p))
# p_filtered      = p_filtered[wavelet_t<0.45]

def arbitraryfilter(fft_f,fft_filter,t_original,s_original,plotbool=False):
    """
    Filter with a Arbitrary Frequency Response
    Procedure taken from https://www.dspguide.com/ch17/1.htm
    Zero Phase Filter.

    Parameters
    ----------
    fft_f : np.array
        Frequency vector of the Fourier Transform of the desired filter.
    fft_filter : np.array
        Fourier Transform of the desired filter. => Zero Phase Filter assumed.
    t_original : np.array
        time vector of the signal.
    s_original : np.array
        frequency vector of the signal.
    plotbool : TYPE, optional
        Decision on plot. The default is False.

    Returns
    -------
    s_filter : np.array
        Filtered signal.

    """

    
    #Calculate impulse response
    filter_t    = irfft(filter_fft+1j*np.zeros_like(filter_fft),n=np.size(t_original))
    print(np.shape(filter_t))
    #Shift+Truncate => Truncation to next even number down
    maxind      = int(np.floor(np.size(filter_t)/2))    
    filter_t_s  = np.zeros(maxind*2,dtype=np.complex128)
    filter_t_s[0:maxind] = filter_t[-maxind:]
    filter_t_s[maxind:] = filter_t[0:maxind]
    
    #Window
    filter_t_s_w = filter_t_s*tukey(maxind*2)
    
    #Convolve
    s_filter = np.convolve(s_original,filter_t_s_w,mode='same')

    
    if plotbool:
        
        fft_signal  = rfft(s_original)

        fig,[ax1,ax2]=plt.subplots(2,1) 
        
        ax1.plot(fft_f,np.abs(filter_fft),label='Filter Functsion')
        ax1.plot(fft_f,np.abs(fft_signal),label='Signal')
        ax1.legend()
        ax1.grid()
        ax1.set_xlabel('Frequency [Hz]')
        ax1.set_xlabel('Amplitude [Pa]')


        #ax2.plot(wavelet_t,filter_t,label='Raw Filter')
        #ax2.plot(wavelet_t[:maxind*2],filter_t_s,'label=Shifted and Truncated Filter')
        ax2.plot(wavelet_t[:maxind*2],filter_t_s_w,label='Windowed Filter')
        ax2.plot(wavelet_t,wavelet_p,label='Original Trace')
        ax2.plot(wavelet_t,s_filter,label='Filtered Trace')
        ax2.legend()
        ax2.grid()
        ax2.set_xlabel('Time [s]')
        ax2.set_xlabel('Amplitude [Pa]')
        plt.tight_layout()
        
    return s_filter

    
    
filter_fft        = GeomSpread_factor*filter_abs*direct#*filterfunc

wavelet_filtered = arbitraryfilter(f,filter_fft,wavelet_t,wavelet_p,plotbool=True)

#spl = iafunc.peakSPL(t_filtered,p_filtered)
    
     