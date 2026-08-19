# -*- coding: utf-8 -*-
"""
MarineEarsReporting.py

Functions to generate the HDF5 data report for the underwater noise registery 
of the BSH/MarineEars

Author:             Nikolas Römer-Stange
Initial Draft:      2021.10.27
Last Update:        2020.11.02

Dependencies:       Defined in section 01 of the script.

Units:              SI unless stated differently.

Nomenclature:       Defined upon first call of variable

References:         See subrutines/functions
    
Comments:           The data report is adjusted to the data requirements for 
                    seismic survey (2021-09-24)

"""

#%% Import of dependencies

import numpy as np
import h5py
from scipy.fft import rfft,rfftfreq
import datetime
import yaml
import pathlib
from segysak.segy import segy_loader
import pynmea2 as pn # Package to read nmea strings
import os
import copy
import pandas

#%% Navigation/ICES

def ImportNMEA(directory):
    """
    ImportNMEA reads all RMC strings of all *.nmea files in acertain directory 
    and returns time, lat and long in 15min intervals

    Parameters
    ----------
    directory : string
        (Relative) Path to the *.nmea files

    Returns
    -------
    nav : float array with three columns and variable number of lines
        first column is posix timestamp, second decimal degree latitude, 
        third is decimal degree longitude

    """
    #Initialisation
    oldtime = datetime.datetime(1900,1,1,0,0,0)
    nav = np.zeros((1,3),dtype=np.float64())
    
    for file in os.listdir(directory):
         filename = os.fsdecode(file)
         if filename.endswith(".nmea"): 
            deck_finput = open(directory+filename)
            for i, l in enumerate(deck_finput):
                try: #This nesting prevents the loop to be interrupted, when the line can not be read
                    msg  = pn.parse(l, check=True) #Try to read the line
                    if msg.sentence_type == 'RMC': #Continue if the line is a RMC string
                        currenttime = datetime.datetime.combine(msg.datestamp,msg.timestamp)
                        if ((currenttime.strftime('%M:%S.%f')=='00:00.000000')or\
                            (currenttime.strftime('%M:%S.%f')=='15:00.000000')or\
                            (currenttime.strftime('%M:%S.%f')=='30:00.000000')or\
                            (currenttime.strftime('%M:%S.%f')=='45:00.000000')or\
                            currenttime-oldtime>datetime.timedelta(minutes=15)):
                            #print(currenttime.strftime('%Y-%m-%d %H:%M:%S.%f'))
                            oldtime = copy.copy(currenttime)
                            nav = np.append(nav,np.array([[currenttime.timestamp(),msg.latitude,msg.longitude]]),axis=0)
                except ValueError: #This happens, when the line can not be read
                    print('     Error in string in line: ',i) # It is an error message
            
            
            deck_finput.close() #Close the input file   
    nav = np.delete(nav,0,0)
    return nav  

def LatCode(lat):
    latcode_1 = np.floor((lat-36.0)*2)+1
    latcode_1 = '%02.0f'%latcode_1
    return latcode_1

def LonCode(lon):
    loncode_1 = np.floor((lon+50)/10+1)
    loncode_1 = chr(int(loncode_1+64))
    loncode_2 = np.floor((lon+50)%10)
    if loncode_1=='A':
        loncode_2=loncode_2-5
    loncode_2 = '%1.0f'%loncode_2
    return loncode_1+loncode_2

def SubRect(lat,lon):
    sublat_1 = (lat%0.5)*60
    sublon_1 = (lon%1)*60
    
    if sublon_1>=40:
        sublon_2=6
    elif sublon_1>=20:
        sublon_2=3
    else:
        sublon_2=0
        
    if sublat_1>=20:
        sublat_2=1
    elif sublat_1>=10:
        sublat_2=2
    else:
        sublat_2=3
    return str(sublat_2+sublon_2)
            
def IcesSubtRect(lat,lon):
    if ((lat>=36)&(lat<=85.5)&(lon>=-44)&(lon<=68.5)):
        firstlat   = LatCode(lat)
        secondlong = LonCode(lon)
        thirdsub   = SubRect(lat,lon)
    else:
        print('Position ',lat,'N, ',lon,'E outside ICES Statistical rectangles grid.')
        firstlat=np.nan
        secondlong=np.nan
        thirdsub=np.nan
    return firstlat+secondlong+thirdsub

def WrapperICES(latarray,lonarray):
    ices = np.zeros_like(latarray,dtype='object_')
    for ind,step in enumerate(latarray):
        ices[ind]=IcesSubtRect(step,lonarray[ind])
    return ices

#%%
        
def ReadGeneralYaml(configfile):
    yaml_file       = open(configfile,'r')
    yaml_input      = yaml.safe_load(yaml_file)
    generalinput    = yaml_input['general']
    airguninput     = yaml_input['airguns']
    boomerinput     = yaml_input['boomer']
    sbpinput        = yaml_input['sub_bottom_prifiler']
    sparkerrinput   = yaml_input['sparker']

    return generalinput,airguninput,boomerinput,sbpinput,sparkerrinput

def WriteGeneralInfo(hf,generalinput):
    #Get current time
    now = datetime.datetime.now()
    
    #Write general info 
    hf.create_dataset('author', data=generalinput['author'])
    hf.create_dataset('comment', data=generalinput['comment'])
    hf.create_dataset('date_of_creation', data=now.strftime("%Y-%m-%d"))
    hf.create_dataset('journey_name', data=generalinput['journey_name'])
    hf.create_dataset('measurement_purpose', data=generalinput['measurement_purpose'])
    hf.create_dataset('measuring_institution', data=generalinput['measuring_institution'])
    hf.create_dataset('point_of_contact', data=generalinput['point_of_contact'])
    
    #Initialize a SPL level for later use in transects
    source_name = []
    source_SPL  = []
    
    return  source_name,source_SPL

#%% 
def GetSourceInfo(source_yaml_input,subgrouphdf5,iterationstepper):
    #Open the wavelet 
    yaml_file = open("../Data/Sources/"+source_yaml_input[iterationstepper]['ident']+'_Measurements.yaml','r')
    yaml_input = yaml.safe_load(yaml_file)
    measurements  = yaml_input[source_yaml_input[iterationstepper]['ident']]
    wavelet_t = np.array(measurements['t'])
    wavelet_p = np.array(measurements['p'])
    wavelet_dt = np.unique(np.round(np.diff(wavelet_t),decimals=15))[0]
    
    #Fourier transform
    f        = rfftfreq(wavelet_t.size, wavelet_dt)
    FFT      = rfft(wavelet_p)
    FFT_A    = np.abs(FFT)/FFT.size
    FFT_A    = FFT_A/np.max(FFT_A)
    FFT_I    = 20*np.log10(FFT_A/np.max(FFT_A))
    
    #Min+MaxFrequency
    f_temp = f[FFT_I>=np.max(FFT_I)-6]
    subgrouphdf5.create_dataset('frequency_min_theo',data=np.min(f_temp),dtype='f')
    subgrouphdf5.create_dataset('frequency_max_theo',data=np.max(f_temp),dtype='f')
    
    #Pulse Length
    wavesum = np.cumsum(wavelet_p**2)
    wavesum = wavesum/np.max(wavesum)*100
    pulse_length = np.max(wavelet_t[wavesum<97.5])-np.max(wavelet_t[wavesum<2.5])
    subgrouphdf5.create_dataset('pulse_length',data=pulse_length,dtype='f')     
    #SPL
    SPL=20.0*np.log10(np.max(np.abs(wavelet_p))/10.0**(-6.0))
    subgrouphdf5.create_dataset('source_level_Lpeak',data=SPL,dtype='f') 
    #SEL
    SEL = 10.0*np.log10(np.sum(wavelet_p**2.0)*wavelet_dt*10.0**12.0)
    subgrouphdf5.create_dataset('source_level_SEL',data=SEL,dtype='f')
    
    return SPL
        

def GetSourceExtendedInfo(source_yaml_input,subgrouphdf5,iterationstepper):
    #Import SEGY file
    segyfilename = pathlib.Path('../Data/ExampleProfiles/'+source_yaml_input[iterationstepper]['example_file'])
    P2D = segy_loader(segyfilename)   
    dt = P2D.sample_rate
    twt = np.array(P2D.twt)/1000
    cdp = np.array(P2D.cdp)
    data= np.array(P2D.data).T
    
    #Fourier transform
    f        = rfftfreq(twt.size, dt/10**3)
    FFT      = np.zeros((f.size,cdp.size))
    for step in range(cdp.size):
        FFT[:,step]      = np.abs(rfft(data[:,step]))/f.size
    FFT      = np.mean(FFT,1)
    FFT_A    = FFT/np.max(FFT)
    FFT_I    = 20*np.log10(FFT_A/np.max(FFT_A))

    #Min+MaxFrequency
    f_temp = f[FFT_I>=np.max(FFT_I)-6]
    subgrouphdf5.create_dataset('frequency_min_observed',data=np.min(f_temp),dtype='f')
    subgrouphdf5.create_dataset('frequency_max_observed',data=np.max(f_temp),dtype='f')
    
    #Pulse Length
    wavesum = np.mean(data,axis=1)
    wavesum = np.cumsum(wavesum**2)
    wavesum = wavesum/np.max(wavesum)*100
    pulse_length = np.max(twt[wavesum<97.5])-np.max(twt[wavesum<2.5])
    subgrouphdf5.create_dataset('max_signal_length',data=pulse_length,dtype='f')    

            
    
def WriteAirgunInfo(airguninput,hf,source_name,source_SPL):
    #Check whether there are airguns to be added and stop if nothing is to be added
    if airguninput['airgun_count']==0:
        print('         There are no airguns to be added.')
        
    #Otherwise: Add the group.
    else:
        #Create the airgun group and add the count
        hfair = hf.create_group('airguns')
        hfair.create_dataset('airgun_count',data=airguninput['airgun_count'],dtype='i')
        
        #Make a list of the airgun configurations and step through the list
        airgunlist = list(airguninput.keys())[1:]
        for ind,step in enumerate(airgunlist):
            print('                  Adding configuration: ', step)
            
            #Create the group
            hfairsub = hfair.create_group(step)
            hfairsub.create_dataset('comment',data=airguninput[step]['comment'])
            hfairsub.create_dataset('pulse_rate',data=airguninput[step]['pulse_rate'],dtype='f')
            
            #Open the source config 
            yaml_file = open("../Data/Sources/"+airguninput[step]['ident']+'_Specs.yaml','r')
            yaml_input = yaml.safe_load(yaml_file)
            airgunspecs  = yaml_input[airguninput[step]['ident']]
            
            #Write the specs
            hfairsub.create_dataset('beam_width',data=airgunspecs['beam_width'],dtype='f')
            hfairsub.create_dataset('manufacturer',data=airgunspecs['manufacturer'])
            hfairsub.create_dataset('name',data=airgunspecs['name'])
            hfairsub.create_dataset('serial_number',data=airgunspecs['serial_number'])
            hfairsub.create_dataset('pressure',data=airgunspecs['pressure'],dtype='f')
            hfairsub.create_dataset('source_depth',data=airgunspecs['source_depth'],dtype='f')
            hfairsub.create_dataset('volume_generator',data=airgunspecs['volume_generator'],dtype='f')
            hfairsub.create_dataset('volume_injector',data=airgunspecs['volume_injector'],dtype='f')
            
            #Calculate and write the more sophisticated specs       
            SPL = GetSourceInfo(airguninput,hfairsub,step)
            source_name.append(step)
            source_SPL.append(SPL)
                       
            #Calculate and write the specs based on example profile
            GetSourceExtendedInfo(airguninput,hfairsub,step)
    
    return source_name,source_SPL

def WriteOthersInfo(source,sourceinput,hf,source_name,source_SPL):
    #Check whether there are sources to be added and stop if nothing is to be added
    countname = source+'_count'
    if sourceinput[countname]==0:
        print('         There are no ',source, ' sources to be added.')
        
    #Otherwise: Add the group.
    else:
        #Create the airgun group and add the count
        if source == 'sbp':
            hfsou = hf.create_group('sub_bottom_prifiler')
        else:
            hfsou = hf.create_group(source)
            
        hfsou.create_dataset(countname,data=sourceinput[countname],dtype='i')
        
        #Make a list of the airgun configurations and step through the list
        sourcelist = list(sourceinput.keys())[1:]
        for ind,step in enumerate(sourcelist):
            print('                  Adding configuration: ', step)
            
            #Create the group
            hfsousub = hfsou.create_group(step)
            hfsousub.create_dataset('comment',data=sourceinput[step]['comment'])
            hfsousub.create_dataset('pulse_rate',data=sourceinput[step]['pulse_rate'],dtype='f')
            
            #Open the source config 
            yaml_file = open("../Data/Sources/"+sourceinput[step]['ident']+'_Specs.yaml','r')
            yaml_input = yaml.safe_load(yaml_file)
            sourcespecs  = yaml_input[sourceinput[step]['ident']]
            
            #Write the specs
            hfsousub.create_dataset('beam_width',data=sourcespecs['beam_width'],dtype='f')
            hfsousub.create_dataset('manufacturer',data=sourcespecs['manufacturer'])
            hfsousub.create_dataset('name',data=sourcespecs['name'])
            hfsousub.create_dataset('serial_number',data=sourcespecs['serial_number'])
            hfsousub.create_dataset('electric_energy_max',data=sourcespecs['electric_energy_max'],dtype='f')
            hfsousub.create_dataset('electric_power_max',data=sourcespecs['electric_energy_max']*sourceinput[step]['pulse_rate'],dtype='f')          
            hfsousub.create_dataset('source_depth',data=sourcespecs['source_depth'],dtype='f')
            if source=='boomer':
                hfsousub.create_dataset('plate_area',data=sourcespecs['plate_area'],dtype='f')
                
            #Calculate and write the more sophisticated specs       
            SPL = GetSourceInfo(sourceinput,hfsousub,step)
            source_name.append(step)
            source_SPL.append(SPL)
                       
            #Calculate and write the specs based on example profile
            GetSourceExtendedInfo(sourceinput,hfsousub,step)
    
    return source_name,source_SPL


#%% Transects
def FindStrongestSPL(iterationcounter,which,profiletable,sourcetable,sourceSPL,maxvalue,hdfsubsub):
    filtersource = [col for col in profiletable if col.startswith(which)]
    findmask = np.array(profiletable[filtersource].iloc[iterationcounter]=='x')
    if any(findmask):
        if which == "sbp":
            hdfsubsub.create_dataset('sub_bottom_profiler',data=np.string_(filtersource)[findmask])
        else:
            hdfsubsub.create_dataset(which,data=np.string_(filtersource)[findmask])
        for sourceind,sourcestep in enumerate(np.string_(filtersource)[findmask]):
            sourcemask = (np.string_(sourcetable)==sourcestep)
            if any(sourcemask):
                if (maxvalue<np.array(sourceSPL)[sourcemask]):
                    maxvalue=np.array(sourceSPL)[sourcemask]
    return maxvalue


def AirgunLoudness(SPL):
    l = 'NA'
    if SPL>209:
        l='very_low'
    if SPL>233:
        l='low'
    if SPL>243:
        l='medium'
    if SPL>253:
        l='high'
    return l

def OtherLoudness(SPL):
    l = 'NA'
    if SPL>186:
        l='very_low'
    if SPL>210:
        l='low'
    if SPL>220:
        l='medium'
    if SPL>230:
        l='high'
    return l
   
   

def WriteTransects(hf,navdirectory,profilelist,source_name,source_SPL):
        
    #Create the airgun group and add the count
    hftrans = hf.create_group('transects')
        
    #Import the NMEA record of the cruise
    nav = ImportNMEA(navdirectory)
    
    #Import the EXCEL profilelist
    df = pandas.read_excel(profilelist)
    FORMAT = df.columns
    table = df[FORMAT]

    
    # Starting to fill the transects
    hftrans.create_dataset('transect_count',data=table['Name'].size)
    
    #Initialize ICES Subrectangles record
    ices_airgun_edukt = np.zeros((1,4),dtype=float)
    ices_other_edukt  = np.zeros((1,4),dtype=float)
    
    
    #Loop through transects
    for ind,step in enumerate(table['Name']):
        #General info of every transect
        outname='transect_'+'%03i'%ind
        hftrans_sub = hftrans.create_group(outname)
        hftrans_sub.create_dataset('comment',data=table['Comment'][ind])
        hftrans_sub.create_dataset('name',data=step)
        
        #Positioning and timing of every transect
        mintime=table['StartTime'][ind].to_pydatetime().timestamp()
        maxtime=table['EndTime'][ind].to_pydatetime().timestamp()
        mask = ((nav[:,0]>=mintime)&(nav[:,0]<=maxtime))
        hftrans_sub.create_dataset('datetime',data=nav[mask,0])
        hftrans_sub.create_dataset('location_count',data=np.sum(mask))
        hftrans_sub.create_dataset('location',data=np.round(nav[mask,1:3],decimals=5))
        
        #Create temporary variable for max spl values on profiles
        airgun_max_spl = 0
        other_max_spl  = 0
    
    
        #Write sound sources of every transect
        hftrans_subsub = hftrans_sub.create_group('sound_sources')
            
        airgun_max_spl  = FindStrongestSPL(ind,'airgun',table,source_name,source_SPL,airgun_max_spl,hftrans_subsub)
        other_max_spl   = FindStrongestSPL(ind,'boomer',table,source_name,source_SPL,other_max_spl,hftrans_subsub)
        other_max_spl   = FindStrongestSPL(ind,'sparker',table,source_name,source_SPL,other_max_spl,hftrans_subsub)
        other_max_spl   = FindStrongestSPL(ind,'sbp',table,source_name,source_SPL,other_max_spl,hftrans_subsub)
     
        #print(other_max_spl)
        #safe into an array to safe preliminary info for ices subrectangle record
        temp = np.append(nav[mask,:],airgun_max_spl*np.ones((np.sum(mask),1),dtype=float),axis=1)
        ices_airgun_edukt = np.append(ices_airgun_edukt,temp,axis=0)
        temp = np.append(nav[mask,:],other_max_spl*np.ones((np.sum(mask),1),dtype=float),axis=1)
        ices_other_edukt  = np.append(ices_other_edukt,temp,axis=0)
        del temp
                    
    
    
    ices_airgun_edukt = np.delete(ices_airgun_edukt,0,0)
    ices_other_edukt = np.delete(ices_other_edukt,0,0)
    
    ices_airgun = np.zeros((ices_airgun_edukt.shape[0],3),dtype='|S10')
    ices_other = np.zeros((ices_other_edukt.shape[0],3),dtype='|S10')
    
    #Check file for GIS Import
    gischeck       = open('Airgun_ICES_CoordinatesDate.txt',"w") # Output
    gischeck_o     = open('Other_ICES_CoordinatesDate.txt',"w") # Output
    
    for posind,posstep in enumerate(ices_airgun_edukt[:,0]):
        if ices_airgun_edukt[posind,3] > 0.0:
            ices_airgun[posind,0] = datetime.datetime.strftime(\
                 datetime.datetime.fromtimestamp(posstep),'%Y-%m-%d')
            ices_airgun[posind,1] = IcesSubtRect(\
                 ices_airgun_edukt[posind,1],ices_airgun_edukt[posind,2])
            ices_airgun[posind,2] = AirgunLoudness(\
                 ices_airgun_edukt[posind,3])
            gischeck.write(datetime.datetime.strftime(\
                 datetime.datetime.fromtimestamp(posstep),'%Y-%m-%d %H:%M')+'\t'+\
                 str(ices_airgun_edukt[posind,1])+'\t'+str(ices_airgun_edukt[posind,2])+'\t'+\
                 IcesSubtRect(ices_airgun_edukt[posind,1],ices_airgun_edukt[posind,2])+'\t'+\
                 str(ices_airgun_edukt[posind,3])+'\n')
        if ices_other_edukt[posind,3] >0.0:    
            ices_other[posind,0] = datetime.datetime.strftime(\
                 datetime.datetime.fromtimestamp(ices_other_edukt[posind,0]),'%Y-%m-%d')
            ices_other[posind,1] = IcesSubtRect(\
                 ices_other_edukt[posind,1],ices_other_edukt[posind,2])
            ices_other[posind,2] = OtherLoudness(\
                 ices_other_edukt[posind,3])
            gischeck_o.write(datetime.datetime.strftime(\
                 datetime.datetime.fromtimestamp(posstep),'%Y-%m-%d %H:%M')+'\t'+\
                 str(ices_other_edukt[posind,1])+'\t'+str(ices_other_edukt[posind,2])+'\t'+\
                 IcesSubtRect(ices_other_edukt[posind,1],ices_other_edukt[posind,2])+'\t'+\
                 str(ices_other_edukt[posind,3])+'\n')
        
    gischeck.close()
    gischeck_o.close()
    
    unique_ices_airgun = np.unique(ices_airgun, axis=0)
    unique_ices_other = np.unique(ices_other, axis=0)
    
    hftrans.create_dataset('ICES_SubRectangles_Airgun',data=unique_ices_airgun)
    hftrans.create_dataset('ICES_SubRectangles_Other',data=unique_ices_other)
      
def WriteTransects_Denmark(navdirectory,profilelist,source_name,source_SPL):
        
    #Import the NMEA record of the cruise
    nav = ImportNMEA(navdirectory)
    
    #Import the EXCEL profilelist
    df = pandas.read_excel(profilelist)
    FORMAT = df.columns
    table = df[FORMAT]

    
    #Initialize ICES Subrectangles record
    ices_airgun_edukt = np.zeros((1,4),dtype=float)
    ices_other_edukt  = np.zeros((1,4),dtype=float)
    
    
    #Loop through transects
    for ind,step in enumerate(table['Name']):
        #General info of every transect
        
        #Positioning and timing of every transect
        mintime=table['StartTime'][ind].to_pydatetime().timestamp()
        maxtime=table['EndTime'][ind].to_pydatetime().timestamp()
        mask = ((nav[:,0]>=mintime)&(nav[:,0]<=maxtime))
        
        #Create temporary variable for max spl values on profiles
        airgun_max_spl = 0
        other_max_spl  = 0
    
       
        #print(other_max_spl)
        #safe into an array to safe preliminary info for ices subrectangle record
        temp = np.append(nav[mask,:],airgun_max_spl*np.ones((np.sum(mask),1),dtype=float),axis=1)
        ices_airgun_edukt = np.append(ices_airgun_edukt,temp,axis=0)
        temp = np.append(nav[mask,:],other_max_spl*np.ones((np.sum(mask),1),dtype=float),axis=1)
        ices_other_edukt  = np.append(ices_other_edukt,temp,axis=0)
        del temp
                    
    
    
    ices_airgun_edukt = np.delete(ices_airgun_edukt,0,0)
    ices_other_edukt = np.delete(ices_other_edukt,0,0)
    
    ices_airgun = np.zeros((ices_airgun_edukt.shape[0],3),dtype='|S10')
    ices_other = np.zeros((ices_other_edukt.shape[0],3),dtype='|S10')
    
    #Check file for GIS Import
    gischeck       = open('Airgun_ICES_CoordinatesDate.txt',"w") # Output
    gischeck_o     = open('Other_ICES_CoordinatesDate.txt',"w") # Output
    
    for posind,posstep in enumerate(ices_airgun_edukt[:,0]):
        if ices_airgun_edukt[posind,3] > 0.0:
            ices_airgun[posind,0] = datetime.datetime.strftime(\
                 datetime.datetime.fromtimestamp(posstep),'%Y-%m-%d')
            ices_airgun[posind,1] = IcesSubtRect(\
                 ices_airgun_edukt[posind,1],ices_airgun_edukt[posind,2])
            ices_airgun[posind,2] = AirgunLoudness(\
                 ices_airgun_edukt[posind,3])
            gischeck.write(datetime.datetime.strftime(\
                 datetime.datetime.fromtimestamp(posstep),'%Y-%m-%d %H:%M')+'\t'+\
                 str(ices_airgun_edukt[posind,1])+'\t'+str(ices_airgun_edukt[posind,2])+'\t'+\
                 IcesSubtRect(ices_airgun_edukt[posind,1],ices_airgun_edukt[posind,2])+'\t'+\
                 str(ices_airgun_edukt[posind,3])+'\n')
        if ices_other_edukt[posind,3] >0.0:    
            ices_other[posind,0] = datetime.datetime.strftime(\
                 datetime.datetime.fromtimestamp(ices_other_edukt[posind,0]),'%Y-%m-%d')
            ices_other[posind,1] = IcesSubtRect(\
                 ices_other_edukt[posind,1],ices_other_edukt[posind,2])
            ices_other[posind,2] = OtherLoudness(\
                 ices_other_edukt[posind,3])
            gischeck_o.write(datetime.datetime.strftime(\
                 datetime.datetime.fromtimestamp(posstep),'%Y-%m-%d %H:%M')+'\t'+\
                 str(ices_other_edukt[posind,1])+'\t'+str(ices_other_edukt[posind,2])+'\t'+\
                 IcesSubtRect(ices_other_edukt[posind,1],ices_other_edukt[posind,2])+'\t'+\
                 str(ices_other_edukt[posind,3])+'\n')
        
    gischeck.close()
    gischeck_o.close()
    
    
#%% Wrapper Function
def MarineEarsReport(configfile):
    try:
        print('01)      Import parameters from YAML:')
        generalinput,airguninput,boomerinput,sbpinput,sparkerrinput = \
            ReadGeneralYaml(configfile)
        print('         Parameters imported.')
        print('02)      Open hdf5 for writing:')
        hf = h5py.File(generalinput['out_name'], 'w')
        
        source_name,source_SPL = WriteGeneralInfo(hf,generalinput)
        print('         General info written into hdf5.')
        
        print('03)      Starting to add the airgun specifications:')
        source_name,source_SPL = WriteAirgunInfo(airguninput,hf,source_name,source_SPL)
        print('         Adding airguns finished.')
        
        print('04)      Starting to add the boomer specifications:')
        source_name,source_SPL = WriteOthersInfo('boomer',boomerinput,hf,source_name,source_SPL)
        print('         Adding boomers finished.')
        
        
        print('05)      Starting to add the sparker specifications:')
        source_name,source_SPL = WriteOthersInfo('sparker',sparkerrinput,hf,source_name,source_SPL)
        print('         Adding sparkers finished.')
        
        print('06)      Starting to add the subbottom profiler specifications:')
        source_name,source_SPL = WriteOthersInfo('sbp',sbpinput,hf,source_name,source_SPL)
        print('         Adding sparkers finished.')
        
        print('07)      Work on transects:')
        WriteTransects(hf,generalinput['navdirectory'],generalinput['profilelist'],source_name,source_SPL)
        print('         Transects finished.')
        
        print('08)      Wrap up:')
        hf.close()
        print('         All done.')
    except Exception as e:
        hf.close()
        print(e)
        print('!!!!! An error occured while writing. Closing file. Check input and error messages above.')