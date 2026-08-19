# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
ImpactAssessment_Report.py

Functions to generate the report for the underwater noise impact 
assessment caused by seismics.

Author:             Nikolas Römer-Stange
Initial Draft:      2021.11.02
Last Update:        2024.12.04

Dependencies:       Defined in section 01 of the script.

Units:              SI unless stated differently.

Nomenclature:       Defined upon first call of variable

References:         See subrutines/functions
    
Comments:           Noise impact assessment as required by Danish
                    and German authorities.

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

#Custom functions
import ImpactAssessment_Functions as iafunc


#%% 02) Function to generate the report for the seismics

def GenerateReport(yaml_file):
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
      
        
    #Make the Spreading plot
    iafunc.SpreadingPlot(general_input['spreading'],general_input['spreading_dBred'],\
        general_input['spreading_wd'], general_input['spreading_safety'],\
        general_input['kraken_count'],general_input['kraken_model'],\
        general_input['kraken_depth'],general_input['kraken_selection'])
        
    string1 = 'To account for the geometric spreading, the '+\
              general_input['spreading']+'-model is used \
              (see e.g. Elmer et al., 2007; Bradley and Stern, 2008; Duncan and Parsons, 2011). '
    if ((general_input['spreading'] == 'Spherical') or \
        (general_input['spreading'] == 'Cylindrical') or\
        (general_input['spreading'] == 'Elmer')):
        string2 = 'For this model, only the distance to the source is relevant.'
    elif (general_input['spreading'] == 'ShallowWater'):
        string2 = 'For this geometric spreading model, the waterdepth '+\
                  str(general_input['spreading_wd']) + ' m and the safety threshold '+\
                  str(general_input['spreading_safety']) + ' dB are relevant.'
    elif (general_input['spreading'] == 'Kraken'):
        string2 = 'Modelling has been performed for this spreading model.'
                  
    #Make the absorption plot
    iafunc.AbsorptionPlot(general_input['abs_model'])
  
    #Make the weights plot
    iafunc.FilterFunctionPlot(general_input['group'])
    
    #SPL distance plot
    ptslimit_spl, ttslimit_spl = iafunc.SPL_Loop(general_input)
    print(ptslimit_spl, ttslimit_spl)
    
    #SEL distance plot
    ptslimit_sel, ttslimit_sel = iafunc.SEL_Loop(general_input)
    print(ptslimit_spl, ttslimit_spl)
    
    #SPLrms plot
    blimit = iafunc.SPLrms_Loop(general_input)
    
    #fig = plt.gcf()
    #fig.set_size_inches(203.2/24.4, 76.3/25.4,forward=True)
    #fig.set_dpi(300)
    #plt.tight_layout()
    #fig.savefig('test2png.png', dpi=300)
    
    #Make the actual document
    geometry_options = {"hmargin": "1cm","vmargin":"2cm","paper": "a4paper"}
    doc = Document(geometry_options=geometry_options)
    
    # Add document header
    header = PageStyle("header")
    
    now = datetime.datetime.now()
    with header.create(Head("L")):
        header.append("Date: ")
        header.append(LineBreak())
        header.append(now.strftime("%Y-%m-%d"))
        
    with header.create(Head("R")):
        header.append("University of Bremen")
        header.append(LineBreak())
        header.append('Marine Technology / Environmental Research')
        
    with header.create(Foot("C")):
        header.append(simple_page_number())

    doc.preamble.append(header)
    doc.change_document_style("header")

    # Add Heading
    with doc.create(MiniPage(align='c')):
        doc.append(LargeText(bold("Underwater Noise Impact Assessment")))
        doc.append(LineBreak())
        doc.append(LineBreak())
        doc.append(MediumText(bold("Sound source: "+sspecs['manufacturer']+" "+sspecs['name'])))
        doc.append(LineBreak())
        doc.append(MediumText(bold("Animal group: "+general_input['group'])))
        
    #Add Animal Info    
    with doc.create(Section('General information')):
        with doc.create(Subsection('Animal group '+general_input['group'])):
            doc.append(groupinfo)
            doc.append(LineBreak())
            doc.append('Permanent and Temporary Threshold Shift (PTS, TTS) limits are determined \
                        with the more conservative measure of the dual exposure metrics \
                      of unweighted zero-peak Sound Pressure Level as compiled from\
                      Southall et al. (2007, 2019); BOEM (2014a,b); Statoil ASA (2015); \
                      Tougaard et al. (2014, 2015, 2016); Schack et al. (2019): ')
            with doc.create(Alignat(numbering=False, escape=False)) as agn:
                agn.append(r'SPL_{uw,0-peak}(PTS) &=')
                agn.extend([yaml_input['PTS']['SPL'][general_input['group']], '\,dB\,re\,1\,\mu\,Pa'])
            with doc.create(Alignat(numbering=False, escape=False)) as agn:
                agn.append(r'SPL_{uw,0-peak}(TTS) &=')
                agn.extend([yaml_input['TTS']['SPL'][general_input['group']], '\,dB\,re\,1\,\mu\,Pa'])
            doc.append(' and the weighted Sound Exposure Level:')
            with doc.create(Alignat(numbering=False, escape=False)) as agn:
                agn.append(r'SEL_{w}(PTS) &=')
                agn.extend([yaml_input['PTS']['SEL'][general_input['group']], '\,dB\,re\,1\,\mu\,Pa^2 s'])
            with doc.create(Alignat(numbering=False, escape=False)) as agn:
                agn.append(r'SEL_{w}(TTS) &=')
                agn.extend([yaml_input['TTS']['SEL'][general_input['group']], '\,dB\,re\,1\,\mu\,Pa^2 s'])
                
            doc.append('The noise limits for behavioral effects are based on the weighted Root \
                        Mean Square Sound Pressure Level with an averaging time of 125 ms (Tougaard et al., 2014, 2015, 2016):')
            with doc.create(Alignat(numbering=False, escape=False)) as agn:
                agn.append(r'SPL_{w,rms} &=')
                agn.extend([yaml_input['Behaviour']['SPLrms'][general_input['group']], '\,dB\,re\,1\,\mu\,Pa'])
                
            doc.append(NewPage())
      #Add source info  
        with doc.create(Subsection('Sound Source '+sspecs['manufacturer']+" "+sspecs['name'])):
            doc.append('The source wavelet is shown in the upper panel of \
                        Figure 1. In the second the source spectrum is given. \
                        This is contrasted with the hearing threshold of the ' +general_input['group'] + ' in the third panel. '+
                        'The signal is characterized by a frequency content of '+\
                        '%3.1e'%wavelet_fmin+' Hz to '+'%3.1e'%wavelet_fmax+\
                        ' Hz at -6 dB relative to the maximum power, a pulse length of '+\
                        '%3.1f'%(wavelet_pulse_length*1000)+' ms and the following UNweighted \
                        source intensities at a reference distance of 1m:')
                       
            with doc.create(Alignat(numbering=False, escape=False)) as agn:
                agn.append(r'SPL_{uw,0-peak}(1\,m) &=')
                agn.extend([np.round(wavelet_SPL), '\,dB\,re\,1\,\mu\,Pa'])
            with doc.create(Alignat(numbering=False, escape=False)) as agn:
                agn.append(r'SEL_{uw}(1\,m) &=')
                agn.extend([np.round(wavelet_SEL), '\,dB\,re\,1\,\mu\,Pa^2 s'])
                
            with doc.create(Figure(position='h!')) as source_pic:
                source_pic.add_image('Seismics_Source', width='300pt')
                source_pic.add_caption('Source characteristics of the ' + sspecs['manufacturer']+" "+sspecs['name']+'; '\
                    'First panel: Time series of the wavelet record; \
                    Second panel: Source power spectrum; ' \
                    'Third panel: Hearing threshold of the '+general_input['group']+'.')
            doc.append('For a single seismic source, the directivity is estimated \
                        with the Lloyds Mirror effect based on Carey (2009). \
                        The most relevant factors are the dominant frequency of ' +\
                        '%.1f'%sspecs['fd']+' Hz and the source tow depth of ' +\
                        '%.1f'%sspecs['towdepth']+' m. The resulting and applied directivity \
                        are shown in Figure 2.') 
            with doc.create(Figure(position='h!')) as source_pic:
                source_pic.add_image('Seismics_Source_Directivity', width='300pt')
            with doc.create(Figure(position='h!')) as source_pic:
                source_pic.add_image('Seismics_Source_DirectivityApply', width='300pt')
                source_pic.add_caption('Upper panel: Calculated directivity, \
                    Lower panel: Directivity function applied to the exposure metric calculation.')
            doc.append(NewPage())    
                    
    with doc.create(Section('Calculation Basics')):
        doc.append('For the calculation of the metrics, knowledge about the \
                    geometric spreading as shwon in Figure 2, absorption and the filter function \
                    for the functional hearing group are relevant.')
        with doc.create(Subsection('Geometric Spreading')):
            doc.append(string1+string2)
            
            with doc.create(Figure(position='h!')) as source_pic:
                source_pic.add_image('Seismics_Spreading', width='300pt')
                source_pic.add_caption('Geometric spreading approximations; \
                    The used approximation is drawn as a blue solid line.')
                    
        with doc.create(Subsection('Absorption')):
            doc.append('The absorption is calculated according to Ainslie and McColm (1998)\
                        for the '+general_input['abs_model']+\
                        ' Sea and shwon in Figure 3.')
            with doc.create(Figure(position='h!')) as source_pic:
                source_pic.add_image('Seismics_Absorption', width='300pt')
                source_pic.add_caption('Absorption model.')
                
        with doc.create(Subsection('Filter Function')):
            if ((general_input['group']=='Fish') or\
        (general_input['group']=='ST') or \
        (general_input['group']=='Human')):
                doc.append('No frequency dependent filter function is applied.')
            else:
                doc.append('The frequency dependent filter function shwon in Figure 4 is calculated \
                            from the weight function as specified in Southall et al. (2019).')
                with doc.create(Figure(position='h!')) as source_pic:
                    source_pic.add_image('Seismics_Filter', width='300pt')
                    source_pic.add_caption('Frequency dependent filter function \
                                            for the functional hearing group '+ \
                                            general_input['group']+'.')         
    doc.append(NewPage())
    doc.append(NewPage())
                    
    with doc.create(Section('Safety distances')):
        doc.append('The dual exposure metric unweigthed SPL and weighted SEL \
                    are calculated to determine the safety distances. To give \
                    conservative distance estimates, the larger distance of the\
                    dual metric should be considered for both PTS and TTS.')
        with doc.create(Subsection('SPL limits PTS and TTS')):
            doc.append('Based on the source signal, absproption and geometrical spreading, \
                    the PTS onset distance limit is ' + str(np.round(ptslimit_spl)) +' m, \
                    while the TTS onset distance limit is ' + str(np.round(ttslimit_spl)) +' m, \
                    as shown in Figure 5.')
            with doc.create(Figure(position='h!')) as source_pic:
                source_pic.add_image('Seismics_SPLdist', width='300pt')
                source_pic.add_caption('Reduction of the unweighted zero-peak SPL as a function of \
                    distance to the source with an overlay of the noise limits.')
                
        with doc.create(Subsection('SEL limits PTS and TTS')):
            doc.append('Based on the source signal, absproption and geometrical spreading, \
                    the PTS onset distance limit is ' + str(np.round(ptslimit_sel)) +' m, \
                    while the TTS onset distance limit is ' + str(np.round(ttslimit_sel)) +' m.\
                    A stationary animal and varying passing distances, which are given on the x-\
                  axis of Figure 6, are considered. For the calculation of the cummulative SEL\
                  , a profile duration \
                  of '+str(general_input['profile_dur'])+ ' h with a survey speed of '+\
                  str(general_input['profile_spd']) +'kn and a shot repitition rate of '+\
                  str(general_input['profile_sr']) + 's are considered.')
            with doc.create(Figure(position='h!')) as source_pic:
                source_pic.add_image('Seismics_SELdist', width='300pt')
                source_pic.add_caption('Reduction of the weighted and cummulative SEL \
                    as a function of the distance to a profile line assuming a \
                    stationary animal.')
            
        with doc.create(Subsection('RMS-SPL limits Behavioural Impacts')):
            doc.append('Based on the source signal, absproption and geometrical spreading, \
                    the onset distance limit for behavioural effects is ' + str(np.round(blimit)) +' m\
                    as shown in Figure 6.')
            with doc.create(Figure(position='h!')) as source_pic:
                source_pic.add_image('Seismics_SplRmsdist', width='300pt')
                source_pic.add_caption('Reduction of the weighted RMS SPL with an \
                    averaging preiod of 125 ms as a function of \
                    distance to the source with an overlay of the noise limit.')
                
           



    doc.generate_pdf('Seismics_UnderwaterNoiseImpactAssessment_'+general_input['group']+'_'+sspecs['manufacturer']+sspecs['name'], clean_tex=False)
    
    
#%% 03) Function to generate a report for the Single Beam EchoSounder

def SingleBeamReport(yaml_file):
    #Open the parameter file and determine all parameters of interest
    yaml_file_imp = open(yaml_file)
    yaml_input = yaml.safe_load(yaml_file_imp)
    yaml_file_imp.close()
    
    #Import animal information
    animal_file = open("../Data/FilterFunctions/SoundLimits.yaml",'r')
    animal_input = yaml.safe_load(animal_file)
    animal_file.close()
    groupinfo = animal_input['Description'][yaml_input['group']]
    
    pts_spl_d ,tts_spl_d, pts_sel_d, tts_sel_d, b_d, eq_d = \
        iafunc.SingleBeam_Plotting(yaml_file)
        
    fw,_ = iafunc.WeightFunction(yaml_input['group'],yaml_input['f']) #Weight in dB for functional hearing group
    absorp = iafunc.SimpleAbsorption(yaml_input['abs_model'],yaml_input['f']) #Absorption

    #Make the Spreading plot
    iafunc.SpreadingPlot(yaml_input['spreading'],yaml_input['spreading_dBred'],\
        yaml_input['spreading_wd'], yaml_input['spreading_safety'],\
        yaml_input['kraken_count'],yaml_input['kraken_model'],\
        yaml_input['kraken_depth'],yaml_input['kraken_selection'])
        
    string1 = 'To account for the geometric spreading, the '+\
              yaml_input['spreading']+'-model is used \
              (see e.g. Elmer et al., 2007; Bradley and Stern, 2008; Duncan and Parsons, 2011). '
    if ((yaml_input['spreading'] == 'Spherical') or \
        (yaml_input['spreading'] == 'Cylindrical') or\
        (yaml_input['spreading'] == 'Elmer')):
        string2 = 'For this model, only the distance to the source is relevant.'
    elif (yaml_input['spreading'] == 'ShallowWater'):
        string2 = 'For this geometric spreading model, the waterdepth '+\
                  str(yaml_input['spreading_wd']) + ' m and the safety threshold '+\
                  str(yaml_input['spreading_safety']) + ' dB are used.'
  

    #Make the actual document
    geometry_options = {"hmargin": "1cm","vmargin":"2cm","paper": "a4paper"}
    doc = Document(geometry_options=geometry_options)
    
    # Add document header
    header = PageStyle("header")
    
    now = datetime.datetime.now()
    with header.create(Head("L")):
        header.append("Date: ")
        header.append(LineBreak())
        header.append(now.strftime("%Y-%m-%d"))
        
    with header.create(Head("R")):
        header.append("University of Bremen")
        header.append(LineBreak())
        header.append('Marine Technology / Environmental Research')
        
    with header.create(Foot("C")):
        header.append(simple_page_number())

    doc.preamble.append(header)
    doc.change_document_style("header")

    # Add Heading
    with doc.create(MiniPage(align='c')):
        doc.append(LargeText(bold("Underwater Noise Impact Assessment")))
        doc.append(LineBreak())
        doc.append(LineBreak())
        doc.append(MediumText(bold("Sound source: "+yaml_input['name'])))
        doc.append(LineBreak())
        doc.append(MediumText(bold("Animal group: "+yaml_input['group'])))
        
        
    with doc.create(Section('General information')):
        with doc.create(Subsection('Animal group '+yaml_input['group'])):
            doc.append(groupinfo)
            doc.append(LineBreak())
            doc.append('Permanent and Temporary Threshold Shift (PTS, TTS) limits are determined \
                        with the more conservative measure of the dual exposure metrics \
                      of unweighted zero-peak Sound Pressure Level as compiled from\
                      Southall et al. (2007, 2019); BOEM (2014a,b); Statoil ASA (2015); \
                      Tougaard et al. (2014, 2015, 2016); Schack et al. (2019): ')
            with doc.create(Alignat(numbering=False, escape=False)) as agn:
                agn.append(r'SPL_{uw,0-peak}(PTS) &=')
                agn.extend([animal_input['PTS']['SPL'][yaml_input['group']], '\,dB\,re\,1\,\mu\,Pa'])
            with doc.create(Alignat(numbering=False, escape=False)) as agn:
                agn.append(r'SPL_{uw,0-peak}(TTS) &=')
                agn.extend([animal_input['TTS']['SPL'][yaml_input['group']], '\,dB\,re\,1\,\mu\,Pa'])
            doc.append(' and the weighted Sound Exposure Level:')
            with doc.create(Alignat(numbering=False, escape=False)) as agn:
                agn.append(r'SEL_{w}(PTS) &=')
                agn.extend([animal_input['PTS']['SEL'][yaml_input['group']], '\,dB\,re\,1\,\mu\,Pa^2 s'])
            with doc.create(Alignat(numbering=False, escape=False)) as agn:
                agn.append(r'SEL_{w}(TTS) &=')
                agn.extend([animal_input['TTS']['SEL'][yaml_input['group']], '\,dB\,re\,1\,\mu\,Pa^2 s'])
                
            doc.append('The noise limits for behavioral effects are based on the weighted Root \
                        Mean Square Sound Pressure Level with an averaging time of 125 ms (Tougaard et al., 2014, 2015, 2016):')
            with doc.create(Alignat(numbering=False, escape=False)) as agn:
                agn.append(r'SPL_{w,rms} &=')
                agn.extend([animal_input['Behaviour']['SPLrms'][yaml_input['group']], '\,dB\,re\,1\,\mu\,Pa'])
                
            doc.append('According to NMFS (2024) sound exposures below the Effective Quite threshold \
                        do not contribute to PTS or TTS regardless of duration or cumulative exposure. The limit is determined by:')
            with doc.create(Alignat(numbering=False, escape=False)) as agn:
                agn.append(r'SPL_{uw,0-peak} &=')
                agn.extend([animal_input['EffectiveQuite']['SPL'][yaml_input['group']], '\,dB\,re\,1\,\mu\,Pa'])
                    
        with doc.create(Subsection('Sound Source '+yaml_input['name'])):
            doc.append('The sound source is characterized by a frequency of ' +\
                str(yaml_input['f']) + ' Hz, a beam width of '+ \
                str(yaml_input['phi']) + ' deg resulting in the beam pattern shown in Figure 1, is towed at a depth of '+\
                str(yaml_input['towingdepth']) + ' m and emitting the following \
                zero-peak Sound Pressure Level at a reference distance of 1 m:')
            with doc.create(Alignat(numbering=False, escape=False)) as agn:
                agn.append(r'SPL_{0-peak}(1\,m) &=')
                agn.extend([np.round(yaml_input['SPL']), '\,dB\,re\,1\,\mu\,Pa'])
            with doc.create(Figure(position='h!')) as source_pic:
                source_pic.add_image('SingleBeam_Directivity1D', width='300pt')
                source_pic.add_caption('Directivity pattern of the sound source (calculation basics see Lurton, 2016).')    

                    
    with doc.create(Section('Calculation Basics')):
        doc.append('For the calculation of the metrics, knowledge about the \
                    geometric spreading as shwon in Figure 2, absorption and the filter function \
                    for the functional hearing group are relevant. At the frequency used, \
                    the absorption accounts for a loss of ' + '%4.1f'%(-absorp*1000) + \
                    'dB/km. The weight function of the functional hearing group \
                    is quantified at ' + '%4.1f'%fw +' dB. ')
                   
        doc.append(string1+string2)
        with doc.create(Figure(position='h!')) as source_pic:
            source_pic.add_image('Seismics_Spreading', width='300pt')
            source_pic.add_caption('Geometric spreading approximations; \
                The used approximation is drawn as a blue solid line.')
                
    #doc.append(NewPage())
                    
    with doc.create(Section('Safety distances')):
        doc.append('Based on Lurton (2016) and Tougaard (2014, 2015, 2016), \
            the underwater noise metrics are calculated and shwon in Figure 3.\
            The generated underwater sound falls below the Effective Quiet threshold \
            at a distance of '+ '%6.1f'%eq_d +' m. \
            The dual exposure metric unweigthed SPL and weighted SEL \
            are calculated to determine the safety distances for PTS and TTS. To give \
            conservative distance estimates, the larger distance of the\
            dual metric should be considered for both PTS and TTS. The safety distances\
            for PTS are '+ '%6.1f'%pts_spl_d +' m based on SPL and '+ '%6.1f'%pts_sel_d +\
            ' m based on SEL. For TTS the safety distance is ' +'%6.1f'%tts_spl_d +\
            ' m based on SPL and '+ '%6.1f'%tts_sel_d + ' m based on SEL. \
            To avoid behavioural effects, a safety distance of ' + '%6.1f'%b_d +\
            ' m should be kept.')
        
        with doc.create(Figure(position='h!')) as source_pic:
            source_pic.add_image('SingleBeam_NoiseLimits', width='500pt')
            source_pic.add_caption('Underwater noise metrics. \
            The onset limits of negative impacts for marine animals are \
            given as red lines in the plots. In the SPL-plot, the Effective Quite limit is given as a blue line.\
            When no lines are visible, \
            the measures are not exceeded for this device.')
        
        
                    
    doc.generate_pdf('SingleBeam_UnderwaterNoiseImpactAssessment_'+yaml_input['group']+'_'+yaml_input['name'], clean_tex=False)



#%% 04) Function to generate a report for the Multi Beam EchoSounder

def MultiBeamReport(yaml_file):
    #Open the parameter file and determine all parameters of interest
    yaml_file_imp = open(yaml_file)
    yaml_input = yaml.safe_load(yaml_file_imp)
    yaml_file_imp.close()
    
    #Import animal information
    animal_file = open("../Data/FilterFunctions/SoundLimits.yaml",'r')
    animal_input = yaml.safe_load(animal_file)
    animal_file.close()
    groupinfo = animal_input['Description'][yaml_input['group']]
    
    pts_spl_d ,tts_spl_d, pts_sel_d, tts_sel_d, b_d = \
        iafunc.MultiBeam_Plotting(yaml_file)
        
    fw,_ = iafunc.WeightFunction(yaml_input['group'],yaml_input['f']) #Weight in dB for functional hearing group
    absorp = iafunc.SimpleAbsorption(yaml_input['abs_model'],yaml_input['f']) #Absorption

    #Make the Spreading plot
    iafunc.SpreadingPlot(yaml_input['spreading'],yaml_input['spreading_dBred'],\
        yaml_input['spreading_wd'], yaml_input['spreading_safety'],\
        yaml_input['kraken_count'],yaml_input['kraken_model'],\
        yaml_input['kraken_depth'],yaml_input['kraken_selection'])
        
    string1 = 'To account for the geometric spreading, the '+\
              yaml_input['spreading']+'-model is used \
              (see e.g. Elmer et al., 2007; Bradley and Stern, 2008; Duncan and Parsons, 2011). '
    if ((yaml_input['spreading'] == 'Spherical') or \
        (yaml_input['spreading'] == 'Cylindrical') or\
        (yaml_input['spreading'] == 'Elmer')):
        string2 = 'For this model, only the distance to the source is relevant.'
    elif (yaml_input['spreading'] == 'ShallowWater'):
        string2 = 'For this geometric spreading model, the waterdepth '+\
                  str(yaml_input['spreading_wd']) + ' m and the safety threshold '+\
                  str(yaml_input['spreading_safety']) + ' dB are used.'
  

    #Make the actual document
    geometry_options = {"hmargin": "1cm","vmargin":"2cm","paper": "a4paper"}
    doc = Document(geometry_options=geometry_options)
    
    # Add document header
    header = PageStyle("header")
    
    now = datetime.datetime.now()
    with header.create(Head("L")):
        header.append("Date: ")
        header.append(LineBreak())
        header.append(now.strftime("%Y-%m-%d"))
        
    with header.create(Head("R")):
        header.append("University of Bremen")
        header.append(LineBreak())
        header.append('Marine Technology / Environmental Research')
        
    with header.create(Foot("C")):
        header.append(simple_page_number())

    doc.preamble.append(header)
    doc.change_document_style("header")

    # Add Heading
    with doc.create(MiniPage(align='c')):
        doc.append(LargeText(bold("Underwater Noise Impact Assessment")))
        doc.append(LineBreak())
        doc.append(LineBreak())
        doc.append(MediumText(bold("Sound source: "+yaml_input['name'])))
        doc.append(LineBreak())
        doc.append(MediumText(bold("Animal group: "+yaml_input['group'])))
        
        
    with doc.create(Section('General information')):
        with doc.create(Subsection('Animal group '+yaml_input['group'])):
            doc.append(groupinfo)
            doc.append(LineBreak())
            doc.append('Permanent and Temporary Threshold Shift (PTS, TTS) limits are determined \
                        with the more conservative measure of the dual exposure metrics \
                      of unweighted zero-peak Sound Pressure Level as compiled from\
                      Southall et al. (2007, 2019); BOEM (2014a,b); Statoil ASA (2015); \
                      Tougaard et al. (2014, 2015, 2016); Schack et al. (2019): ')
            with doc.create(Alignat(numbering=False, escape=False)) as agn:
                agn.append(r'SPL_{uw,0-peak}(PTS) &=')
                agn.extend([animal_input['PTS']['SPL'][yaml_input['group']], '\,dB\,re\,1\,\mu\,Pa'])
            with doc.create(Alignat(numbering=False, escape=False)) as agn:
                agn.append(r'SPL_{uw,0-peak}(TTS) &=')
                agn.extend([animal_input['TTS']['SPL'][yaml_input['group']], '\,dB\,re\,1\,\mu\,Pa'])
            doc.append(' and the weighted Sound Exposure Level:')
            with doc.create(Alignat(numbering=False, escape=False)) as agn:
                agn.append(r'SEL_{w}(PTS) &=')
                agn.extend([animal_input['PTS']['SEL'][yaml_input['group']], '\,dB\,re\,1\,\mu\,Pa^2 s'])
            with doc.create(Alignat(numbering=False, escape=False)) as agn:
                agn.append(r'SEL_{w}(TTS) &=')
                agn.extend([animal_input['TTS']['SEL'][yaml_input['group']], '\,dB\,re\,1\,\mu\,Pa^2 s'])
                
            doc.append('The noise limits for behavioral effects are based on the weighted Root \
                        Mean Square Sound Pressure Level with an averaging time of 125 ms (Tougaard et al., 2014, 2015, 2016):')
            with doc.create(Alignat(numbering=False, escape=False)) as agn:
                agn.append(r'SPL_{w,rms} &=')
                agn.extend([animal_input['Behaviour']['SPLrms'][yaml_input['group']], '\,dB\,re\,1\,\mu\,Pa'])
                    
        with doc.create(Subsection('Sound Source '+yaml_input['name'])):
            doc.append('The sound source is characterized by a frequency of ' +\
                str(yaml_input['f']) + ' Hz, a beam width of '+ \
                str(yaml_input['phi']) + ' deg resulting in the beam pattern shown in Figure 1, is towed at a depth of '+\
                str(yaml_input['towingdepth']) + ' m and emitting the following \
                zero-peak Sound Pressure Level at a reference distance of 1 m:')
            with doc.create(Alignat(numbering=False, escape=False)) as agn:
                agn.append(r'SPL_{0-peak}(1\,m) &=')
                agn.extend([np.round(yaml_input['SPL']), '\,dB\,re\,1\,\mu\,Pa'])
            with doc.create(Figure(position='h!')) as source_pic:
                source_pic.add_image('MultiBeam_Directivity1D', width='300pt')
                source_pic.add_caption('Directivity pattern of the sound source; The pattern is approximated with the element transducer pattern (across track) and the beam pattern (along track, calculation basics see Lurton, 2016).')    

                    
    with doc.create(Section('Calculation Basics')):
        doc.append('For the calculation of the metrics, knowledge about the \
                    geometric spreading as shwon in Figure 2, absorption and the filter function \
                    for the functional hearing group are relevant. At the frequency used, \
                    the absorption accounts for a loss of ' + '%4.1f'%(-absorp*1000) + \
                    'dB/km. The weight function of the functional hearing group \
                    is quantified at ' + '%4.1f'%fw +' dB. ')
                   
        doc.append(string1+string2)
        with doc.create(Figure(position='h!')) as source_pic:
            source_pic.add_image('Seismics_Spreading', width='300pt')
            source_pic.add_caption('Geometric spreading approximations; \
                The used approximation is drawn as a blue solid line.')
                
    #doc.append(NewPage())
                    
    with doc.create(Section('Safety distances')):
        doc.append('Based on Lurton (2016) and Tougaard (2014, 2015, 2016), \
            the underwater noise metrics are calculated and shwon in Figure 3.\
            The dual exposure metric unweigthed SPL and weighted SEL \
            are calculated to determine the safety distances for PTS and TTS. To give \
            conservative distance estimates, the larger distance of the\
            dual metric should be considered for both PTS and TTS. The safety distances\
            for PTS are '+ '%6.1f'%pts_spl_d +' m based on SPL and '+ '%6.1f'%pts_sel_d +\
            ' m based on SEL. For TTS the safety distance is ' +'%6.1f'%tts_spl_d +\
            ' m based on SPL and '+ '%6.1f'%tts_sel_d + ' m based on SEL. \
            To avoid behavioural effects, a safety distance of ' + '%6.1f'%b_d +\
            ' m should be kept.')
        
        with doc.create(Figure(position='h!')) as source_pic:
            source_pic.add_image('MultiBeam_NoiseLimits', width='500pt')
            source_pic.add_caption('Underwater noise metrics. \
            The onset limits of negative impacts for marine animals are \
            given as red lines in the plots. When no red lines are visible, \
            the measure is not exceeded for this device.')
        
        
                    
    doc.generate_pdf('MultiBeam_UnderwaterNoiseImpactAssessment_'+yaml_input['group']+'_'+yaml_input['name'], clean_tex=False)


    