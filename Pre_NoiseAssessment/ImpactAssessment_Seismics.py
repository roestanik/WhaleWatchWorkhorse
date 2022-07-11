# -*- coding: utf-8 -*-
"""
ImpactAssessment_Seismics.py

Functions to generate plots and data for the underwater noise impact 
assessment caused by seismics.

Author:             Nikolas Römer-Stange
Initial Draft:      2021.11.02
Last Update:        2020.11.02

Dependencies:       Defined in section 01 of the script.

Units:              SI unless stated differently.

Nomenclature:       Defined upon first call of variable

References:         See subrutines/functions
    
Comments:           Noise impact assessment as required by Danish
                    and German authorities.

"""

#%% Import of dependencies
import ImpactAssessment_SeismicsReport as rep


#%%
yaml_file = "Seismics_Paremeters_SercelMicroGI_VHF.yaml"
rep.GenerateReport(yaml_file)

#%%
yaml_file = "SingleBeam_Parameters_EK60_VHF.yaml"
rep.SingleBeamReport(yaml_file)

#%%
yaml_file = "MultiBeam_Parameters_Reson7125_VHF.yaml"
rep.MultiBeamReport(yaml_file)

