# -*- coding: utf-8 -*-
"""
ImpactAssessment_Seismics.py

Functions to generate plots and data for the underwater noise impact 
assessment caused by seismics.

Author:             Nikolas Römer-Stange
Initial Draft:      2021.11.02
Last Update:        2025.08.08

Dependencies:       Defined in section 01 of the script.

Units:              SI unless stated differently.

Nomenclature:       Defined upon first call of variable

References:         See subrutines/functions
    
Comments:           Noise impact assessment as required by Danish
                    and German authorities.

"""

#%% Import of dependencies
import ImpactAssessment_SeismicsReport as rep

import matplotlib.pyplot as plt
#%% Seismics
# # Micro-GI reduced to 60 bar
# yaml_file = "Seismics_Paremeters_SercelMicroGI-Red_VHF.yaml"
# rep.GenerateReport(yaml_file)
# yaml_file = "Seismics_Paremeters_SercelMicroGI-Red_PIW.yaml"
# rep.GenerateReport(yaml_file)
# yaml_file = "Seismics_Paremeters_SercelMicroGI-Red_Fish.yaml"
# rep.GenerateReport(yaml_file)
# plt.close('all')

# # Micro-GI
#yaml_file = "Seismics_Paremeters_SercelMicroGI_VHF.yaml"
#rep.GenerateReport(yaml_file)
# yaml_file = "Seismics_Paremeters_SercelMicroGI_PIW.yaml"
# rep.GenerateReport(yaml_file)
# yaml_file = "Seismics_Paremeters_SercelMicroGI_Fish.yaml"
# rep.GenerateReport(yaml_file)
# plt.close('all')

# # # Sparker
# yaml_file = "Seismics_Paremeters_AADura_VHF.yaml"
# rep.GenerateReport(yaml_file)
# yaml_file = "Seismics_Paremeters_AADura_PIW.yaml"
# rep.GenerateReport(yaml_file)
# yaml_file = "Seismics_Paremeters_AADura_Fish.yaml"
# rep.GenerateReport(yaml_file)
# plt.close('all')



#%% Single Beam Echo Sounder 

#SES 2000 compact
# yaml_file = "SingleBeam_Parameters_InnomarSES2000compact_VHF.yaml"
# rep.SingleBeamReport(yaml_file)
# yaml_file = "SingleBeam_Parameters_InnomarSES2000compact_PIW.yaml"
# rep.SingleBeamReport(yaml_file)



# # # SES 2000 medium
yaml_file = "SingleBeam_Parameters_InnomarSES2000compact_VHF_Kathrine.yaml"
rep.SingleBeamReport(yaml_file)
yaml_file = "SingleBeam_Parameters_InnomarSES2000compact_PIW_Kathrine.yaml"
rep.SingleBeamReport(yaml_file)
#yaml_file = "SingleBeam_Parameters_InnomarSES2000_Fish.yaml"
#rep.SingleBeamReport(yaml_file)
#plt.close('all')

# # #EK80
# yaml_file = "SingleBeam_Parameters_EK80-200_VHF.yaml"
# rep.SingleBeamReport(yaml_file)
# yaml_file = "SingleBeam_Parameters_EK80-200_PIW.yaml"
# rep.SingleBeamReport(yaml_file)
# yaml_file = "SingleBeam_Parameters_EK80-200_Fish.yaml"
# rep.SingleBeamReport(yaml_file)
# plt.close('all')
# yaml_file = "SingleBeam_Parameters_EK80-120_VHF.yaml"
# rep.SingleBeamReport(yaml_file)
# yaml_file = "SingleBeam_Parameters_EK80-120_PIW.yaml"
# rep.SingleBeamReport(yaml_file)
# yaml_file = "SingleBeam_Parameters_EK80-120_Fish.yaml"
# rep.SingleBeamReport(yaml_file)
# plt.close('all')
# yaml_file = "SingleBeam_Parameters_EK80-70_VHF.yaml"
# rep.SingleBeamReport(yaml_file)
# yaml_file = "SingleBeam_Parameters_EK80-70_PIW.yaml"
# rep.SingleBeamReport(yaml_file)
# yaml_file = "SingleBeam_Parameters_EK80-70_Fish.yaml"
# rep.SingleBeamReport(yaml_file)
# plt.close('all')
# yaml_file = "SingleBeam_Parameters_EK80-38_VHF.yaml"
# rep.SingleBeamReport(yaml_file)
# yaml_file = "SingleBeam_Parameters_EK80-38_PIW.yaml"
# rep.SingleBeamReport(yaml_file)
# yaml_file = "SingleBeam_Parameters_EK80-38_Fish.yaml"
# rep.SingleBeamReport(yaml_file)
# plt.close('all')

# # # CMax
# yaml_file = "SingleBeam_Parameters_Cmax-100_VHF.yaml"
# rep.SingleBeamReport(yaml_file)
# yaml_file = "SingleBeam_Parameters_Cmax-100_PIW.yaml"
# rep.SingleBeamReport(yaml_file)
# yaml_file = "SingleBeam_Parameters_Cmax-100_Fish.yaml"
# rep.SingleBeamReport(yaml_file)
# plt.close('all')

# yaml_file = "SingleBeam_Parameters_Cmax-325_VHF.yaml"
# rep.SingleBeamReport(yaml_file)
# yaml_file = "SingleBeam_Parameters_Cmax-325_PIW.yaml"
# rep.SingleBeamReport(yaml_file)
# yaml_file = "SingleBeam_Parameters_Cmax-325_Fish.yaml"
# rep.SingleBeamReport(yaml_file)
# plt.close('all')

# #%% Multi Beam Echo Sounder 
# yaml_file = "MultiBeam_Parameters_EM2040-200-red_VHF.yaml"
# rep.MultiBeamReport(yaml_file)
# yaml_file = "MultiBeam_Parameters_EM2040-200-red_PIW.yaml"
# rep.MultiBeamReport(yaml_file)
# yaml_file = "MultiBeam_Parameters_EM2040-200-red_Fish.yaml"
# rep.MultiBeamReport(yaml_file)
# plt.close('all')

# yaml_file = "MultiBeam_Parameters_EM2040-200_VHF.yaml"
# rep.MultiBeamReport(yaml_file)
# yaml_file = "MultiBeam_Parameters_EM2040-200_PIW.yaml"
# rep.MultiBeamReport(yaml_file)
# yaml_file = "MultiBeam_Parameters_EM2040-200_Fish.yaml"
# rep.MultiBeamReport(yaml_file)
# plt.close('all')

# yaml_file = "MultiBeam_Parameters_EM2040-400-red_VHF.yaml"
# rep.MultiBeamReport(yaml_file)
# yaml_file = "MultiBeam_Parameters_EM2040-400-red_PIW.yaml"
# rep.MultiBeamReport(yaml_file)
# yaml_file = "MultiBeam_Parameters_EM2040-400-red_Fish.yaml"
# rep.MultiBeamReport(yaml_file)
# plt.close('all')

# yaml_file = "MultiBeam_Parameters_EM2040-400_VHF.yaml"
# rep.MultiBeamReport(yaml_file)
# yaml_file = "MultiBeam_Parameters_EM2040-400_PIW.yaml"
# rep.MultiBeamReport(yaml_file)
# yaml_file = "MultiBeam_Parameters_EM2040-400_Fish.yaml"
# rep.MultiBeamReport(yaml_file)
# plt.close('all')

