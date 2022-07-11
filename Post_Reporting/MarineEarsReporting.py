# -*- coding: utf-8 -*-
"""
MarineEarsReporting.py

Script to generate the HDF5 data report for the underwater noise registery 
of the BSH/MarineEars

Author:             Nikolas Römer-Stange
Initial Draft:      2021.10.27
Last Update:        2020.11.01

Dependencies:       Defined in section 01 of the script.

Units:              SI unless stated differently.

Nomenclature:       Defined upon first call of variable

References:         See subrutines/functions
    
Comments:           The data report is adjusted to the data requirements for 
                    seismic survey (2021-09-24)

"""

import MarineEarsReportingFuncs as helper
configfile      = "He569_MarineEars.yaml" #General configs to be adopted for each expedition
helper.MarineEarsReport(configfile)
