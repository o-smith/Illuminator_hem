#!/usr/bin/env python3
import scipy.interpolate

"""
This module provides objects to model waste water heat recovery systems of different types.
"""

class WWHRS_InstantaneousSystemB:
    """ A class to represent instantaneous waste water heat recovery systems with arrangement B
    
    For System B WWHRS, output of the heat exchanger is fed to the shower only
    """

    def __init__(self, flow_rates, efficiencies, cold_water_source, utilisation_factor):
        self.__cold_water_source = cold_water_source
        self.__flow_rates = flow_rates
        self.__efficiencies = efficiencies
        self.__utilisation_factor = utilisation_factor

    def return_temperature(self, temp_target, flowrate_waste_water=None, flowrate_cold_water=None):
        # TODO The cold water flow rate depends on the temperature returned from
        #      this function, which may create a circularity in the calculation.
        #      May need to integrate System B into shower module and/or combine
        #      equations.
        
        # Get cold feed temperature
        temp_cold = self.__cold_water_source.temperature()
        
        # TODO If flowrates have been provided for waste and cold water:
        #    - Calc heat recovered from waste water. Need to do this per shower
        #      individually? Need WWHRS_Connection object?
        wwhrs_efficiency = self.get_efficiency_from_flowrate(flowrate_waste_water) * \
        self.__utilisation_factor
        
        #    - Calc temp of pre-heated water based on heat recovered and flow rates
        temp_cold = temp_cold + ((wwhrs_efficiency/100) * (temp_target - temp_cold))
        # TODO different flow rates on either side of heat excahnger have been accounted for in the PCDF efficiency.
        # This efficiency is based on several assumed parameters and is therefore not robust/flexible to change in those
        # parameters. e.g. waste water temperature, certain flow rates. See CALCM:03 WWHRS tech document for SAP 10.
                
        # Return temp of cold water (pre-heated if relevant)
        return(temp_cold)

    def get_efficiency_from_flowrate(self, flowrate):
        # Get the interpolated efficiency from the flowrate of waste water
        y_interp = scipy.interpolate.interp1d(self.__flow_rates, self.__efficiencies)
        
        return(y_interp(flowrate))
        
class WWHRS_InstantaneousSystemC:
    """ A class to represent instantaneous waste water heat recovery systems with arrangement C
    
    For System C WWHRS, output of the heat exchanger is fed to the hot water system only
    """
    
    def __init__(self, flow_rates, efficiencies, cold_water_source, utilisation_factor):
        self.__cold_water_source = cold_water_source
        self.__stored_temperature = self.__cold_water_source.temperature()
        self.__flow_rates = flow_rates
        self.__efficiencies = efficiencies
        self.__utilisation_factor = utilisation_factor

    def set_temperature_for_return(self, water_temperature):
        self.__stored_temperature = water_temperature

    def temperature(self):
        return (self.__stored_temperature)

    def return_temperature(self, temp_target, flowrate_waste_water=None, flowrate_cold_water=None):
        # TODO The cold water flow rate depends on the temperature returned from
        #      this function, which may create a circularity in the calculation.
        #      May need to integrate System B into shower module and/or combine
        #      equations.
        
        # Get cold feed temperature
        temp_cold = self.__cold_water_source.temperature()
        
        # TODO If flowrates have been provided for waste and cold water:
        #    - Calc heat recovered from waste water. Need to do this per shower
        #      individually? Need WWHRS_Connection object?
        wwhrs_efficiency = self.get_efficiency_from_flowrate(flowrate_waste_water) * \
        self.__utilisation_factor
        
        
        #    - Calc temp of pre-heated water based on heat recovered and flow rates
        temp_cold = temp_cold + ((wwhrs_efficiency/100) * (temp_target - temp_cold))
        # TODO different flow rates on either side of heat excahnger have been accounted for in the PCDF efficiency.
        # This efficiency is based on several assumed parameters and is therefore not robust/flexible to change in those
        # parameters. e.g. waste water temperature, certain flow rates. See CALCM:03 WWHRS tech document for SAP 10.
        
        # Return temp of cold water (pre-heated if relevant)
        return(temp_cold)

    def get_efficiency_from_flowrate(self, flowrate):
        # Get the interpolated efficiency from the flowrate of waste water
        y_interp = scipy.interpolate.interp1d(self.__flow_rates, self.__efficiencies)
        
        return(y_interp(flowrate))

class WWHRS_InstantaneousSystemA:
    """ A class to represent instantaneous waste water heat recovery systems with arrangement A
    
    For System A WWHRS, output of the heat exchanger is fed to both the shower
    and the hot water system
    """
    
    def __init__(self, flow_rates, efficiencies, cold_water_source, utilisation_factor):
        self.__cold_water_source = cold_water_source
        self.__stored_temperature = self.__cold_water_source.temperature()
        self.__flow_rates = flow_rates
        self.__efficiencies = efficiencies
        self.__utilisation_factor = utilisation_factor

    def set_temperature_for_return(self, water_temperature):
        self.__stored_temperature = water_temperature

    def temperature(self):
        return (self.__stored_temperature)

    def return_temperature(self, temp_target, flowrate_waste_water=None, flowrate_cold_water=None):
        
        # Get cold feed temperature
        temp_cold = self.__cold_water_source.temperature()
        
        wwhrs_efficiency = self.get_efficiency_from_flowrate(flowrate_waste_water) * \
        self.__utilisation_factor
        
        #    - Calc temp of pre-heated water based on heat recovered and flow rates
        temp_cold = temp_cold + ((wwhrs_efficiency/100) * (temp_target - temp_cold))
        
        # Return temp of cold water (pre-heated if relevant)
        return(temp_cold)

    def get_efficiency_from_flowrate(self, flowrate):
        # Get the interpolated efficiency from the flowrate of waste water
        y_interp = scipy.interpolate.interp1d(self.__flow_rates, self.__efficiencies)
        
        return(y_interp(flowrate))
