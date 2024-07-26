#!/usr/bin/env python3

"""
This module contains miscellaneous free functions related to water heat demand.
"""

# Local imports
import core.material_properties as material_properties

def frac_hot_water(temp_target, temp_hot, temp_cold):
    """ Calculate the fraction of hot water required when mixing hot and cold
    water to achieve a target temperature

    Arguments:
    temp_target -- temperature to be achieved, in any units
    temp_hot    -- temperature of hot water to be mixed, in same units as temp_target
    temp_cold   -- temperature of cold water to be mixed, in same units as temp_target
    """
    return (temp_target - temp_cold) / (temp_hot - temp_cold)

def water_demand_to_kWh(litres_demand, demand_temp, cold_temp):
        """
        Calculates the kWh energy content of the hot water demand. 
         
        Arguments:
        litres_demand  -- hot water demand in litres
        demand_temp    -- temperature of hot water inside the pipe, in degrees C
        cold_temp      -- temperature outside the pipe, in degrees C
        """
        kWh_demand = (material_properties.WATER.volumetric_energy_content_kWh_per_litre(demand_temp, cold_temp) * litres_demand)

        return(kWh_demand)
