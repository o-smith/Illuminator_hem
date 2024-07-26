#!/usr/bin/env python3

"""
This module provides objects to model other hot water uses (any draw from the tap)
"""

# Local imports
import core.units as units
from core.material_properties import WATER
from core.water_heat_demand.misc import frac_hot_water


class OtherHotWater:
    """ An object to model all other hot water use """

    def __init__(self, flowrate, cold_water_source):
        """ Construct a OtherHotWater object

        Arguments:
        flowrate        -- tap/outlet flow rate, in litres/minute
        cold_water_feed -- reference to ColdWaterSource object representing the
                           cold water feed attached to the shower
        """
        self.__flowrate          = flowrate
        self.__cold_water_source = cold_water_source

    def get_flowrate(self):
        return self.__flowrate
        
    def get_cold_water_source(self):
        return(self.__cold_water_source)
        
    def hot_water_demand(self, temp_target, temp_hot_water, total_demand_duration):
        """ Calculate volume of hot water required

        (and volume of warm water draining to WWHRS, if applicable)

        Arguments:
        temp_target           -- temperature of warm water delivered at tap/outlet head, in Celcius
        total_demand_duration -- cumulative running time of this event during the current
                                 timestep, in minutes
        """
        temp_cold = self.__cold_water_source.temperature()

        # TODO Account for behavioural variation factor fbeh
        vol_warm_water = self.__flowrate * total_demand_duration
        # ^^^ litres = litres/minute * minutes
        vol_hot_water  = vol_warm_water * frac_hot_water(temp_target, temp_hot_water, temp_cold)

        return vol_hot_water, vol_warm_water

