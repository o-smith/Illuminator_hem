#!/usr/bin/env python3

"""
This module provides objects to model showers of different types.
"""

# Local imports
import core.units as units
from core.material_properties import WATER
from core.water_heat_demand.misc import frac_hot_water


class Bath:
    """ An object to model a bath """

    def __init__(self, size, cold_water_source, flowrate):
        """ Construct a Bath object

        Arguments:
        size            -- bath size in litres - may not be needed but here but retained for flexibility
        cold_water_feed -- reference to ColdWaterSource object representing the
                           cold water feed attached to the shower
        flowrate        -- tap/outlet flow rate, in litres/minute
        """
        self.__bathsize          = size # TODO entire capacity or typical usage?
        self.__cold_water_source = cold_water_source
        self.__flowrate = flowrate

    def get_size(self):
        return self.__bathsize

    def get_cold_water_source(self):
        return(self.__cold_water_source)

    def get_flowrate(self):
        return self.__flowrate

    def hot_water_demand(self, temp_target, temp_hot_water):
        """ Calculate volume of hot water required

        (and volume of warm water draining to WWHRS, if applicable)

        Arguments:
        temp_target           -- temperature of warm water delivered at tap, in Celcius
        """
        temp_cold = self.__cold_water_source.temperature()

        vol_warm_water = self.__bathsize # in litres. we may wish to modify the volume of water
        # compared to the size of the bath

        vol_hot_water  = vol_warm_water * frac_hot_water(temp_target, temp_hot_water, temp_cold)

        return vol_hot_water, vol_warm_water # litres

