#!/usr/bin/env python3

"""
This module provides objects to model showers of different types.
"""

# Local imports
import core.units as units
from core.material_properties import WATER
from core.water_heat_demand.misc import frac_hot_water, water_demand_to_kWh
import core.heating_systems.wwhrs as wwhrs


class MixerShower:
    """ An object to model mixer showers i.e. those that mix hot and cold water """

    def __init__(self, flowrate, cold_water_source, wwhrs=None):
        """ Construct a MixerShower object

        Arguments:
        flowrate        -- shower's flow rate, in litres/minute
        cold_water_feed -- reference to ColdWaterSource object representing the
                           cold water feed attached to the shower
        """
        self.__flowrate          = flowrate
        self.__cold_water_source = cold_water_source
        self.__wwhrs = wwhrs

    def get_cold_water_source(self):
        return(self.__cold_water_source)
        
    def hot_water_demand(self, temp_target, temp_hot_water, total_shower_duration):
        """ Calculate volume of hot water required

        (and volume of warm water draining to WWHRS, if applicable)

        Arguments:
        temp_target           -- temperature of warm water delivered at shower head, in Celcius
        total_shower_duration -- cumulative running time of this shower during the current
                                 timestep, in minutes
        """
        temp_cold = self.__cold_water_source.temperature()

        # TODO Account for behavioural variation factor fbeh
        vol_warm_water = self.__flowrate * total_shower_duration
        # ^^^ litres = litres/minute * minutes
        vol_hot_water  = vol_warm_water * frac_hot_water(temp_target, temp_hot_water, temp_cold)
        # first calculate the volume of hot water needed if heating from cold water source

        if self.__wwhrs is not None:
            # Assumed temperature entering WWHRS
            temp_drain = 35.0

            # Get the actual return temperature given the temperature and flowrate of the waste water.
            wwhrs_return_temperature = self.__wwhrs.return_temperature(
                temp_drain,
                self.__flowrate,
                None,
                )

            if isinstance(self.__wwhrs, wwhrs.WWHRS_InstantaneousSystemB): # just returns hot water to the shower
                # return the required volume of hot water once the recovered heat has been accounted for.
                vol_hot_water  = vol_warm_water * frac_hot_water(temp_target, temp_hot_water, wwhrs_return_temperature)

            elif isinstance(self.__wwhrs, wwhrs.WWHRS_InstantaneousSystemC): # just returns hot water to the hot water source
                # Set the actual return temperature given the temperature and flowrate of the waste water.
                self.__wwhrs.set_temperature_for_return(wwhrs_return_temperature)
                    
            elif isinstance(self.__wwhrs, wwhrs.WWHRS_InstantaneousSystemA): #  returns hot water to the hot water source and shower
                # Set the actual return temperature given the temperature and flowrate of the waste water.
                self.__wwhrs.set_temperature_for_return(wwhrs_return_temperature)
                
                # return the required volume of hot water once the recovered heat has been accounted for.
                vol_hot_water  = vol_warm_water * frac_hot_water(temp_target, temp_hot_water, wwhrs_return_temperature)

        return vol_hot_water, vol_warm_water


class InstantElecShower:
    """ An object to model instantaneous electric showers

    i.e. those with an electric heating element that heats cold water to the
    desired temperature on-demand
    """

    def __init__(self, rated_power, cold_water_source, elec_supply_conn):
        """ Construct an InstantElecShower object

        Arguments:
        rated_power      -- shower's rated electrical power, in kW
        cold_water_feed  -- reference to ColdWaterSource object representing
                            the cold water feed attached to the shower
        elec_supply_conn -- reference to EnergySupplyConnection object representing
                            the electricity supply attached to the shower
        """
        self.__pwr              = rated_power
        # TODO Does the above account for target temperature? Presumably pwr
        #      stays constant while flow rate changes? Is this how modern
        #      electric showers work, or do they modulate their power output?
        self.__cold_water_source = cold_water_source
        self.__elec_supply_conn  = elec_supply_conn

    def get_cold_water_source(self):
        return(self.__cold_water_source)

    def hot_water_demand(self, temp_target, temp_hot_water, total_shower_duration):
        """ Calculate electrical energy required

        (and volume of warm water draining to WWHRS, if applicable)

        Arguments:
        temp_target           -- temperature of warm water delivered at shower head, in Celcius
        total_shower_duration -- cumulative running time of this shower during
                                 the current timestep, in minutes
        """
        temp_cold = self.__cold_water_source.temperature()

        # TODO Account for behavioural variation factor fbeh
        elec_demand    = self.__pwr * (total_shower_duration / units.minutes_per_hour)
        # ^^^ kWh = kW * hours
        vol_warm_water = elec_demand \
                       / WATER.volumetric_energy_content_kWh_per_litre(temp_target, temp_cold)
        vol_hot_water_equiv \
            = vol_warm_water * frac_hot_water(temp_target, temp_hot_water, temp_cold)

        self.__elec_supply_conn.demand_energy(elec_demand)

        # Instantaneous electric shower heats its own water, so no demand on
        # the water heating system.
        # TODO For now, we return the equivalent volume of hot water so that it
        #      can feed into the internal gains calculation. Ideally, we should
        #      find a way to do this that preserves the same interface between
        #      shower types, but at the time of writing the code that calls this
        #      function handles shower types differently anyway.
        return vol_hot_water_equiv, vol_warm_water
        # TODO Should this return hot water demand or send message to HW system?
        #      The latter would allow for different showers to be connected to
        #      different HW systems, but complicates the implementation of the
        #      HW system as it will have to deal with calls from several
        #      different objects and work out when to amalgamate the figures to
        #      do its own calculation, rather than being given a single overall
        #      figure for each timestep.
        # TODO Also send vol_warm_water to connected WWHRS object? Account for
        #      heat loss between shower head and drain?
