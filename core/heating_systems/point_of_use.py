#!/usr/bin/env python3

"""
This module provides object(s) to model the behaviour of point of use heaters
"""

import core.water_heat_demand.misc as misc


class PointOfUse:
    """ Class to represent point of use water heaters """

    def __init__(
            self,
            efficiency,
            energy_supply_conn,
            simulation_time,
            cold_feed,
            temp_hot_water
            ):
        """ Construct an InstantElecHeater object

        Arguments:
        efficiency         -- efficiency of the heater, between 0 and 1
        energy_supply_conn -- reference to EnergySupplyConnection object
        simulation_time    -- reference to SimulationTime object
        cold_feed            -- reference to ColdWaterSource object
        """
        self.__efficiency         = efficiency
        self.__energy_supply_conn = energy_supply_conn
        self.__simulation_time    = simulation_time
        self.__cold_feed          = cold_feed
        self.__temp_hot_water     = temp_hot_water

    def get_cold_water_source(self):
        return self.__cold_feed

    def get_temp_hot_water(self):
        return self.__temp_hot_water

    def demand_hot_water(self, volume_demanded_target):
        demand_temp = self.__temp_hot_water
        # TODO set required temperature rather than hard coding - also elsewhere in the code
        if 'temp_hot_water' in volume_demanded_target:
            volume_demanded = volume_demanded_target['temp_hot_water']['warm_vol']
        else:
            volume_demanded = 0.0

        water_energy_demand = misc.water_demand_to_kWh(volume_demanded, demand_temp, self.__cold_feed.temperature())
        
        energy_used = self.demand_energy(water_energy_demand)

        # Assumption is that system specified has sufficient capacity to meet any realistic demand
        return energy_used

    def demand_energy(self, energy_demand):
        """ Demand energy (in kWh) from the heater """
        # Energy that heater is able to supply is limited by power rating
        fuel_demand = energy_demand * (1/self.__efficiency)

        self.__energy_supply_conn.demand_energy(fuel_demand)
        return fuel_demand
