#!/usr/bin/env python3

"""
This module provides objects to represent the internal gains.
"""
# Standard library imports
from math import ceil, floor

# Local imports
import core.units as units
from core.schedule import expand_schedule, expand_events


class InternalGains:
    """ An object to represent internal gains """

    def __init__(self, total_internal_gains, simulation_time, start_day, time_series_step):
        """ Construct a InternalGains object

        Arguments:
        total_internal_gains -- list of internal gains, in W/m2 (one entry per hour)
        simulation_time      -- reference to SimulationTime object
        start_day            -- first day of the time series, day of the year, 0 to 365 (single value)
        time_series_step     -- timestep of the time series data, in hours
        """
        self.__total_internal_gains = total_internal_gains
        self.__simulation_time  = simulation_time
        self.__start_day = start_day
        self.__time_series_step = time_series_step

    def total_internal_gain(self, zone_area):
        """ Return the total internal gain for the current timestep in W"""
        return self.__total_internal_gains[self.__simulation_time.time_series_idx(self.__start_day, self.__time_series_step)] * zone_area


class ApplianceGains:
    """ An object to represent internal gains and energy consumption from appliances"""

    def __init__(self, total_energy_supply, energy_supply_conn, gains_fraction, simulation_time, start_day, time_series_step):
        """ Construct a InternalGains object

        Arguments:
        total_energy_supply      -- list of energy supply from appliances, in W / m2 (one entry per hour)
        energy_supply_connection -- reference to EnergySupplyConnection object representing
                                    the electricity supply attached to the appliance
        gains_fraction           -- fraction of energy supply which is counted as an internal gain
        simulation_time          -- reference to SimulationTime object
        start_day                -- first day of the time series, day of the year, 0 to 365 (single value)
        time_series_step         -- timestep of the time series data, in hours
        """
        self.__total_energy_supply = total_energy_supply
        self.__energy_supply_conn = energy_supply_conn
        self.__gains_fraction = gains_fraction
        self.__simulation_time  = simulation_time
        self.__start_day = start_day
        self.__time_series_step = time_series_step

    def total_internal_gain(self, zone_area):
        """ Return the total internal gain for the current timestep, in W """
        # Forward elctricity demand (in kWh) to relevant EnergySupply object
        total_energy_supplied = self.__total_energy_supply[self.__simulation_time.time_series_idx(self.__start_day, self.__time_series_step)] 
        total_energy_supplied_W = total_energy_supplied * zone_area # convert to W
        total_energy_supplied_kWh = total_energy_supplied_W / units.W_per_kW * self.__simulation_time.timestep() # convert to kWh

        self.__energy_supply_conn.demand_energy(total_energy_supplied_kWh)

        return total_energy_supplied_W * self.__gains_fraction


class EventApplianceGains:
    """ An object to represent internal gains and energy consumption from appliances"""

    def __init__(self, energy_supply_conn, simulation_time, appliance_data, TFA):
        """ Construct a InternalGains object

        Arguments:
        energy_supply_connection -- reference to EnergySupplyConnection object representing
                                    the electricity supply attached to the appliance
        simulation_time          -- reference to SimulationTime object
        appliance_data           -- dictionary of appliance gains data from project dict, including:
                                    gains_fraction    -- proportion of appliance demand turned into heat gains
                                    start_day         -- first day of the time series, day of the year, 0 to 365 (single value)
                                    time_series_step  -- timestep of the time series data, in hours
                                    Standby           -- appliance power consumption when not in use in Watts
                                    Events            -- list of appliance usage events, which are dictionaries,
                                                         containing demand_W, start and duration, 
                                                         with start and duration in hours
                                    loadshifting      -- (optional) dictionary definining loadshifting parameters
                                                         max_shift_hrs
                                                         demand_limit_weighted
                                                         weight_timeseries
                                                         demand_timeseries
        TFA                      -- total floor area of dwelling
        """

        self.__energy_supply_conn = energy_supply_conn
        self.__gains_fraction = appliance_data["gains_fraction"]
        self.__simulation_time  = simulation_time
        self.__start_day = appliance_data["start_day"]
        self.__time_series_step = appliance_data["time_series_step"]
        self.__total_floor_area = TFA
        self.__standby_power = appliance_data["Standby"]
        self.__usage_events = appliance_data["Events"]
        if "loadshifting" in appliance_data:
            self.__max_shift = appliance_data["loadshifting"]["max_shift_hrs"] / self.__time_series_step
            self.__demand_limit = appliance_data["loadshifting"]["demand_limit_weighted"]
            self.__weight_timeseries = appliance_data["loadshifting"]["weight_timeseries"]
            self.__otherdemand_timeseries = appliance_data["loadshifting"]["demand_timeseries"]
        else:
            self.__max_shift = 0
            self.__demand_limit = None
            self.__weight_timeseries = None
            self.__otherdemand_timeseries = None
        
        self.__total_power_supply = self.__total_power_supply()

            
        
    def __total_power_supply(self):
        #initialize list with standby power on all timesteps
        total_power_supply = [self.__standby_power for x in range(\
                                ceil(self.__simulation_time.total_steps()
                                     * self.__simulation_time.timestep()
                                     / self.__time_series_step))]
        
        #focus on 2 factors - per appliance shiftability
        #and overall shifting
        
        
        for eventdict in self.__usage_events:
            s, a = self.shift_event(eventdict, self.__otherdemand_timeseries, self.__weight_timeseries)
            for i, x in enumerate(a):
                total_power_supply[floor(s + i) % len(total_power_supply)] += x
        return total_power_supply
    
    def shift_event(self, eventdict, demand_timeseries, weight_timeseries):
        #demand limit could also use ie a linear function instead of a hard limit...
        s, a = self.event_to_schedule(eventdict)
        if self.__max_shift > 0:
            start_shift = self.shift_recursive(s, a, demand_timeseries, weight_timeseries, self.__demand_limit ,self.__max_shift,[], 0)
        else:
            start_shift = 0
        return s + start_shift, a
    
    def shift_recursive(self, s, a, demand_timeseries, weight_timeseries, demandlimit, max_shift, pos_list, start_shift):
        """ 
        shifts an event forward in time one timestep at a time,
        until either the total weighted demand on that timestep is below demandlimit
        or the event has been shifted beyond the maximum allowed number of timesteps
        away from its original position. In the latter case, move the event to the
        most favourable time within the allowed window
        """
        pos_list.append(0)
        for i, x in enumerate(a):
            idx = (floor(s + i) +start_shift) % len(demand_timeseries)
            otherdemand = demand_timeseries[idx] * weight_timeseries[idx]
            newdemand = x * weight_timeseries[idx]
            pos_list[-1] += newdemand + otherdemand
            if newdemand + otherdemand > demandlimit:
                #check if start shift is too high? and if its past limit look up results of
                #each prev shift and choose the best one
                start_shift += 1
                if start_shift <= max_shift:
                    start_shift = self.shift_recursive(s, a, demand_timeseries, weight_timeseries, demandlimit, max_shift, pos_list, start_shift)
                else:
                    #choose the timestep within the allowed window with the lowest demand.
                    #also add the entire length of the series - this will be removed again by the modulo operator,
                    #but prevents an infinite loop from occuring
                    start_shift = pos_list.index(min(pos_list)) + len(demand_timeseries)
        return start_shift
    
    def event_to_schedule(self, eventdict):
        demand_W_event = eventdict["demand_W"]
        start = eventdict["start"]
        duration = eventdict["duration"]
        startoffset = start % self.__time_series_step
        
        #if the event overruns the end of the timestep it starts in,
        #power needs to be allocated to two (or more) timesteps
        #according to the length of time within each timestep the appliance is being used for
        integralx = 0.0
        res = [0 for x in range(ceil(duration / self.__time_series_step))]
        while integralx < duration:
            segment = min(self.__time_series_step - startoffset, duration - integralx)
            idx = floor((integralx) / self.__time_series_step)
            #subtract standby power from the added event power
            #as it is already accounted for when the list is initialised
            res[idx] += (demand_W_event - self.__standby_power) * segment
            integralx += segment
        return start, res

    def total_internal_gain(self, zone_area):
        """ Return the total internal gain for the current timestep, in W """
        # Forward elctricity demand (in kWh) to relevant EnergySupply object
        total_power_supplied = self.__total_power_supply[self.__simulation_time.time_series_idx(self.__start_day, self.__time_series_step)] 
        total_power_supplied_zone = total_power_supplied * zone_area / self.__total_floor_area 
        total_energy_supplied_kWh = total_power_supplied_zone / units.W_per_kW * self.__simulation_time.timestep() # convert to kWh

        self.__energy_supply_conn.demand_energy(total_energy_supplied_kWh)

        return total_power_supplied_zone * self.__gains_fraction
