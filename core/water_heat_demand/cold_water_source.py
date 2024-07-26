#!/usr/bin/env python3

"""
This module provides objects to represent the source(s) of cold water.
"""

class ColdWaterSource:
    """ An object to represent a source of cold water """

    def __init__(self, cold_water_temps, simulation_time, start_day, time_series_step):
        """ Construct a ColdWaterSource object

        Arguments:
        cold_water_temps -- list of cold water temperatures, in deg C (one entry per hour)
        simulation_time  -- reference to SimulationTime object
        start_day        -- first day of the time series, day of the year, 0 to 365 (single value)
        time_series_step -- timestep of the time series data, in hours
        """
        self.__cold_water_temps = cold_water_temps
        self.__simulation_time  = simulation_time
        self.__start_day = start_day
        self.__time_series_step = time_series_step

    def temperature(self):
        """ Return the cold water temperature for the current timestep """
        return self.__cold_water_temps[self.__simulation_time.time_series_idx(self.__start_day, self.__time_series_step)]
