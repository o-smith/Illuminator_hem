#!/usr/bin/env python3

"""
This module contains object(s) to track and control information on the
simulation timestep.
"""

# Standard library imports
import math

# Local imports
import core.units as units


class SimulationTime:
    """ An iterator object to track properties relating to the simulation timestep

    This object is a "single source of truth" for information on timesteps, and it controls
    incrementing the timestep. It can be queried by other objects that have references to it.
    """
    # TODO Re-write this class in terms of datetime and timedelta objects
    # TODO Add options to return timestep (timedelta) in particular units
    #      (e.g. seconds, minutes or hours)
    # TODO Account for GMT/BST switchover

    # Define hours that start each month (and end next month). Note there are 13
    # values so that end of final month is handled correctly.
    # E.g. Jan is hours 0-743
    __MONTH_START_END_HOUR = [0, 744, 1416, 2160, 2880, 3624, 4344, 5088, 5832, 6552, 7296, 8016, 8760]

    def __init__(self, starttime, endtime, step):
        """ Construct a SimulationTime object

        Arguments:
        starttime -- The start time of the simulation, in hours from an arbitrary zero point
        endtime   -- The end time of the simulation, in hours from the same
                     arbitrary zero point as starttime
        step      -- The time increment for each step of the calculation, in hours

        Other variables:
        current   -- The current simulation time, in hours from the same
                     arbitrary zero point as starttime
        total     -- Number of timesteps in simulation
        idx       -- Number of timesteps already run (i.e. zero-based ordinal
                     enumeration of current timestep)
        first     -- True if there have been no iterations yet, False otherwise
        """

        self.__step    = step
        self.__end     = endtime
        self.__current = starttime
        self.__starttime = starttime

        self.__total   = math.ceil((endtime - starttime) / step)
        self.__idx     = 0

        self.__first   = True

    def __iter__(self):
        """ Return a reference to this object when an iterator is required """
        return self

    def __next__(self):
        """ Increment simulation timestep """
        if self.__first:
            # If we are on the first iteration, don't increment counters, but
            # set flag to False for next iteration
            self.__first = False
        else:
            # Increment counters
            self.__idx     = self.__idx     + 1
            self.__current = self.__current + self.__step

        if self.__current >= self.__end:
            # If we have reached the end of the simulation, stop iteration
            raise StopIteration

        return self.__idx, self.__current, self.timestep()
    
    def reset(self):
        """ Reset the simulation time to the start time """
        self.__current = self.__starttime
        self.__idx = 0
        self.__first = True

    def current(self):
        """ Return current simulation time """
        return self.__current

    def index(self):
        """ Return ordinal enumeration of current timestep """
        return self.__idx

    def current_hour(self):
        """ Return current hour """
        # Round down to remove fractions of hour
        return int(math.floor(self.__current))

    def hour_of_day(self):
        """ Return hour of day (00:00-01:00 is hour zero) """
        # TODO Assumes that self.__current == 0 is midnight - make this more flexible
        # Remainder from division by hours in day gives time relative to start of day
        time_of_day = self.__current % units.hours_per_day
        # Round down to remove fractions of hour
        return int(math.floor(time_of_day))

    def current_day(self):
        """ Return current day (day 0 is 1st Jan) """
        # Divide current time in hours by hours in day and round down to get current day
        # TODO Assumes that day 0 is 0 <= self.__current < 24 - make this more flexible
        return int(math.floor(self.__current / units.hours_per_day))

    def time_series_idx(self, start_day, time_series_step):
        """ Calculate array lookup index """
        # Index in array of time-series data is current hour (relative to start
        # of year) adjusted for the time-series step (in hours) and for the start day
        # (relative to start of year) of time-series data
        return math.floor((self.current() - start_day * units.hours_per_day) / time_series_step)

    def time_series_idx_days(self, start_day, time_series_step):
        """ Calculate array lookup index """
        # Index in array of time-series data is current day (relative to start
        # of year) adjusted for the time-series step (in days) and for the start day
        # (relative to start of year) of time-series data
        
        # TODO: (Potential) Decide from which hour of the day the system should be targeting next day charge level
        # Currently 9pm
        if self.current_hour() >= 21:
            return math.floor(self.current_day() + 1 - start_day)
        else:
            return math.floor(self.current_day() - start_day)

    def total_steps(self):
        """ Return the total number of timesteps in simulation """
        return self.__total

    def timestep(self):
        """ Return the length of the current timestep, in hours """
        return self.__step

    def current_month(self):
        """ Return current month (0 for January, 11 for December) """
        current_hr = self.current_hour()
        for i, end_hr in enumerate(self.__MONTH_START_END_HOUR):
            # Find first month end which is greater than the current hour, and return
            if current_hr < end_hr:
                return i - 1

    def current_month_start_end_hour(self):
        """ Return the hours upon which the current month starts and ends """
        month_idx = self.current_month()
        return self.__MONTH_START_END_HOUR[month_idx], self.__MONTH_START_END_HOUR[month_idx + 1]
