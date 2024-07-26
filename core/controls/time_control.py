#!/usr/bin/env python3

"""
This module provides objects to model time controls.
"""

# Standard library imports
import sys
from heapq import nsmallest
from math import ceil

# Local imports
import core.units as units


class OnOffTimeControl:
    """ An object to model a time-only control with on/off (not modulating) operation """

    def __init__(self, schedule, simulation_time, start_day, time_series_step):
        """ Construct an OnOffTimeControl object

        Arguments:
        schedule         -- list of boolean values where true means "on" (one entry per hour)
        simulation_time  -- reference to SimulationTime object
        start_day        -- first day of the time series, day of the year, 0 to 365 (single value)
        time_series_step -- timestep of the time series data, in hours
        """
        self.__schedule        = schedule
        self.__simulation_time = simulation_time
        self.__start_day = start_day
        self.__time_series_step = time_series_step

    def is_on(self):
        """ Return true if control will allow system to run """
        return self.__schedule[self.__simulation_time.time_series_idx(self.__start_day, self.__time_series_step)]


class ToUChargeControl:
    """ An object to model a control that governs electrical charging of a heat storage device 
        that can respond to signals from the grid, for example when carbon intensity is low """

    def __init__(self, schedule, simulation_time, start_day, time_series_step, charge_level):
        """ Construct a ToUChargeControl object

        Arguments:
        schedule         -- list of boolean values where true means "on" (one entry per hour)
        simulation_time  -- reference to SimulationTime object
        start_day        -- first day of the time series, day of the year, 0 to 365 (single value)
        time_series_step -- timestep of the time series data, in hours
        charge_level     -- Proportion of the charge targeted for each day
        """
        self.__schedule        = schedule
        self.__simulation_time = simulation_time
        self.__start_day = start_day
        self.__time_series_step = time_series_step
        self.__charge_level = charge_level

    def is_on(self):
        """ Return true if control will allow system to run """
        return self.__schedule[self.__simulation_time.time_series_idx(self.__start_day, self.__time_series_step)]

    def target_charge(self):
        """ Return the charge level value from the list given in inputs; one value per day """
        return self.__charge_level[
            self.__simulation_time.time_series_idx_days(self.__start_day, self.__time_series_step)
        ]


class OnOffCostMinimisingTimeControl:

    def __init__(
            self,
            schedule,
            simulation_time,
            start_day,
            time_series_step,
            time_on_daily,
            ):
        """ Construct an OnOffCostMinimisingControl object

        Arguments:
        schedule         -- list of cost values (one entry per time_series_step)
        simulation_time  -- reference to SimulationTime object
        start_day        -- first day of the time series, day of the year, 0 to 365 (single value)
        time_series_step -- timestep of the time series data, in hours
        time_on_daily    -- number of "on" hours to be set per day
        """
        self.__simulation_time = simulation_time
        self.__start_day = start_day
        self.__time_series_step = time_series_step
        self.__time_on_daily = time_on_daily

        timesteps_per_day = int(units.hours_per_day / time_series_step)
        timesteps_on_daily = int(time_on_daily / time_series_step)
        time_series_len_days = ceil(len(schedule) * time_series_step / units.hours_per_day)

        # For each day of schedule, find the specified number of hours with the lowest cost
        self.__schedule = []
        for day in range(0, time_series_len_days):
            # Get part of the schedule for current day
            schedule_day_start = day * timesteps_per_day
            schedule_day_end = schedule_day_start + timesteps_per_day
            schedule_day = schedule[schedule_day_start : schedule_day_end]

            # Find required number of timesteps with lowest costs
            schedule_day_cost_lowest = sorted(set(nsmallest(timesteps_on_daily, schedule_day)))

            # Initialise boolean schedule for day
            schedule_onoff_day = [False] * timesteps_per_day

            # Set lowest cost times to True, then next lowest etc. until required
            # number of timesteps have been set to True
            assert len(schedule_onoff_day) == len(schedule_day)
            timesteps_to_be_allocated = timesteps_on_daily
            for cost in schedule_day_cost_lowest:
                for idx, entry in enumerate(schedule_day):
                    if timesteps_to_be_allocated < 1:
                        break
                    if entry == cost:
                        schedule_onoff_day[idx] = True
                        timesteps_to_be_allocated -= 1

            # Add day of schedule to overall
            self.__schedule.extend(schedule_onoff_day)

    def is_on(self):
        """ Return true if control will allow system to run """
        return self.__schedule[self.__simulation_time.time_series_idx(self.__start_day, self.__time_series_step)]


class SetpointTimeControl:
    """ An object to model a control with a setpoint which varies per timestep """

    def __init__(
            self,
            schedule,
            simulation_time,
            start_day,
            time_series_step,
            setpoint_min=None,
            setpoint_max=None,
            default_to_max=None,
            duration_advanced_start=0.0,
            ):
        """ Construct a SetpointTimeControl object

        Arguments:
        schedule         -- list of float values (one entry per hour)
        simulation_time  -- reference to SimulationTime object
        start_day        -- first day of the time series, day of the year, 0 to 365 (single value)
        time_series_step -- timestep of the time series data, in hours
        setpoint_min -- min setpoint allowed
        setpoint_max -- max setpoint allowed
        default_to_max -- if both min and max limits are set but setpoint isn't,
                          whether to default to min (False) or max (True) 
        duration_advanced_start -- how long before heating period the system
                                   should switch on, in hours
        """
        self.__schedule        = schedule
        self.__simulation_time = simulation_time
        self.__start_day = start_day
        self.__time_series_step = time_series_step
        self.__setpoint_min = setpoint_min
        self.__setpoint_max = setpoint_max
        self.__default_to_max = default_to_max
        self.__timesteps_advstart \
            = round(duration_advanced_start / self.__simulation_time.timestep())

    def in_required_period(self):
        """ Return true if current time is inside specified time for heating/cooling
        
        (not including timesteps where system is only on due to min or max
        setpoint or advanced start)
        """
        schedule_idx = self.__simulation_time.time_series_idx(
            self.__start_day,
            self.__time_series_step,
            )
        setpnt = self.__schedule[schedule_idx]
        return (setpnt is not None)

    def is_on(self):
        """ Return true if control will allow system to run """
        schedule_idx = self.__simulation_time.time_series_idx(
            self.__start_day,
            self.__time_series_step,
            )
        setpnt = self.__schedule[schedule_idx]

        if setpnt is None:
            # Look ahead for duration of warmup period: system is on if setpoint
            # is not None heating period if found
            for timesteps_ahead in range(1, 1 + self.__timesteps_advstart):
                if len(self.__schedule) <= schedule_idx + timesteps_ahead:
                    # Stop looking ahead if we have reached the end of the schedule
                    break
                if self.__schedule[schedule_idx + timesteps_ahead] is not None:
                    # If heating period starts within duration of warmup period
                    # from now, system is on
                    return True

        # For this type of control, system is always on if min or max are set
        if setpnt is None and self.__setpoint_min is None and self.__setpoint_max is None:
            return False
        else:
            return True

    def setpnt(self):
        """ Return setpoint for the current timestep """
        schedule_idx = self.__simulation_time.time_series_idx(
            self.__start_day,
            self.__time_series_step,
            )
        setpnt = self.__schedule[schedule_idx]

        if setpnt is None:
            # Look ahead for duration of warmup period and use setpoint from
            # start of heating period if found
            for timesteps_ahead in range(1, 1 + self.__timesteps_advstart):
                if len(self.__schedule) <= schedule_idx + timesteps_ahead:
                    # Stop looking ahead if we have reached the end of the schedule
                    break
                if self.__schedule[schedule_idx + timesteps_ahead] is not None:
                    # If heating period starts within duration of warmup period
                    # from now, use setpoint from start of heating period
                    setpnt = self.__schedule[schedule_idx + timesteps_ahead]
                    break

        if setpnt is None:
            # If no setpoint value is in the schedule, use the min/max if set
            if self.__setpoint_max is None and self.__setpoint_min is None:
                pass # Use setpnt None
            elif self.__setpoint_max is not None and self.__setpoint_min is None:
                setpnt = self.__setpoint_max
            elif self.__setpoint_min is not None and self.__setpoint_max is None:
                setpnt = self.__setpoint_min
            else: # min and max both set
                if self.__default_to_max is None:
                    sys.exit('ERROR: Setpoint not set but min and max both set, '
                             'and which to use by default not specified')
                elif self.__default_to_max:
                    setpnt = self.__setpoint_max
                else:
                    setpnt = self.__setpoint_min
        else:
            # If there is a maximum limit, take the lower of this and the schedule value
            if self.__setpoint_max is not None:
                setpnt = min(self.__setpoint_max, setpnt)
            # If there is a minimum limit, take the higher of this and the schedule value
            if self.__setpoint_min is not None:
                setpnt = max(self.__setpoint_min, setpnt)
        return setpnt
