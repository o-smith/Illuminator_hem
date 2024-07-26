#!/usr/bin/env python3

"""
This module provides ways to define schedules which can be expressed concisely
and built from sub-schedules (e.g. construct a weekly schedule from daily
schedules) in input files.
"""

# Standard library imports
import sys
from math import floor

def expand_schedule(sched_type, sched_dict, sched_main, nullable):
    """ Construct a schedule from direct entries or sub-schedules.

    Arguments:
    sched_type -- the type of value contained within the schedule (e.g. bool or
                  float), which cannot be a string or a dict
    sched_dict -- dictionary of schedules (lists) where each schedule element
                  can be either:
                  - a string referencing the name of another schedule in
                    sched_dict
                  - a dict with 'value' and 'repeat' fields, denoting that the
                    value in the 'value' field is repeated the number of times
                    given in the 'repeat field'
                  - a value of the type given in the sched_type argument
    sched_main -- name of main top-level schedule in sched_dict where processing
                  should start
    nullable -- flag denoting whether null values are allowed (True) or not (False)
    """
    if sched_type == dict or sched_type == str:
        # Exit with error if specified value type is dict or string, as these have special meanings
        sys.exit("Schedule type cannot be dict or string")
        # TODO Exit just the current case instead of whole program entirely?

    def process_schedule_entry(sched_entry):
        """ Process a single schedule entry """
        if isinstance(sched_entry, str):
            # If entry is a string, look up sub-schedule with name given by the string
            return process_schedule_entries(sched_dict[sched_entry])
        elif isinstance(sched_entry, dict):
            # If entry is a dict, repeat 'value' element number of times given in 'repeat' element
            # Note that the variable val below is a list
            val = process_schedule_entry(sched_entry['value'])
            return val * sched_entry['repeat']
        elif (isinstance(sched_entry, sched_type) or (isinstance(sched_entry, int) and sched_type is float)) \
        or (nullable and sched_entry is None):
            # If entry is a value of the expected type (e.g. bool or float), store as-is
            # Note: must return a list here, to be consistent with the other returns from this func
            return [sched_entry]
        else:
            # If entry is of an unexpected type, exit with error message
            sys.exit( "Invalid type (" + str(type(sched_entry)) + ") in schedule entry. Expected " \
                    + str(sched_type)
                    )

    def process_schedule_entries(sched):
        """ Process all entries in a schedule (list) """
        sched_expanded = []
        for sched_entry in sched:
            sched_expanded.extend(process_schedule_entry(sched_entry))
        return sched_expanded

    return process_schedule_entries(sched_dict[sched_main])


def expand_events(event_list, sim_timestep, tot_timesteps, name, event_type, schedule):
    """ Construct or update a schedule from a list of events, appending the event type to each event and
        ensuring events are ordered by 'start' time within each timestep.
    
    Arguments:
    event_list        -- list of event dictionaries, where the 'start' element gives
                         the start time of the event, in hours from the start of the simulation
    sim_timestep      -- length of simulation timestep, in hours
    tot_timesteps     -- total number of timesteps in the simulation
    event_type        -- type of the events being processed (e.g., "Shower", "Bath", "Others")
    schedule          -- the existing schedule dictionary to update
    """
    
    for event in event_list:
        starting_timestep = floor(event['start'] / sim_timestep)
        
        if starting_timestep >= 0 and starting_timestep < tot_timesteps:
            event_with_type_name = event.copy()
            event_with_type_name['type'] = event_type
            event_with_type_name['name'] = name

            if schedule[starting_timestep] is None:
                schedule[starting_timestep] = [event_with_type_name]
            else:
                # Insert the event into the correct position to maintain order by 'start' time
                inserted = False
                for i, existing_event in enumerate(schedule[starting_timestep]):
                    if existing_event['start'] > event_with_type_name['start']:
                        schedule[starting_timestep].insert(i, event_with_type_name)
                        inserted = True
                        break
                if not inserted:
                    schedule[starting_timestep].append(event_with_type_name)

    return schedule