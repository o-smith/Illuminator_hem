#!/usr/bin/env python3

"""
This module contains the hot water demand calculations.
"""

# Standard library imports
import sys

# Local imports
from core.pipework import Pipework_Simple, Location
import core.units as units
import core.water_heat_demand.misc as misc
from core.water_heat_demand.shower import MixerShower, InstantElecShower
from core.water_heat_demand.bath import Bath
from core.water_heat_demand.other_hot_water_uses import OtherHotWater


class DHWDemand:

    # TODO    Enhance analysis for overlapping events
    # Part of draft code for future overlapping analysis of events            
    # For pipework losses count only none overlapping events
    # Time of finalisation of the previous hot water event
    #__time_end_previous_event = 0.0

    def __init__(
            self,
            showers_dict,
            baths_dict,
            other_hw_users_dict,
            hw_pipework_list,
            cold_water_sources,
            wwhrs,
            energy_supplies,
            event_schedules,
            ):
        """ Construct a DHWDemand object """
        self.__event_schedules = event_schedules

        def dict_to_shower(name, data):
            """ Parse dictionary of shower data and return approprate shower object """
            cold_water_source = cold_water_sources[data['ColdWaterSource']]
            # TODO Need to handle error if ColdWaterSource name is invalid.

            shower_type = data['type']
            if shower_type == 'MixerShower':
                wwhrs_instance = None
                if 'WWHRS' in data:
                    wwhrs_instance = wwhrs[data['WWHRS']] # find the instance of WWHRS linked to by the shower

                shower = MixerShower(data['flowrate'], cold_water_source, wwhrs_instance)
            elif shower_type == 'InstantElecShower':
                energy_supply = energy_supplies[data['EnergySupply']]
                # TODO Need to handle error if EnergySupply name is invalid.
                energy_supply_conn = energy_supply.connection(name)

                shower = InstantElecShower(
                    data['rated_power'],
                    cold_water_source,
                    energy_supply_conn,
                    )
            else:
                sys.exit(name + ': shower type (' + shower_type + ') not recognised.')
                # TODO Exit just the current case instead of whole program entirely?
            return shower

        self.__showers = {}
        no_of_showers = 0
        for name, data in showers_dict.items():
            self.__showers[name] = dict_to_shower(name, data)
            # Count number of showers that draw from HW system
            if data['type'] != 'InstantElecShower':
                no_of_showers += 1


        def dict_to_baths(name, data):
            """ Parse dictionary of bath data and return approprate bath object """
            cold_water_source = cold_water_sources[data['ColdWaterSource']]
            # TODO Need to handle error if ColdWaterSource name is invalid.

            bath = Bath(data['size'], cold_water_source, data['flowrate'])

            return bath

        self.__baths = {}
        for name, data in baths_dict.items():
            self.__baths[name] = dict_to_baths(name, data)

        def dict_to_other_water_events(name, data):
            """ Parse dictionary of bath data and return approprate other event object """
            cold_water_source = cold_water_sources[data['ColdWaterSource']]
            # TODO Need to handle error if ColdWaterSource name is invalid.

            other_event = OtherHotWater(data['flowrate'], cold_water_source)

            return other_event

        self.__other_hw_users = {}
        for name, data in other_hw_users_dict.items():
            self.__other_hw_users[name] = dict_to_other_water_events(name, data)

        total_no_of_hot_water_tapping_points = \
            no_of_showers + len(self.__baths.keys()) + len(self.__other_hw_users.keys())

        def dict_to_water_distribution_system(data: dict):
            # Calculate average length of pipework between HW system and tapping point
            length_average = data["length"] / total_no_of_hot_water_tapping_points
            pipework = Pipework_Simple(
                Location.from_string(data["location"]),
                data["internal_diameter_mm"] / units.mm_per_m,
                length_average,
                "water",
            )

            return(pipework)

        self.__hw_distribution_pipework = []

        for data in hw_pipework_list:
            pipework = dict_to_water_distribution_system(data)
            self.__hw_distribution_pipework.append(pipework)

    def hot_water_demand(self, t_idx, temp_hot_water):
        """ Calculate the hot water demand for the current timestep

        Arguments:
        t_idx -- timestep index/count
        """
        hw_demand_vol = 0.0
        hw_demand_vol_target = {}
        hw_energy_demand = 0.0
        hw_duration = 0.0
        all_events = 0.0
        #none_overlapping_events = 0.0
        vol_hot_water_equiv_elec_shower = 0.0

        """
        Events have been organised now so that they are structured by timple step t_idx and 
        sorted for each time step from start to end. 
        """
        # TODO No overlapping is currently considered in terms of interaction with hot water
        #      source (tank). The first event that starts is Served before the second event 
        #      is considered even if this starts before the previous event has finished.
        
        usage_events = self.__event_schedules[t_idx]

        if usage_events is not None:
            for event in usage_events:
                # is_IES = False
                if event['type'] == "Shower":
                    for name, shower in self.__showers.items():
                        if name is event['name']:
                            # If shower is used in the current timestep, get details of use
                            # and calculate HW demand from shower
                            the_cold_water_temp = shower.get_cold_water_source()
                            cold_water_temperature = the_cold_water_temp.temperature()

                            shower_temp = event['temperature']
                            label_temp = str(shower_temp)  # Convert to string for dictionary key
                            shower_duration = event['duration']

                            hw_demand_i, hw_demand_target_i = shower.hot_water_demand(shower_temp, temp_hot_water, shower_duration)

                            if not isinstance(shower, InstantElecShower):
                                event['warm_volume'] = hw_demand_target_i
                                # don't add hw demand and pipework loss from electric shower
                                hw_demand_vol += hw_demand_i
                                hw_energy_demand += misc.water_demand_to_kWh(
                                    hw_demand_i,
                                    temp_hot_water,
                                    cold_water_temperature
                                    )
                                hw_duration += event['duration'] # shower minutes duration
                                all_events += 1
                                
                                if label_temp in hw_demand_vol_target:
                                    hw_demand_vol_target[label_temp]['warm_vol'] += hw_demand_target_i
                                else:
                                    hw_demand_vol_target[label_temp] = {'warm_temp': shower_temp, 'warm_vol': hw_demand_target_i}
                            else:
                                # If electric shower, function returns equivalent
                                # amount of hot water for internal gains calculation
                                vol_hot_water_equiv_elec_shower += hw_demand_i
                                # is_IES = True
                            
                elif event['type'] == "Other":
                    for name, other in self.__other_hw_users.items():
                        if name is event['name']:
                            # If other is used in the current timestep, get details of use
                            # and calculate HW demand from other
                            the_cold_water_temp = other.get_cold_water_source()
                            cold_water_temperature = the_cold_water_temp.temperature()

                            other_temp = event['temperature']
                            label_temp = str(other_temp)  # Convert to string for dictionary key
                            other_duration = event['duration']
                            hw_demand_i, hw_demand_target_i = other.hot_water_demand(other_temp, temp_hot_water, other_duration)
                            event['warm_volume'] = hw_demand_target_i
                            if label_temp in hw_demand_vol_target:
                                hw_demand_vol_target[label_temp]['warm_vol'] += hw_demand_target_i
                            else:
                                hw_demand_vol_target[label_temp] = {'warm_temp': other_temp, 'warm_vol': hw_demand_target_i}
                                
                            hw_demand_vol += hw_demand_i
                            # Check if it makes sense to call again the hot_water_demand function instead of sending hw_demand_i previously calculated
                            # other.hot_water_demand(other_temp, temp_hot_water, other_duration)[0],
                            hw_energy_demand += misc.water_demand_to_kWh(
                                hw_demand_i,
                                temp_hot_water,
                                cold_water_temperature
                                )
                            hw_duration += event['duration'] # other minutes duration
                            all_events += 1
                    
                elif event['type'] == "Bath":
                    for name, bath in self.__baths.items():
                        if name is event['name']:
                            # If bath is used in the current timestep, get details of use
                            # and calculate HW demand from bath
                            the_cold_water_temp = bath.get_cold_water_source()
                            cold_water_temperature = the_cold_water_temp.temperature()

                            # Assume flow rate for bath event is the same as other hot water events
                            peak_flowrate = bath.get_flowrate()
                            bath_temp = event['temperature']
                            label_temp = str(bath_temp)  # Convert to string for dictionary key
                            hw_demand_i, hw_demand_target_i = bath.hot_water_demand(bath_temp, temp_hot_water)
                            event['warm_volume'] = hw_demand_target_i
                            if bath_temp in hw_demand_vol_target:
                                hw_demand_vol_target[label_temp]['warm_vol'] += hw_demand_target_i
                            else:
                                hw_demand_vol_target[label_temp] = {'warm_temp': bath_temp, 'warm_vol': hw_demand_target_i}

                            hw_demand_vol += hw_demand_i
                            # litres bath  / litres per minute flowrate = minutes
                            bath_duration = bath.get_size() / peak_flowrate
                            event['duration'] = bath_duration
                            # Check if it makes sense to call again the hot_water_demand function instead of sending hw_demand_i previously calculated
                            # bath.hot_water_demand(bath_temp, temp_hot_water)[0],
                            hw_energy_demand += misc.water_demand_to_kWh(
                                hw_demand_i,
                                temp_hot_water,
                                cold_water_temperature
                                )
                            hw_duration += bath_duration
                            all_events += 1
                
                # TODO    Enhance analysis for overlapping events
                # Part of draft code for future overlapping analysis of events            
                # For pipework losses count only none overlapping events
                #if not is_IES:
                #    time_start_current_event = event['start']
                #    if time_start_current_event > self.__time_end_previous_event:
                #        none_overlapping_events += 1
                #    # 0.0 can be modified for additional minutes when pipework could be considered still warm/hot    
                #    self.__time_end_previous_event = deepcopy(time_start_current_event + (event['duration'] + 0.0) / 60.0)
                            
        hw_vol_at_tapping_points = hw_demand_vol + vol_hot_water_equiv_elec_shower

        vol_hot_water_left_in_pipework = 0.0
        for pipework in self.__hw_distribution_pipework:
            vol_hot_water_left_in_pipework += pipework.volume_litres()
        hw_demand_vol += all_events * vol_hot_water_left_in_pipework
        # TODO     Refine pipework losses by considering overlapping of events
        #          and shared pipework between serving tap points
        #          none_overlapping_events calculated above is a lower bound(ish)
        #          approximation for this
        if usage_events:
            for event in usage_events:
                event['pipework_volume'] = vol_hot_water_left_in_pipework

        if hw_demand_vol > 0.0:
            hw_demand_vol_target['temp_hot_water'] = {'warm_temp': temp_hot_water, 'warm_vol': hw_demand_vol}

        # Return:
        # - litres hot water per timestep (demand on hw system)
        # - dictionary with litrese of warm water required at different temperature levels
        # - litres hot water per timestep (output at tapping points)
        # - minutes demand per timestep,
        # - number of events in timestep
        # - hot water energy demand (kWh)
        # - usage_events updated to reflect pipework volumes and bath durations 
        return \
            hw_demand_vol, \
            hw_demand_vol_target, \
            hw_vol_at_tapping_points, \
            hw_duration, \
            all_events, \
            hw_energy_demand, \
            usage_events, \
            vol_hot_water_equiv_elec_shower

    def calc_pipework_losses(
            self,
            delta_t_h,
            hw_duration,
            no_of_hw_events,
            demand_water_temperature,
            internal_air_temperature,
            external_air_temperature,
            ):

        if not self.__hw_distribution_pipework:
            # Return heat loss in kWh for the timestep
            return 0.0, 0.0

        cool_down_loss_internal = 0
        cool_down_loss_external = 0

        for pipework in self.__hw_distribution_pipework:
            if pipework.get_location() == Location.INTERNAL:
                cool_down_loss_internal += pipework.cool_down_loss(
                    demand_water_temperature,
                    internal_air_temperature
                    )
            elif pipework.get_location() == Location.EXTERNAL:
                cool_down_loss_external += pipework.cool_down_loss(
                    demand_water_temperature,
                    external_air_temperature
                    )
            else:
                #fallback else block with a failed exception
                raise ValueError(f"Unexpected location value: {pipework._Pipework_Simple__location}")

        pipework_heat_loss_internal \
            = no_of_hw_events * cool_down_loss_internal

        pipework_heat_loss_external \
            = no_of_hw_events * cool_down_loss_external

        # Return heat loss in kWh for the timestep
        return pipework_heat_loss_internal, pipework_heat_loss_external
