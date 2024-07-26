#!/usr/bin/env python3

"""
This module provides the high-level control flow for the core calculation, and
initialises the relevant objects in the core model.
"""

# Standard library imports
import sys
from math import ceil, floor

# Local imports
import core.units as units
from core.simulation_time import SimulationTime
from core.external_conditions import ExternalConditions
from core.schedule import expand_schedule, expand_events
from core.ductwork import Ductwork, DuctType
from core.controls.time_control import \
    OnOffTimeControl, SetpointTimeControl, ToUChargeControl, \
    OnOffCostMinimisingTimeControl
from core.cooling_systems.air_conditioning import AirConditioning
from core.energy_supply.energy_supply import EnergySupply
from core.energy_supply.elec_battery import ElectricBattery
from core.energy_supply.pv import PhotovoltaicSystem
from core.heating_systems.emitters import Emitters
from core.heating_systems.heat_pump import HeatPump, HeatPump_HWOnly, SourceType
from core.heating_systems.storage_tank import \
    ImmersionHeater, SolarThermalSystem, StorageTank, PVDiverter
from core.heating_systems.instant_elec_heater import InstantElecHeater
from core.heating_systems.elec_storage_heater import ElecStorageHeater
from core.heating_systems.boiler import Boiler, BoilerServiceWaterCombi
from core.heating_systems.heat_battery import HeatBattery
from core.heating_systems.heat_network import HeatNetwork
from core.space_heat_demand.zone import Zone
from core.space_heat_demand.building_element import \
    BuildingElementOpaque, BuildingElementTransparent, BuildingElementGround, \
    BuildingElementAdjacentZTC, BuildingElementAdjacentZTU_Simple, \
    BuildingElement, HeatFlowDirection
from core.space_heat_demand.ventilation import \
    InfiltrationVentilation, Window, Leaks, \
    MechanicalVentilation, AirTerminalDevices, Vent, CombustionAppliances,\
    air_change_rate_to_flow_rate
from core.space_heat_demand.thermal_bridge import \
    ThermalBridgeLinear, ThermalBridgePoint
from core.water_heat_demand.cold_water_source import ColdWaterSource
from core.water_heat_demand.dhw_demand import DHWDemand
from core.space_heat_demand.internal_gains import InternalGains, ApplianceGains, EventApplianceGains
import core.water_heat_demand.misc as misc
import core.heating_systems.wwhrs as wwhrs
from core.heating_systems.point_of_use import PointOfUse
from core.units import Kelvin2Celcius, seconds_per_hour

# Constants
# Set heating setpoint to absolute zero to ensure no heating demand
temp_setpnt_heat_none = Kelvin2Celcius(0.0)
# Set cooling setpoint to Planck temperature to ensure no cooling demand
temp_setpnt_cool_none = Kelvin2Celcius(1.4e32)

class Project:
    """ An object to represent the overall model to be simulated """

    frac_dhw_energy_internal_gains = 0.25 # used in two functions later for internal gains from pipework

    def __init__(
            self,
            proj_dict,
            print_heat_balance,
            detailed_output_heating_cooling,
            use_fast_solver,
            ):
        """ Construct a Project object and the various components of the simulation

        Arguments:
        proj_dict -- dictionary of project data, containing nested dictionaries
                     and lists of input data for system components, external
                     conditions, occupancy etc.
        print_heat_balance -- flag to idindicate whether to print the heat balance outputs
        detailed_output_heating_cooling -- flag to indicate whether detailed output should be
                                           provided for heating and cooling (where possible)
        use_fast_solver -- flag to indicate whether to use the optimised solver (results
                           may differ slightly due to reordering of floating-point ops)

        Other (self.__) variables:
        simtime            -- SimulationTime object for this Project
        external_conditions -- ExternalConditions object for this Project
        cold_water_sources -- dictionary of ColdWaterSource objects with names as keys
        energy_supplies    -- dictionary of EnergySupply objects with names as keys
        controls           -- dictionary of control objects (of varying types) with names as keys
        hot_water_sources  -- dictionary of hot water source objects (of varying types)
                              with names as keys
        showers            -- dictionary of shower objects (of varying types) with names as keys
        space_heat_systems -- dictionary of space heating system objects (of varying
                              types) with names as keys
        zones              -- dictionary of Zone objects with names as keys
        """
        self.__detailed_output_heating_cooling = detailed_output_heating_cooling

        self.__simtime = SimulationTime(
            proj_dict['SimulationTime']['start'],
            proj_dict['SimulationTime']['end'],
            proj_dict['SimulationTime']['step'],
            )

        # TODO Some inputs are not currently used, so set to None here rather
        #      than requiring them in input file.
        # TODO Read timezone from input file. For now, set timezone to 0 (GMT)
        # Let direct beam conversion input be optional, this will be set if comes from weather file.
        if proj_dict['ExternalConditions']['direct_beam_conversion_needed']:
            dir_beam_conversion = proj_dict['ExternalConditions']['direct_beam_conversion_needed']
        else:
            dir_beam_conversion = False

        def convert_shading(shading_segments):
            """ Function to convert orientation from -180 to +180 (anticlockwise) to 0-360 (clockwise) """
            for element in shading_segments:
                element["start"] = self.__init_orientation(element["start360"])
                element["end"] = self.__init_orientation(element["end360"])
            return shading_segments

        self.__external_conditions = ExternalConditions(
            self.__simtime,
            proj_dict['ExternalConditions']['air_temperatures'],
            proj_dict['ExternalConditions']['wind_speeds'],
            proj_dict['ExternalConditions']['wind_directions'],
            proj_dict['ExternalConditions']['diffuse_horizontal_radiation'],
            proj_dict['ExternalConditions']['direct_beam_radiation'],
            proj_dict['ExternalConditions']['solar_reflectivity_of_ground'],
            proj_dict['ExternalConditions']['latitude'],
            proj_dict['ExternalConditions']['longitude'],
            0, #proj_dict['ExternalConditions']['timezone'],
            0, #proj_dict['ExternalConditions']['start_day'],
            365, #proj_dict['ExternalConditions']['end_day'],
            1, #proj_dict['ExternalConditions']['time_series_step'],
            None, #proj_dict['ExternalConditions']['january_first'],
            None, #proj_dict['ExternalConditions']['daylight_savings'],
            None, #proj_dict['ExternalConditions']['leap_day_included'],
            dir_beam_conversion,
            convert_shading(proj_dict['ExternalConditions']['shading_segments']),
            )

        if 'flat' in proj_dict['General']['build_type']:
            storey_of_dwelling = proj_dict['General']['storey_of_dwelling']
        else:
            storey_of_dwelling = None

        # Loop trough zones to sum up volume to avoid infiltration redundant input.
        total_volume = 0
        for zones, zone_data in proj_dict['Zone'].items():
            total_volume += zone_data['volume']

        self.__cold_water_sources = {}
        for name, data in proj_dict['ColdWaterSource'].items():
            self.__cold_water_sources[name] \
                = ColdWaterSource(data['temperatures'], self.__simtime, data['start_day'], data['time_series_step'])

        self.__energy_supplies = {}
        energy_supply_unmet_demand = EnergySupply('unmet_demand', self.__simtime)
        self.__energy_supplies['_unmet_demand'] = energy_supply_unmet_demand
        diverters = {}
        for name, data in proj_dict['EnergySupply'].items():
            if 'ElectricBattery' in data:
                electric_battery = ElectricBattery(
                    data['ElectricBattery']['capacity'],
                    data['ElectricBattery']['charge_discharge_efficiency_round_trip'],
                    data['ElectricBattery']['battery_age'],
                    data['ElectricBattery']['minimum_charge_rate_one_way_trip'],
                    data['ElectricBattery']['maximum_charge_rate_one_way_trip'],
                    data['ElectricBattery']['maximum_discharge_rate_one_way_trip'],
                    data['ElectricBattery']['battery_location'],
                    self.__simtime,
                    self.__external_conditions
                )
            else:
                electric_battery = None
            self.__energy_supplies[name] = EnergySupply(
                data['fuel'],
                self.__simtime,
                electric_battery,
                data.get('priority'),
                data['is_export_capable']
                )
            # TODO Consider replacing fuel type string with fuel type object

            if 'diverter' in data:
                diverters[name] = data['diverter']

        def dict_to_ctrl(name, data):
            """ Parse dictionary of control data and return appropriate control object """
            ctrl_type = data['type']
            if ctrl_type == 'OnOffTimeControl':
                sched = expand_schedule(bool, data['schedule'], "main", False)
                ctrl = OnOffTimeControl(
                    schedule=sched,
                    simulation_time=self.__simtime,
                    start_day=data['start_day'],
                    time_series_step=data['time_series_step']
                )
            elif ctrl_type == 'SetpointTimeControl':
                sched = expand_schedule(float, data['schedule'], "main", True)

                setpoint_min = None
                setpoint_max = None
                default_to_max = None
                advanced_start = 0.0
                if 'setpoint_min' in data:
                    setpoint_min = data['setpoint_min']
                if 'setpoint_max' in data:
                    setpoint_max = data['setpoint_max']
                if 'default_to_max' in data:
                    default_to_max = data['default_to_max']
                if 'advanced_start' in data:
                    advanced_start = data['advanced_start']

                ctrl = SetpointTimeControl(
                    schedule=sched,
                    simulation_time=self.__simtime,
                    start_day=data['start_day'],
                    time_series_step=data['time_series_step'],
                    setpoint_min=setpoint_min,
                    setpoint_max=setpoint_max,
                    default_to_max=default_to_max,
                    duration_advanced_start=advanced_start,
                )
            elif ctrl_type == 'ToUChargeControl':
                sched = expand_schedule(bool, data['schedule'], "main", False)

                # Simulating manual charge control
                # Set charge_level to 1.0 (max) for each day of simulation (plus 1)
                charge_level = [1.0] * ceil((self.__simtime.total_steps() * self.__simtime.timestep())/24 + 1)
                # If charge_level is present in the input file overwrite initial vector
                # User can specify a vector with all days (plus 1), or as a single float value to be used for each day
                if 'charge_level' in data:
                    # If the input is a vector, use the vector
                    if isinstance(data['charge_level'], (list, tuple)):
                        charge_level=data['charge_level']
                    # Else, if input is a single value, use that value for each day of simulation
                    else:
                        charge_level = [data['charge_level']] * ceil((self.__simtime.total_steps() * self.__simtime.timestep())/24 + 1)

                ctrl = ToUChargeControl(
                    schedule=sched,
                    simulation_time=self.__simtime,
                    start_day=data['start_day'],
                    time_series_step=data['time_series_step'],
                    charge_level=charge_level
                )
            elif ctrl_type == 'OnOffCostMinimisingTimeControl':
                sched = expand_schedule(float, data['schedule'], "main", False)
                ctrl = OnOffCostMinimisingTimeControl(
                    sched,
                    self.__simtime,
                    data['start_day'],
                    data['time_series_step'],
                    data['time_on_daily'],
                    )
            else:
                sys.exit(name + ': control type (' + ctrl_type + ') not recognised.')
                # TODO Exit just the current case instead of whole program entirely?
            return ctrl

        self.__controls = {}
        for name, data in proj_dict['Control'].items():
            self.__controls[name] = dict_to_ctrl(name, data)

        def dict_to_wwhrs(name, data):
            """ Parse dictionary of WWHRS source data and return approprate WWHRS source object """
            wwhrs_source_type = data['type']
            if wwhrs_source_type == 'WWHRS_InstantaneousSystemB':
                cold_water_source = self.__cold_water_sources[data['ColdWaterSource']]
                # TODO Need to handle error if ColdWaterSource name is invalid.

                the_wwhrs = wwhrs.WWHRS_InstantaneousSystemB(
                    data['flow_rates'],
                    data['efficiencies'],
                    cold_water_source,
                    data['utilisation_factor']
                    )
            else:
                if wwhrs_source_type == 'WWHRS_InstantaneousSystemC':
                    cold_water_source = self.__cold_water_sources[data['ColdWaterSource']]
                    # TODO Need to handle error if ColdWaterSource name is invalid.
    
                    the_wwhrs = wwhrs.WWHRS_InstantaneousSystemC(
                        data['flow_rates'],
                        data['efficiencies'],
                        cold_water_source,
                        data['utilisation_factor']
                        )
                else:
                    if wwhrs_source_type == 'WWHRS_InstantaneousSystemA':
                        cold_water_source = self.__cold_water_sources[data['ColdWaterSource']]
                        # TODO Need to handle error if ColdWaterSource name is invalid.
        
                        the_wwhrs = wwhrs.WWHRS_InstantaneousSystemA(
                            data['flow_rates'],
                            data['efficiencies'],
                            cold_water_source,
                            data['utilisation_factor']
                            )
                    else:
                        sys.exit(name + ': WWHRS (' + wwhrs_source_type + ') not recognised.')
                        # TODO Exit just the current case instead of whole program entirely?
            return the_wwhrs
            
        if 'WWHRS' in proj_dict:
            self.__wwhrs = {}
            for name, data in proj_dict['WWHRS'].items():
                self.__wwhrs[name] = dict_to_wwhrs(name, data)
        else:
            self.__wwhrs = None

        def dict_to_event_schedules(data, name, event_type, existing_schedule):
            """ Process list of events (for hot water draw-offs, appliance use etc.) """
            sim_timestep = self.__simtime.timestep()
            tot_timesteps = self.__simtime.total_steps()
            
            # Initialise schedule of events with no events taking place if no existing schedule is provided
            if existing_schedule is None:
                schedule = {t_idx: None for t_idx in range(tot_timesteps)}
            else:
                schedule = existing_schedule
            
            return expand_events(data, sim_timestep, tot_timesteps, name, event_type, schedule)

        self.__event_schedules = None
        for sched_type, schedules in proj_dict['Events'].items():
            for name, data in schedules.items():
                self.__event_schedules = dict_to_event_schedules(data, name, sched_type, self.__event_schedules)

        # TODO - this assumes there is only one hot water source, and if any
        # hot water source is point of use, they all are. In future, allow more
        # than one hot water source and assign hot water source to each outlet?
        if proj_dict['HotWaterSource']['hw cylinder']['type'] == 'PointOfUse':
            hw_pipework_dict = {}
        else:
            hw_pipework_dict = proj_dict['HotWaterDemand']['Distribution']

        if 'Bath' in proj_dict['HotWaterDemand']:
            baths = proj_dict['HotWaterDemand']['Bath']
        else:
            baths = {}

        if 'Shower' in proj_dict['HotWaterDemand']:
            showers = proj_dict['HotWaterDemand']['Shower']
        else:
            showers = {}

        if 'Other' in proj_dict['HotWaterDemand']:
            others = proj_dict['HotWaterDemand']['Other']
        else:
            others = {}

        self.__dhw_demand = DHWDemand(
            showers,
            baths,
            others,
            hw_pipework_dict,
            self.__cold_water_sources,
            self.__wwhrs,
            self.__energy_supplies,
            self.__event_schedules,
            )

        def dict_to_building_element(name, data):
            building_element_type = data['type']

            # Calculate r_c from u_value if only the latter has been provided
            data['r_c'] = self.__init_resistance_or_uvalue(name, data)

            if building_element_type == 'BuildingElementOpaque':

                if data['pitch'] < BuildingElement._BuildingElement__PITCH_LIMIT_HORIZ_CEILING:
                    is_unheated_pitched_roof = data['is_unheated_pitched_roof']
                else:
                    is_unheated_pitched_roof = False

                building_element = BuildingElementOpaque(
                    data['area'],
                    is_unheated_pitched_roof,
                    data['pitch'],
                    data['a_sol'],
                    data['r_c'],
                    data['k_m'],
                    data['mass_distribution_class'],
                    self.__init_orientation(data['orientation360']),
                    data['base_height'],
                    data['height'],
                    data['width'],
                    self.__external_conditions,
                    )
            elif building_element_type == 'BuildingElementTransparent':
                building_element = BuildingElementTransparent(
                    data['pitch'],
                    data['r_c'],
                    self.__init_orientation(data['orientation360']),
                    data['g_value'],
                    data['frame_area_fraction'],
                    data['base_height'],
                    data['height'],
                    data['width'],
                    data['shading'],
                    self.__external_conditions,
                    )
            elif building_element_type == 'BuildingElementGround':
                if data['floor_type'] == 'Slab_no_edge_insulation':
                    edge_insulation = None
                    height_upper_surface = None
                    thermal_transm_envi_base = None
                    thermal_transm_walls = None
                    area_per_perimeter_vent = None
                    shield_fact_location = None
                    thickness_walls = data['thickness_walls']
                    thermal_resist_insul = None
                    depth_basement_floor = None
                    thermal_resist_walls_base = None
                    height_basement_walls = None

                elif data['floor_type'] == 'Slab_edge_insulation':
                    edge_insulation = data['edge_insulation']
                    height_upper_surface = None
                    thermal_transm_envi_base = None
                    thermal_transm_walls = None
                    area_per_perimeter_vent = None
                    shield_fact_location = None
                    thickness_walls = data['thickness_walls']
                    thermal_resist_insul = None
                    depth_basement_floor = None
                    thermal_resist_walls_base = None
                    height_basement_walls = None

                elif data['floor_type'] == 'Suspended_floor':
                    edge_insulation = None
                    height_upper_surface = data['height_upper_surface']
                    thermal_transm_envi_base = None
                    thermal_transm_walls = data['thermal_transm_walls']
                    area_per_perimeter_vent = data['area_per_perimeter_vent']
                    shield_fact_location = data['shield_fact_location']
                    thickness_walls = data['thickness_walls']
                    thermal_resist_insul = data['thermal_resist_insul']
                    depth_basement_floor = None
                    thermal_resist_walls_base = None
                    height_basement_walls = None

                elif data['floor_type'] == 'Heated_basement':
                    edge_insulation = None
                    height_upper_surface = None
                    thermal_transm_envi_base = None
                    thermal_transm_walls = None
                    area_per_perimeter_vent = None
                    shield_fact_location = None
                    thickness_walls = data['thickness_walls']
                    thermal_resist_insul = None
                    depth_basement_floor = data['depth_basement_floor']
                    thermal_resist_walls_base = data['thermal_resist_walls_base']
                    height_basement_walls = None

                elif data['floor_type'] == 'Unheated_basement':
                    edge_insulation = None
                    height_upper_surface = None
                    thermal_transm_envi_base = data['thermal_transm_envi_base']
                    thermal_transm_walls = data['thermal_transm_walls']
                    area_per_perimeter_vent = None
                    shield_fact_location = None
                    thickness_walls = data['thickness_walls']
                    thermal_resist_insul = None
                    depth_basement_floor = data['depth_basement_floor']
                    thermal_resist_walls_base = data['thermal_resist_walls_base']
                    height_basement_walls = data['height_basement_walls']

                else:
                    sys.exit("Type of Floor ("+str(data['floor_type'])+") is not valid")

                building_element = BuildingElementGround(
                    data['total_area'],
                    data['area'],
                    data['pitch'],
                    data['u_value'],
                    data['r_f'],
                    data['k_m'],
                    data['mass_distribution_class'],
                    data['floor_type'],
                    edge_insulation,
                    height_upper_surface,
                    thermal_transm_envi_base,
                    thermal_transm_walls,
                    area_per_perimeter_vent ,
                    shield_fact_location,
                    thickness_walls,
                    thermal_resist_insul,
                    depth_basement_floor,
                    thermal_resist_walls_base,
                    height_basement_walls,
                    data['perimeter'],
                    data['psi_wall_floor_junc'],
                    self.__external_conditions,
                    self.__simtime,
                    )
            elif building_element_type == 'BuildingElementAdjacentZTC':
                building_element = BuildingElementAdjacentZTC(
                    data['area'],
                    data['pitch'],
                    data['r_c'],
                    data['k_m'],
                    data['mass_distribution_class'],
                    self.__external_conditions,
                    )
            elif building_element_type == 'BuildingElementAdjacentZTU_Simple':
                building_element = BuildingElementAdjacentZTU_Simple(
                    data['area'],
                    data['pitch'],
                    data['r_c'],
                    data['r_u'],
                    data['k_m'],
                    data['mass_distribution_class'],
                    self.__external_conditions,
                    )
            else:
                sys.exit( name + ': building element type ('
                        + building_element_type + ') not recognised.' )
                # TODO Exit just the current case instead of whole program entirely?
            return building_element

        windows_dict = {}
        for zone in proj_dict['Zone'].values():
            for building_element_name, building_element in zone['BuildingElement'].items():
                if building_element['type'] == 'BuildingElementTransparent':
                    #Check control for window open exists
                    if building_element.get('Control_WindowOpenable') is None:
                        on_off_ctrl_obj = None
                    else:
                        on_off_ctrl_obj = self.__controls[building_element['Control_WindowOpenable']]
                    windows_dict[building_element_name] = Window(
                        self.__external_conditions,
                        building_element['free_area_height'],
                        building_element['mid_height'],
                        building_element['max_window_open_area'],
                        building_element['window_part_list'],
                        building_element['orientation360'],
                        building_element['pitch'],
                        proj_dict['InfiltrationVentilation']['altitude'],
                        on_off_ctrl_obj,
                    )

        f_cross = proj_dict['InfiltrationVentilation']['cross_vent_factor']
        shield_class = proj_dict['InfiltrationVentilation']['shield_class']
        terrain_class = proj_dict['InfiltrationVentilation']['terrain_class']
        altitude = proj_dict['InfiltrationVentilation']['altitude']
        if proj_dict['InfiltrationVentilation'].get('Control_WindowAdjust') is not None:
            self.__control_WindowAdjust = self.__controls[proj_dict['InfiltrationVentilation']['Control_WindowAdjust']]
        else:
            self.__control_WindowAdjust = None
        
        vents_dict = {}
        for vent_name, vent_data in proj_dict['InfiltrationVentilation']['Vents'].items():
            vents_dict[vent_name] = Vent(
                self.__external_conditions,
                vent_data['mid_height_air_flow_path'],
                vent_data['area_cm2'],
                vent_data['pressure_difference_ref'],
                vent_data['orientation360'],
                vent_data['pitch'],
                proj_dict['InfiltrationVentilation']['altitude'],
            )

        pitches = []
        areas = []
        for zone in proj_dict['Zone'].values():
                for building_element_name, building_element in zone['BuildingElement'].items():
                    if building_element['type'] == 'BuildingElementOpaque':
                        if BuildingElement.pitch_class(building_element['pitch']) == \
                            HeatFlowDirection.UPWARDS:
                            pitches.append(building_element['pitch'])
                            areas.append(building_element['area'])
        #Work out the average pitch, weighted by area.
        area_tot = sum(areas)
        if len(pitches) > 0:
            weighting = [x/area_tot for x in areas]
            weighted_pitches = list(map(lambda x,y: x * y, weighting, pitches))
            average_pitch = sum(weighted_pitches)
        else:
            #This case doesn't matter as if the area of roof = 0, the leakage coefficient = 0 anyway.
            average_pitch = 0

        surface_area_facades_list = []
        surface_area_roof_list = []
        for zone in proj_dict['Zone'].values():
            for building_element_name, building_element in zone['BuildingElement'].items():
                #If wall
                if BuildingElement.pitch_class(building_element['pitch']) == \
                    HeatFlowDirection.HORIZONTAL:
                    if building_element['type'] == 'BuildingElementOpaque':
                        surface_area_facades_list.append(building_element['area'])
                    elif building_element['type'] == 'BuildingElementTransparent':
                        area = building_element['height'] * building_element['width']
                        surface_area_facades_list.append(area)
                #If roof:
                if BuildingElement.pitch_class(building_element['pitch']) == \
                    HeatFlowDirection.UPWARDS:
                    if building_element['type'] == 'BuildingElementOpaque':
                        surface_area_roof_list.append(building_element['area'])
                    elif building_element['type'] == 'BuildingElementTransparent':
                        area = building_element['height'] * building_element['width']
                        surface_area_roof_list.append(area)

        surface_area_facades = sum(surface_area_facades_list)
        surface_area_roof = sum(surface_area_roof_list)

        leaks_dict = proj_dict['InfiltrationVentilation']['Leaks']
        leaks_dict['area_facades'] = surface_area_facades
        leaks_dict['area_roof'] = surface_area_roof
        leaks_dict['altitude'] = proj_dict['InfiltrationVentilation']['altitude']

        # Empty dictionary for air terminal devices until passive ducts work
        atds_dict = {}
        # if proj_dict['InfiltrationVentilation'].get('AirTerminalDevices') is not None:
        #     for atd_name, atds_data in proj_dict['InfiltrationVentilation']['AirTerminalDevices'].items():
        #         atds_dict[atd_name] = AirTerminalDevices(
        #         atds_data['area_cm2'],
        #         atds_data['pressure_difference_ref'],
        #         )

        self.__mech_vents_dict = {}
        self.__space_heating_ductwork_dict = {}
        self.__space_heating_ductwork_dict = None
        # Check if Mechanical Ventilation exists
        if proj_dict['InfiltrationVentilation'].get('MechanicalVentilation') is not None:
            for mech_vents_name, mech_vents_data in proj_dict['InfiltrationVentilation']['MechanicalVentilation'].items():
                # Assign the appropriate control object
                if mech_vents_data.get('Control') is not None:
                    ctrl_intermittent_MEV = self.__controls[mech_vents_data['Control']]
                else:
                    ctrl_intermittent_MEV = None
    
                energy_supply = self.__energy_supplies[mech_vents_data['EnergySupply']]
                # TODO Need to handle error if EnergySupply name is invalid.
                energy_supply_conn = energy_supply.connection(mech_vents_name)

                if mech_vents_data['vent_type'] == 'MVHR':
                    self.__mech_vents_dict[mech_vents_name] = MechanicalVentilation(
                        self.__external_conditions,
                        mech_vents_data['sup_air_flw_ctrl'],
                        mech_vents_data['sup_air_temp_ctrl'],
                        0, #mech_vents_data['design_zone_cooling_covered_by_mech_vent'],
                        0, #mech_vents_data['design_zone_heating_covered_by_mech_vent'],
                        mech_vents_data['vent_type'],
                        mech_vents_data['SFP'],
                        mech_vents_data['design_outdoor_air_flow_rate'],
                        self.__simtime,
                        energy_supply_conn,
                        total_volume,
                        proj_dict['InfiltrationVentilation']['altitude'],
                        ctrl_intermittent_MEV = ctrl_intermittent_MEV,
                        mvhr_eff = mech_vents_data['mvhr_eff'],
                    )
                elif mech_vents_data['vent_type'] in ("Intermittent MEV", "Centralised continuous MEV", "Decentralised continuous MEV"):
                    self.__mech_vents_dict[mech_vents_name] = MechanicalVentilation(
                        self.__external_conditions,
                        mech_vents_data['sup_air_flw_ctrl'],
                        mech_vents_data['sup_air_temp_ctrl'],
                        0, #mech_vents_data['design_zone_cooling_covered_by_mech_vent'],
                        0, #mech_vents_data['design_zone_heating_covered_by_mech_vent'],
                        mech_vents_data['vent_type'],
                        mech_vents_data['SFP'],
                        mech_vents_data['design_outdoor_air_flow_rate'],
                        self.__simtime,
                        energy_supply_conn,
                        total_volume,
                        proj_dict['InfiltrationVentilation']['altitude'],
                        ctrl_intermittent_MEV,
                    )
                else:
                    sys.exit("Mechanical ventilation type not recognised")

                # TODO not all dwellings have mech vents - update to make mech vents optional
                if mech_vents_data['vent_type'] == "MVHR":
                    self.__space_heating_ductwork_dict = {}
                    if mech_vents_name not in self.__space_heating_ductwork_dict:
                        self.__space_heating_ductwork_dict[mech_vents_name] = []
                    for ductwork_data in mech_vents_data["ductwork"]:
                        if ductwork_data['cross_section_shape'] == 'circular':
                            duct_perimeter = None
                            internal_diameter = ductwork_data['internal_diameter_mm'] / units.mm_per_m
                            external_diameter = ductwork_data['external_diameter_mm'] / units.mm_per_m
                        elif ductwork_data['cross_section_shape'] == 'rectangular':
                            duct_perimeter = ductwork_data['duct_perimeter_mm'] / units.mm_per_m
                            internal_diameter = None
                            external_diameter = None
                        else:
                            sys.exit("Duct shape not valid")

                        ductwork = Ductwork(
                            ductwork_data['cross_section_shape'],
                            duct_perimeter,
                            internal_diameter,
                            external_diameter,
                            ductwork_data['length'],
                            ductwork_data['insulation_thermal_conductivity'],
                            ductwork_data['insulation_thickness_mm'] / units.mm_per_m,
                            ductwork_data['reflective'],
                            ductwork_data['duct_type'],
                            mech_vents_data['mvhr_location'],
                            mech_vents_data['mvhr_eff'],
                            )
                        self.__space_heating_ductwork_dict[mech_vents_name].append(ductwork)

        combustion_appliances_dict = {}
        for combustion_appliances_name, combustion_appliances_data in proj_dict['InfiltrationVentilation']['CombustionAppliances'].items():
            combustion_appliances_dict[combustion_appliances_name] = CombustionAppliances(
                combustion_appliances_data['supply_situation'],
                combustion_appliances_data['exhaust_situation'],
                combustion_appliances_data['fuel_type'],
                combustion_appliances_data['appliance_type'],
            )
        self.__ventilation = InfiltrationVentilation(
            self.__external_conditions,
            self.__simtime,
            f_cross,
            shield_class,
            terrain_class,
            average_pitch,
            windows_dict,
            vents_dict,
            leaks_dict,
            combustion_appliances_dict,
            atds_dict,
            self.__mech_vents_dict,
            self.__detailed_output_heating_cooling,
            altitude,
            total_volume,
            )
        
        if 'required_vent' in proj_dict['Control']:
            req_vent_dict = proj_dict['Control']['required_vent']
            self.__sched_req_vent = expand_schedule(float, req_vent_dict['schedule'], "main", True)
            self.__req_vent_start_day = req_vent_dict['start_day']
            self.__req_vent_time_series_step = req_vent_dict['time_series_step']
        else:
            req_vent_dict = None
            self.__sched_req_vent = None

        def dict_to_thermal_bridging(data):
            # If data is for individual thermal bridges, initialise the relevant
            # objects and return a list of them. Otherwise, just use the overall
            # figure given.
            if isinstance(data, dict):
                thermal_bridging = []
                for tb_name, tb_data in data.items():
                    tb_type = tb_data['type']
                    if tb_type == 'ThermalBridgeLinear':
                        tb = ThermalBridgeLinear(
                                tb_data['linear_thermal_transmittance'],
                                tb_data['length']
                                )
                    elif tb_type == 'ThermalBridgePoint':
                        tb = ThermalBridgePoint(tb_data['heat_transfer_coeff'])
                    else:
                        sys.exit( tb_name + ': thermal bridge type ('
                                + tb_type + ') not recognised.' )
                        # TODO Exit just the current case instead of whole program entirely?
                    thermal_bridging.append(tb)
            else:
                thermal_bridging = data
            return thermal_bridging

        self.__heat_system_name_for_zone = {}
        self.__cool_system_name_for_zone = {}
        opening_area_total = 0.0
        for z_data in proj_dict['Zone'].values():
            for building_element_data in z_data['BuildingElement'].values():
                if building_element_data['type'] == 'BuildingElementTransparent':
                    opening_area_total \
                        += building_element_data['height'] * building_element_data['width']

        def dict_to_zone(name, data):
            # Record which heating and cooling system this zone is heated/cooled by (if applicable)
            if 'SpaceHeatSystem' in data:
                # Check that no heating system has been assigned to more than one zone
                # TODO Adapt this check to work for multiple systems per zone
                if data['SpaceHeatSystem'] in self.__heat_system_name_for_zone.values():
                    sys.exit('Invalid input: SpaceHeatSystem (' + data['SpaceHeatSystem'] 
                           + ') has been assigned to more than one Zone')
                self.__heat_system_name_for_zone[name] = data['SpaceHeatSystem']
            else:
                self.__heat_system_name_for_zone[name] = None
            if 'SpaceCoolSystem' in data:
                # Check that no cooling system has been assigned to more than one zone
                # TODO Adapt this check to work for multiple systems per zone
                if data['SpaceCoolSystem'] in self.__cool_system_name_for_zone.values():
                    sys.exit('Invalid input: SpaceCoolSystem (' + data['SpaceCoolSystem'] 
                           + ') has been assigned to more than one Zone')
                self.__cool_system_name_for_zone[name] = data['SpaceCoolSystem']
            else:
                self.__cool_system_name_for_zone[name] = None

            # Read in building elements and add to list
            building_elements = []
            for building_element_name, building_element_data in data['BuildingElement'].items():
                building_elements.append(
                    dict_to_building_element(building_element_name, building_element_data)
                    )

            # Read in thermal bridging data
            thermal_bridging = dict_to_thermal_bridging(data['ThermalBridging'])

            return Zone(
                data['area'],
                data['volume'],
                building_elements,
                thermal_bridging,
                self.__ventilation,
                self.__external_conditions.air_temp(),
                data['temp_setpnt_init'],
                self.__control_WindowAdjust,
                print_heat_balance = print_heat_balance,
                use_fast_solver = use_fast_solver,
                )

        self.__zones = {}
        self.__energy_supply_conn_unmet_demand_zone = {}
        for name, data in proj_dict['Zone'].items():
            self.__zones[name] = dict_to_zone(name, data)
            self.__energy_supply_conn_unmet_demand_zone[name] \
                = self.__energy_supplies['_unmet_demand'].connection(name)

        self.__total_floor_area = sum(zone.area() for zone in self.__zones.values())
        self.__total_volume = sum(zone.volume() for zone in self.__zones.values())

        def convert_energy_to_Wm2(data):
            # Convert energy supplied to appliances from W to W / m2
            total_energy_supply = []
            for energy_data in expand_schedule(float, data['schedule'], "main", False):
                total_energy_supply.append(energy_data / self.__total_floor_area)
            return total_energy_supply
        
        self.__internal_gains = {}
        if 'InternalGains' in proj_dict:
            for name, data in proj_dict['InternalGains'].items():
                self.__internal_gains[name] = InternalGains(
                                                convert_energy_to_Wm2(data),
                                                self.__simtime,
                                                data['start_day'],
                                                data['time_series_step']
                                                )
        
        # Add internal gains from applicances to the internal gains dictionary and
        # create an energy supply connection for appliances
        for name, data in proj_dict['ApplianceGains'].items():
            energy_supply = self.__energy_supplies[data['EnergySupply']]
            # TODO Need to handle error if EnergySupply name is invalid.
            energy_supply_conn = energy_supply.connection(name)
            if "Events" in data.keys() and "Standby" in data.keys():
                self.__internal_gains[name] = EventApplianceGains(
                                                 energy_supply_conn,
                                                 self.__simtime,
                                                 data,
                                                 self.__total_floor_area
                                                 )
            else:
                # Convert energy supplied to appliances from W to W / m2
                self.__internal_gains[name] = ApplianceGains(
                                                 convert_energy_to_Wm2(data),
                                                 energy_supply_conn,
                                                 data['gains_fraction'],
                                                 self.__simtime,
                                                 data['start_day'],
                                                 data['time_series_step']
                                                 )
            
        # Where wet distribution heat source provide more than one service, some
        # calculations can only be performed after all services have been
        # calculated. Give these systems a timestep_end function and add these
        # systems to the following list, which will be iterated over later.
        self.__timestep_end_calcs = []

        def dict_to_heat_source_wet(name, data):
            heat_source_type = data['type']
            if heat_source_type == 'HeatPump':
                if SourceType.is_exhaust_air(data['source_type']):
                    # Check that ventilation system exists and  is compatible with exhaust air HP
                    if self.__mech_vents_dict is not None:
                        for mech_vent in self.__mech_vents_dict.values():
                            if mech_vent.vent_type == 'Intermittent MEV' or mech_vent.vent_type == 'Decentralised continuous MEV':
                                sys.exit('Exhaust air heat pump does not work with Intermittent MEV or Decentralised continuous MEV.')
                            else:
                                throughput_exhaust_air = \
                                    mech_vent.design_outdoor_air_flow_rate_m3_h
                                # Record heat source as potentially requiring overventilation
                                self.__heat_source_wet_names_requiring_overvent.append(name)
                                # Only pass in proj_obj into HeatPump class for exhaust air hps
                                proj_obj = self
                    else:
                        sys.exit('Ventilation object does not exist')
                else:
                    throughput_exhaust_air = None
                    proj_obj = None

                if SourceType.from_string(data['source_type']) == SourceType.HEAT_NETWORK:
                    energy_supply_HN = self.__energy_supplies[data['EnergySupply_heat_network']]
                    # TODO Check that EnergySupply object representing heat source
                    #      has an appropriate fuel type
                else:
                    energy_supply_HN = None


                if 'boiler' in data:
                    energy_supply_boiler = self.__energy_supplies[data['boiler']['EnergySupply']]
                    energy_supply_aux_boiler = self.__energy_supplies[data['boiler']['EnergySupply_aux']]
                    energy_supply_conn_aux_boiler = energy_supply_aux_boiler.connection('Boiler_auxiliary: ' + name)
                    boiler = Boiler(
                        data['boiler'],
                        energy_supply_boiler,
                        energy_supply_conn_aux_boiler,
                        self.__simtime,
                        self.__external_conditions,
                        )

                    if 'cost_schedule_hybrid'in data['boiler']:
                        cost_schedule_hybrid_hp = data['boiler']['cost_schedule_hybrid']
                    else:
                        cost_schedule_hybrid_hp = None

                    self.__timestep_end_calcs.append(boiler)
                else:
                    cost_schedule_hybrid_hp = None
                    boiler = None
                energy_supply = self.__energy_supplies[data['EnergySupply']]
                energy_supply_conn_name_auxiliary = 'HeatPump_auxiliary: ' + name
                heat_source = HeatPump(
                    data,
                    energy_supply,
                    energy_supply_conn_name_auxiliary,
                    self.__simtime,
                    self.__external_conditions,
                    len(proj_dict['Zone'].items()),
                    throughput_exhaust_air,
                    energy_supply_HN,
                    self.__detailed_output_heating_cooling,
                    boiler,
                    cost_schedule_hybrid_hp,
                    proj_obj,
                    )
                self.__timestep_end_calcs.append(heat_source)
                if 'BufferTank' in data :
                    self.__heat_sources_wet_with_buffer_tank.append(name)
            elif heat_source_type == 'Boiler':
                energy_supply = self.__energy_supplies[data['EnergySupply']]
                energy_supply_aux = self.__energy_supplies[data['EnergySupply_aux']]
                energy_supply_conn_aux = energy_supply_aux.connection('Boiler_auxiliary: ' + name)
                heat_source = Boiler(
                    data,
                    energy_supply,
                    energy_supply_conn_aux,
                    self.__simtime,
                    self.__external_conditions,
                    )
                self.__timestep_end_calcs.append(heat_source)
            elif heat_source_type == 'HIU':
                energy_supply = self.__energy_supplies[data['EnergySupply']]
                energy_supply_conn_name_auxiliary = 'HeatNetwork_auxiliary: ' + name
                energy_supply_conn_name_building_level_distribution_losses \
                    = 'HeatNetwork_building_level_distribution_losses: ' + name
                heat_source = HeatNetwork(
                    data['power_max'],
                    data['HIU_daily_loss'],
                    data['building_level_distribution_losses'],
                    energy_supply,
                    energy_supply_conn_name_auxiliary,
                    energy_supply_conn_name_building_level_distribution_losses,
                    self.__simtime,
                    )
                self.__timestep_end_calcs.append(heat_source)
                # Create list of internal gains for each hour of the year, in W / m2
                internal_gains_HIU = [heat_source.HIU_loss() \
                                        * units.W_per_kW \
                                        / self.__total_floor_area]
                total_internal_gains_HIU = internal_gains_HIU * units.days_per_year * units.hours_per_day
                # Append internal gains object to self.__internal_gains dictionary
                if name in self.__internal_gains.keys():
                    sys.exit('Name of HIU duplicates name of an existing InternalGains object')
                self.__internal_gains[name] = InternalGains(
                    total_internal_gains_HIU,
                    self.__simtime,
                    0, # Start day of internal gains time series
                    1.0, # Timestep of internal gains time series
                    )
            elif heat_source_type == 'HeatBattery':
                energy_supply = self.__energy_supplies[data['EnergySupply']]
                # TODO Need to handle error if EnergySupply name is invalid.
                energy_supply_conn = energy_supply.connection(name)
                charge_control: ToUChargeControl = self.__controls[data['ControlCharge']]
                heat_source = HeatBattery(
                    data,
                    charge_control,
                    energy_supply,
                    energy_supply_conn,
                    self.__simtime,
                    self.__external_conditions,
                    )
                self.__timestep_end_calcs.append(heat_source)
            else:
                sys.exit(name + ': heat source type (' \
                       + heat_source_type + ') not recognised.')
                # TODO Exit just the current case instead of whole program entirely?
            return heat_source

        # If one or more wet distribution heat sources have been provided, add them to the project
        self.__heat_sources_wet = {}
        self.__heat_sources_wet_with_buffer_tank = []
        self.__heat_source_wet_names_requiring_overvent = []
        # If no wet distribution heat sources have been provided, then skip.
        if 'HeatSourceWet' in proj_dict:
            for name, data in proj_dict['HeatSourceWet'].items():
                self.__heat_sources_wet[name] = dict_to_heat_source_wet(name, data)

        def dict_to_heat_source(
                name,
                data,
                cold_water_source,
                temp_setpoint,
                volume,
                daily_losses,
                heat_exchanger_surface_area=None,
                ):
            """ Parse dictionary of heat source data and return approprate heat source object """
            if 'Control' in data.keys():
                ctrl = self.__controls[data['Control']]
                # TODO Need to handle error if Control name is invalid.
            else:
                ctrl = None

            heat_source_type = data['type']
            if heat_source_type == 'ImmersionHeater':
                energy_supply = self.__energy_supplies[data['EnergySupply']]
                # TODO Need to handle error if EnergySupply name is invalid.
                energy_supply_conn_name = name
                energy_supply_conn = energy_supply.connection(name)

                heat_source = ImmersionHeater(
                    data['power'],
                    energy_supply_conn,
                    self.__simtime,
                    ctrl,
                    )
            elif heat_source_type == 'SolarThermalSystem':
                energy_supply = self.__energy_supplies[data['EnergySupply']]
                # TODO Need to handle error if EnergySupply name is invalid.
                energy_supply_conn_name = name
                energy_supply_conn = energy_supply.connection(name)

                heat_source = SolarThermalSystem(
                    data['sol_loc'],
                    data['area_module'],
                    data['modules'],
                    data['peak_collector_efficiency'],
                    data['incidence_angle_modifier'],
                    data['first_order_hlc'],
                    data['second_order_hlc'],
                    data['collector_mass_flow_rate'],
                    data['power_pump'],
                    data['power_pump_control'],
                    energy_supply_conn,
                    data['tilt'],
                    self.__init_orientation(data['orientation360']),
                    data['solar_loop_piping_hlc'],
                    self.__external_conditions,
                    self.__simtime,
                    self,
                    )
                
            elif heat_source_type == 'HeatSourceWet':
                energy_supply_conn_name = data['name'] + '_water_heating'

                heat_source_wet = self.__heat_sources_wet[data['name']]
                if isinstance(heat_source_wet, HeatPump):
                    heat_source = heat_source_wet.create_service_hot_water(
                        energy_supply_conn_name,
                        temp_setpoint,
                        data['temp_flow_limit_upper'],
                        cold_water_source,
                        ctrl,
                        )
                elif isinstance(heat_source_wet, Boiler):
                    heat_source = heat_source_wet.create_service_hot_water_regular(
                        energy_supply_conn_name,
                        temp_setpoint,
                        cold_water_source,
                        ctrl,
                        )
                elif isinstance(heat_source_wet, HeatNetwork):
                    # Add heat network hot water service for feeding hot water cylinder
                    heat_source = heat_source_wet.create_service_hot_water_storage(
                        energy_supply_conn_name,
                        temp_setpoint,
                        ctrl,
                        )
                elif isinstance(heat_source_wet, HeatBattery):
                    heat_source = heat_source_wet.create_service_hot_water_regular(
                        data,
                        energy_supply_conn_name,
                        temp_setpoint,
                        cold_water_source,
                        ctrl,
                        )
                else:
                    sys.exit(name + ': HeatSource type not recognised')
                    # TODO Exit just the current case instead of whole program entirely?
            elif heat_source_type == 'HeatPump_HWOnly':
                energy_supply = self.__energy_supplies[data['EnergySupply']]
                # TODO Need to handle error if EnergySupply name is invalid.
                energy_supply_conn_name = name
                energy_supply_conn = energy_supply.connection(name)

                heat_source = HeatPump_HWOnly(
                    data['power_max'],
                    data['test_data'],
                    data['vol_hw_daily_average'],
                    volume,
                    daily_losses,
                    heat_exchanger_surface_area,
                    data['in_use_factor_mismatch'],
                    data['tank_volume_declared'],
                    data['heat_exchanger_surface_area_declared'],
                    data['daily_losses_declared'],
                    energy_supply_conn,
                    self.__simtime,
                    ctrl,
                    )
            else:
                sys.exit(name + ': heat source type (' + heat_source_type + ') not recognised.')
                # TODO Exit just the current case instead of whole program entirely?
            return heat_source, energy_supply_conn_name

        # List of diverter objects (for end-of-timestep calculations
        self.__diverters = []

        def dict_to_hot_water_source(name, data):
            """ Parse dictionary of HW source data and return approprate HW source object """
            energy_supply_conn_names = []

            hw_source_type = data['type']
            if hw_source_type == 'StorageTank':
                cold_water_source = self.__cold_water_sources[data['ColdWaterSource']]
                # TODO Need to handle error if ColdWaterSource name is invalid.
                # TODO assuming here there is only one WWHRS
                if self.__wwhrs is not None:
                    for wwhrs_name in self.__wwhrs:
                        if isinstance(self.__wwhrs[wwhrs_name], wwhrs.WWHRS_InstantaneousSystemC) \
                        or isinstance(self.__wwhrs[wwhrs_name], wwhrs.WWHRS_InstantaneousSystemA):
                            cold_water_source = self.__wwhrs[wwhrs_name]

                if 'primary_pipework' in data:
                    primary_pipework_lst = data['primary_pipework']
                    for pipework_data in primary_pipework_lst:
                        pipework_data['internal_diameter'] \
                            = pipework_data['internal_diameter_mm'] / units.mm_per_m
                        pipework_data['external_diameter'] \
                            = pipework_data['external_diameter_mm'] / units.mm_per_m
                        pipework_data['insulation_thickness'] \
                            = pipework_data['insulation_thickness_mm'] / units.mm_per_m
                else:
                    primary_pipework_lst = []


                heat_source_dict= {}
                heat_source_names_dict = {}
                for heat_source_name, heat_source_data in data['HeatSource'].items():

                    # heat exchanger area
                    if 'HeatPump_HWOnly' in heat_source_data['type']:
                        heat_exchanger_surface_area = data['heat_exchanger_surface_area']
                    else:
                        heat_exchanger_surface_area = None

                    heat_source, conn_name = dict_to_heat_source(
                        heat_source_name,
                        heat_source_data,
                        cold_water_source,
                        data['setpoint_temp'],
                        data['volume'],
                        data['daily_losses'],
                        heat_exchanger_surface_area,
                        )
                    heat_source_dict[heat_source] = heat_source_data['heater_position'], \
                                                    heat_source_data['thermostat_position']
                    energy_supply_conn_names.append(conn_name)
                    heat_source_names_dict[heat_source_name] = heat_source

                if 'Control_hold_at_setpnt' in data:
                    ctrl_hold_at_setpnt = self.__controls[data['Control_hold_at_setpnt']]
                else:
                    ctrl_hold_at_setpnt = None

                hw_source = StorageTank(
                    data['volume'],
                    data['daily_losses'],
                    data['min_temp'],
                    data['setpoint_temp'],
                    cold_water_source,
                    self.__simtime,
                    heat_source_dict,
                    self,
                    self.__external_conditions,
                    24,
                    primary_pipework_lst,
                    energy_supply_unmet_demand.connection(name),
                    ctrl_hold_at_setpnt
                    )
                energy_supply_conn_names.append(name)

                for heat_source_name, heat_source_data in data['HeatSource'].items():
                    energy_supply_name = heat_source_data['EnergySupply']
                    if energy_supply_name in diverters \
                    and diverters[energy_supply_name]['StorageTank'] == name \
                    and diverters[energy_supply_name]['HeatSource'] == heat_source_name:
                        energy_supply = self.__energy_supplies[heat_source_data['EnergySupply']]
                        pv_diverter = PVDiverter(
                            hw_source,
                            heat_source_names_dict[heat_source_name],
                            )
                        energy_supply.connect_diverter(pv_diverter)
                        self.__diverters.append(pv_diverter)

            elif hw_source_type == 'CombiBoiler':
                energy_supply_conn_name = data['HeatSourceWet'] + '_water_heating'
                energy_supply_conn_names.append(energy_supply_conn_name)
                cold_water_source = self.__cold_water_sources[data['ColdWaterSource']]
                hw_source = self.__heat_sources_wet[data['HeatSourceWet']].create_service_hot_water_combi(
                    data,
                    energy_supply_conn_name,
                    data['setpoint_temp'],
                    cold_water_source
                    )
            elif hw_source_type == 'PointOfUse':
                energy_supply = self.__energy_supplies[data['EnergySupply']]
                # TODO Need to handle error if EnergySupply name is invalid.
                energy_supply_conn_name = name
                energy_supply_conn_names.append(energy_supply_conn_name)
                energy_supply_conn = energy_supply.connection(name)

                cold_water_source = self.__cold_water_sources[data['ColdWaterSource']]
                hw_source = PointOfUse(
                    data['efficiency'],
                    energy_supply_conn,
                    self.__simtime,
                    cold_water_source,
                    data['setpoint_temp']
                )
            elif hw_source_type == 'HIU':
                energy_supply_conn_name = data['HeatSourceWet'] + '_water_heating'
                energy_supply_conn_names.append(energy_supply_conn_name)
                cold_water_source = self.__cold_water_sources[data['ColdWaterSource']]
                hw_source = self.__heat_sources_wet[data['HeatSourceWet']].create_service_hot_water_direct(
                    energy_supply_conn_name,
                    data['setpoint_temp'],
                    cold_water_source,
                    )
            elif hw_source_type == 'HeatBattery':
                # TODO MC - add PCM heat battery in here
                pass
            else:
                sys.exit(name + ': hot water source type (' + hw_source_type + ') not recognised.')
                # TODO Exit just the current case instead of whole program entirely?
            return hw_source, energy_supply_conn_names

        self.__hot_water_sources = {}
        self.__energy_supply_conn_names_for_hot_water_source = {}
        for name, data in proj_dict['HotWaterSource'].items():
            self.__hot_water_sources[name], \
                self.__energy_supply_conn_names_for_hot_water_source[name] \
                = dict_to_hot_water_source(name, data)

        # Some systems (e.g. exhaust air heat pumps) may require overventilation
        # so initialise an empty list to hold the names of these systems
        self.__heat_system_names_requiring_overvent = []

        def dict_to_space_heat_system(name, data):
            space_heater_type = data['type']

            # Setting flag on the existence of a buffer tank in the emitter's loop
            # This is currently only considered in the case of HP as heat source.
            with_buffer_tank = False
            
            # ElecStorageHeater needs extra controllers
            if space_heater_type == 'ElecStorageHeater' and 'ControlCharger' in data.keys():
                charge_control = self.__controls[data['ControlCharger']]

            if 'Control' in data.keys():
                ctrl = self.__controls[data['Control']]
                # TODO Need to handle error if Control name is invalid.
            else:
                ctrl = None

            if space_heater_type == 'InstantElecHeater':
                energy_supply = self.__energy_supplies[data['EnergySupply']]
                # TODO Need to handle error if EnergySupply name is invalid.
                energy_supply_conn_name = name
                energy_supply_conn = energy_supply.connection(name)

                space_heater = InstantElecHeater(
                    data['rated_power'],
                    data['frac_convective'],
                    energy_supply_conn,
                    self.__simtime,
                    ctrl,
                    )
            elif space_heater_type == 'ElecStorageHeater':
                energy_supply = self.__energy_supplies[data['EnergySupply']]
                # TODO Need to handle error if EnergySupply name is invalid.
                energy_supply_conn_name = name
                energy_supply_conn = energy_supply.connection(name)

                space_heater = ElecStorageHeater(
                    data['rated_power'],
                    data['rated_power_instant'],
                    data['air_flow_type'],
                    data['temp_dis_safe'],
                    data['thermal_mass'],
                    data['frac_convective'],
                    data['U_ins'],
                    data['temp_charge_cut'],
                    data['mass_core'],
                    data['c_pcore'],
                    data['temp_core_target'],
                    data['A_core'],
                    data['c_wall'],
                    data['n_wall'],
                    data['thermal_mass_wall'],
                    data['fan_pwr'],
                    data['n_units'],
                    self.__zones[data['Zone']],
                    energy_supply_conn,
                    self.__simtime,
                    ctrl,
                    charge_control,
                )
            elif space_heater_type == 'WetDistribution':
                energy_supply_conn_name = data['HeatSource']['name'] + '_space_heating: ' + name
                heat_source = self.__heat_sources_wet[data['HeatSource']['name']]
                if isinstance(heat_source, HeatPump):
                    # TODO If EAHP, feed zone volume into function below
                    
                    # For HPs, checking if there's a buffer tank to inform both the service space heating
                    # and the emitters of its presence.
                    if data['HeatSource']['name'] in self.__heat_sources_wet_with_buffer_tank:
                        with_buffer_tank = True
                    
                    volume_heated = self.__total_volume_heated_by_system(name)
                    heat_source_service = heat_source.create_service_space_heating(
                        energy_supply_conn_name,
                        data['HeatSource']['temp_flow_limit_upper'],
                        data['temp_diff_emit_dsgn'],
                        ctrl,
                        volume_heated,
                        )
                    if heat_source.source_is_exhaust_air():
                        # Record heating system as potentially requiring overventilation
                        self.__heat_system_names_requiring_overvent.append(name)

                elif isinstance(heat_source, Boiler):
                    heat_source_service = heat_source.create_service_space_heating(
                        energy_supply_conn_name,
                        ctrl,
                        )
                elif isinstance(heat_source, HeatNetwork):
                    heat_source_service = heat_source.create_service_space_heating(
                        energy_supply_conn_name,
                        ctrl,
                        )
                elif isinstance(heat_source, HeatBattery):
                    heat_source_service = heat_source.create_service_space_heating(
                        energy_supply_conn_name,
                        ctrl,
                        )
                else:
                    sys.exit(name + ': HeatSource type not recognised')
                    # TODO Exit just the current case instead of whole program entirely?

                space_heater = Emitters(
                    data['thermal_mass'],
                    data['c'],
                    data['n'],
                    data['temp_diff_emit_dsgn'],
                    data['frac_convective'],
                    heat_source_service,
                    self.__zones[data['Zone']],
                    self.__external_conditions,
                    data['ecodesign_controller'],
                    data['design_flow_temp'],
                    self.__simtime,
                    self.__detailed_output_heating_cooling,
                    with_buffer_tank
                    )
            elif space_heater_type == 'WarmAir':
                energy_supply_conn_name = data['HeatSource']['name'] + '_space_heating: ' + name
                heat_source = self.__heat_sources_wet[data['HeatSource']['name']]
                if isinstance(heat_source, HeatPump):
                    volume_heated = self.__total_volume_heated_by_system(name)
                    space_heater = heat_source.create_service_space_heating_warm_air(
                        energy_supply_conn_name,
                        ctrl,
                        data['frac_convective'],
                        volume_heated,
                        )
                    if heat_source.source_is_exhaust_air():
                        # Record heating system as potentially requiring overventilation
                        self.__heat_system_names_requiring_overvent.append(name)
                else:
                    sys.exit(name + ': HeatSource type not recognised')
                    # TODO Exit just the current case instead of whole program entirely?
            else:
                sys.exit(name + ': space heating system type (' \
                       + space_heater_type + ') not recognised.')
                # TODO Exit just the current case instead of whole program entirely?

            return space_heater, energy_supply_conn_name

        # If one or more space heating systems have been provided, add them to the project
        self.__space_heat_systems = {}
        self.__energy_supply_conn_name_for_space_heat_system = {}
        # If no space heating systems have been provided, then skip. This
        # facilitates running the simulation with no heating systems at all
        if 'SpaceHeatSystem' in proj_dict:
            for name, data in proj_dict['SpaceHeatSystem'].items():
                # Only initialise systems that are actually used
                # TODO Adapt check to work with multiple systems per zone
                if name not in self.__heat_system_name_for_zone.values():
                    # TODO Add warning message here. Not adding now because it may break web interface
                    pass
                else:
                    self.__space_heat_systems[name], \
                        self.__energy_supply_conn_name_for_space_heat_system[name] \
                        = dict_to_space_heat_system(name, data)

        def dict_to_space_cool_system(name, data):
            if 'Control' in data.keys():
                ctrl = self.__controls[data['Control']]
                # TODO Need to handle error if Control name is invalid.
            else:
                ctrl = None

            cooling_system_type = data['type']
            if cooling_system_type == 'AirConditioning':
                energy_supply = self.__energy_supplies[data['EnergySupply']]
                # TODO Need to handle error if EnergySupply name is invalid.
                energy_supply_conn_name = name
                energy_supply_conn = energy_supply.connection(name)

                cooling_system = AirConditioning(
                   data['cooling_capacity'],
                   data['efficiency'],
                   data['frac_convective'],
                   energy_supply_conn,
                   self.__simtime,
                   ctrl,
                   )
            else:
                sys.exit(name + ': CoolSystem type not recognised')

            return cooling_system, energy_supply_conn_name

        self.__space_cool_systems = {}
        self.__energy_supply_conn_name_for_space_cool_system = {}
        # If no space cooling systems have been provided, then skip. This
        # facilitates running the simulation with no cooling systems at all
        if 'SpaceCoolSystem' in proj_dict:
            for name, data in proj_dict['SpaceCoolSystem'].items():
                # Only initialise systems that are actually used
                if name not in self.__cool_system_name_for_zone.values():
                    # TODO Add warning message here. Not adding now because it may break web interface
                    pass
                else:
                    self.__space_cool_systems[name], \
                        self.__energy_supply_conn_name_for_space_cool_system[name] \
                        = dict_to_space_cool_system(name, data)

        def dict_to_on_site_generation(name, data):
            """ Parse dictionary of on site generation data and
                return approprate on site generation object """
            on_site_generation_type = data['type']
            if on_site_generation_type == 'PhotovoltaicSystem':

                energy_supply = self.__energy_supplies[data['EnergySupply']]
                # TODO Need to handle error if EnergySupply name is invalid.
                energy_supply_conn = energy_supply.connection(name)

                pv_system = PhotovoltaicSystem(
                    data['peak_power'],
                    data['ventilation_strategy'],
                    data['pitch'],
                    self.__init_orientation(data['orientation360']),
                    data['base_height'], 
                    data['height'],
                    data['width'],
                    self.__external_conditions,
                    energy_supply_conn,
                    self.__simtime,
                    data["shading"],
                    data['inverter_peak_power'],
                    data['inverter_is_inside'],
                    )
            else:
                sys.exit(name + ': on site generation type ('
                         + on_site_generation_type + ') not recognised.')
                # TODO Exit just the current case instead of whole program entirely?
            return pv_system

        self.__on_site_generation = {}
        # If no on site generation have been provided, then skip.
        if 'OnSiteGeneration' in proj_dict:
            for name, data in proj_dict['OnSiteGeneration'].items():
                self.__on_site_generation[name] = dict_to_on_site_generation(name, data)

    def __init_resistance_or_uvalue(self, name, data):
        """ Return thermal resistance of construction (r_c) based on alternative inputs

        User will either provide r_c directly or provide u_value which needs to be converted
        """
        # If r_c has been provided directly, then use it, and print warning if
        # u_value has been provided in addition
        if 'r_c' in data.keys():
            if 'u_value' in data.keys():
                print( 'Warning: For BuildingElement input object "' \
                     + name + '" both r_c and u_value have been provided. ' \
                     + 'The value for r_c will be used.' \
                     )
            return data['r_c']
        # If only u_value has been provided, use it to calculate r_c
        else:
            if 'u_value' in data.keys():
                return BuildingElement.convert_uvalue_to_resistance(data['u_value'], data['pitch'])
            else:
                sys.exit( 'Error: For BuildingElement input object "' \
                        + name + '" neither r_c nor u_value have been provided.' \
                        )

    def __init_orientation(self, orientation360):
        """ Convert orientation from 0-360 (clockwise) to -180 to +180 (anticlockwise) """
        return 180 - orientation360

    def total_floor_area(self):
        return self.__total_floor_area

    def __total_volume_heated_by_system(self, heat_system_name):
        # TODO Adapt if condition to work with multiple systems per zone
        return sum(
            zone.volume()
            for z_name, zone in self.__zones.items()
            if self.__heat_system_name_for_zone[z_name] == heat_system_name
            )

    def calc_HTC_HLP(self):
        """ Calculate heat transfer coefficient (HTC) and heat loss parameter (HLP)
        according to the SAP10.2 specification """

        HTC_dict = {}
        HLP_dict = {}

        # Calculate the total fabric heat loss, total heat capacity, total ventilation heat
        # loss and total heat transfer coeffient for thermal bridges across all zones
        for z_name, zone in self.__zones.items():
            fabric_heat_loss = zone.total_fabric_heat_loss()
            thermal_bridges = zone.total_thermal_bridges()
            vent_heat_loss = zone.total_vent_heat_loss()

            # Calculate the heat transfer coefficent (HTC), in W / K
            # TODO check ventilation losses are correct
            HTC = fabric_heat_loss + thermal_bridges + vent_heat_loss

            # Calculate the HLP, in W / m2 K
            HLP = HTC / zone.area()

            HTC_dict[z_name] = HTC
            HLP_dict[z_name] = HLP

        total_HTC = sum(HTC_dict.values())
        total_HLP = total_HTC / self.__total_floor_area
        
        return total_HTC, total_HLP, HTC_dict, HLP_dict

    def calc_HCP(self):
        """ Calculate the total heat capacity normalised for floor area """
        # TODO party walls and solid doors should be exluded according to SAP spec - if party walls are
        # assumed to be ZTU building elements this could be set to zero?

        # Initialise variable
        total_heat_capacity = 0

        # Calculate the total heat capacity and total zone area
        for z_name, zone in self.__zones.items():
            total_heat_capacity += zone.total_heat_capacity()

        # Calculate the thermal mass parameter, in kJ / m2 K
        HCP = total_heat_capacity / self.__total_floor_area

        return HCP

    def calc_HLFF(self):
        "Calculate the heat loss form factor, defined as exposed area / floor area"

        total_heat_loss_area = 0
        for z_name, zone in self.__zones.items():
            total_heat_loss_area += zone.total_heat_loss_area()
        HLFF = total_heat_loss_area / self.__total_floor_area
        return HLFF

    def temp_internal_air(self):
        # Initialise internal air temperature and total area of all zones
        internal_air_temperature = 0

        # TODO here we are treating overall indoor temperature as average of all zones
        for z_name, zone in self.__zones.items():
            internal_air_temperature += zone.temp_internal_air() * zone.volume()

        internal_air_temperature /= self.__total_volume # average internal temperature
        return internal_air_temperature

    def __pipework_losses_and_internal_gains_from_hw(
            self,
            delta_t_h,
            vol_hot_water_at_tapping_point,
            hw_duration,
            no_of_hw_events,
            temp_hot_water,
            ):

        pw_losses_internal, pw_losses_external \
            = self.__calc_pipework_losses(
                delta_t_h,
                hw_duration,
                no_of_hw_events,
                temp_hot_water,
                )

        gains_internal_dhw_use \
            = self.frac_dhw_energy_internal_gains \
            * misc.water_demand_to_kWh(
                vol_hot_water_at_tapping_point,
                temp_hot_water,
                self.temp_internal_air(),
                )

        # Return:
        # - losses from internal distribution pipework (kWh)
        # - losses from external distribution pipework (kWh)
        # - internal gains due to hot water use (kWh)
        return pw_losses_internal, pw_losses_external, gains_internal_dhw_use

    def __pipework_losses_and_internal_gains_from_hw_StorageTank(
            self,
            delta_t_h,
            volume_water_remove_from_tank,
            hw_duration,
            no_of_hw_events,
            temp_final_drawoff,
            temp_average_drawoff,
            temp_hot_water,
            vol_hot_water_equiv_elec_shower,
            ):

        pw_losses_internal, pw_losses_external \
            = self.__calc_pipework_losses(
                delta_t_h,
                hw_duration,
                no_of_hw_events,
                temp_final_drawoff,
                )

        gains_internal_dhw_use_StorageTank \
            = self.frac_dhw_energy_internal_gains \
            * misc.water_demand_to_kWh(
                volume_water_remove_from_tank,
                temp_average_drawoff,
                self.temp_internal_air(),
                )

        gains_internal_dhw_use_IES \
            = self.frac_dhw_energy_internal_gains \
            * misc.water_demand_to_kWh(
                vol_hot_water_equiv_elec_shower,
                temp_hot_water,
                self.temp_internal_air(),
                )

        gains_internal_dhw_use = gains_internal_dhw_use_StorageTank + \
                                 gains_internal_dhw_use_IES

        # Return:
        # - losses from internal distribution pipework (kWh)
        # - losses from external distribution pipework (kWh)
        # - internal gains due to hot water use (kWh)
        return pw_losses_internal, pw_losses_external, gains_internal_dhw_use

    def __calc_pipework_losses(self, delta_t_h, hw_duration, no_of_hw_events, temp_hot_water):
        # sum up all hw_demand and allocate pipework losses also.
        # hw_demand is volume.

        demand_water_temperature = temp_hot_water
        internal_air_temperature = self.temp_internal_air()
        external_air_temperature = self.__external_conditions.air_temp()

        return self.__dhw_demand.calc_pipework_losses(
            delta_t_h,
            hw_duration,
            no_of_hw_events,
            demand_water_temperature,
            internal_air_temperature,
            external_air_temperature,
            )

    def __calc_internal_gains_buffer_tank(self):
        """ Calculate the losses in the buffer tank """
        internal_gains_buffer_tanks = 0.0
        for h_name in self.__heat_sources_wet_with_buffer_tank:
            internal_gains_buffer_tanks += self.__heat_sources_wet[h_name].buffer_int_gains()
            
        return internal_gains_buffer_tanks

    def __calc_internal_gains_ductwork(self):
        """ Calculate the losses/gains in the MVHR ductwork """
        if self.__space_heating_ductwork_dict is None:
            return 0.0

        internal_gains_ductwork_watts = 0.0
        for mvhr_ductwork in self.__space_heating_ductwork_dict.values():
            # assume MVHR unit is running 100% of the time
            for duct in mvhr_ductwork:
                if duct.get_duct_type() in (DuctType.INTAKE, DuctType.EXHAUST):
                    # Heat loss from intake or exhaust ducts is to zone, so add
                    # to internal gains (may be negative gains)
                    internal_gains_ductwork_watts += \
                        duct.total_duct_heat_loss(self.temp_internal_air(), self.__external_conditions.air_temp())
                elif duct.get_duct_type() in (DuctType.SUPPLY, DuctType.EXTRACT):
                    # Heat loss from supply and extract ducts is to outside, so
                    # subtract from internal gains
                    internal_gains_ductwork_watts -= \
                        duct.total_duct_heat_loss(self.temp_internal_air(), self.__external_conditions.air_temp())

        return internal_gains_ductwork_watts

    def __space_heat_internal_gains_for_zone(
            self,
            zone: Zone,
            gains_internal_dhw: float,
            internal_gains_ductwork_per_m3: float,
            gains_internal_buffer_tank: float,
            ):
        # Initialise to dhw internal gains split proportionally to zone floor area
        gains_internal_zone = (gains_internal_buffer_tank + gains_internal_dhw) * zone.area() / self.__total_floor_area
        for internal_gains_name, internal_gains_object in self.__internal_gains.items():
            gains_internal_zone \
                += internal_gains_object.total_internal_gain(zone.area())

        # Add gains from ventilation fans (also calculates elec demand from fans)
        # TODO Remove the branch on the type of ventilation (find a better way)
        if self.__mech_vents_dict is not None:
            for mech_vent in self.__mech_vents_dict.values():
                gains_internal_zone += mech_vent.fans(zone.volume(), self.__total_volume)
                gains_internal_zone += internal_gains_ductwork_per_m3 * zone.volume()
        return gains_internal_zone

    def __get_heat_cool_systems_for_zone(self, z_name):
        """ Look up relevant heating and cooling systems for the specified zone """
        # TODO For now, the existing single system inputs are each added to a
        #      list. This will eventually need to handle a list being specified
        #      in the inputs.
        h_name_list = [self.__heat_system_name_for_zone[z_name]]
        c_name_list = [self.__cool_system_name_for_zone[z_name]]

        temp_setpnt_heat_system, temp_setpnt_cool_system, \
            frac_convective_heat_system, frac_convective_cool_system \
            = self.__get_setpoints_and_convective_fractions(h_name_list, c_name_list)

        # Sort heating and cooling systems by setpoint (highest first for
        # heating, lowest first for cooling)
        # TODO In the event of two systems having the same setpoint, make
        #      sure the one listed first by the user takes priority
        h_name_list_sorted \
            = sorted(temp_setpnt_heat_system, key=lambda x: temp_setpnt_heat_system[x], reverse=True)
        c_name_list_sorted \
            = sorted(temp_setpnt_cool_system, key=lambda x: temp_setpnt_cool_system[x], reverse=False)

        return \
            h_name_list_sorted, c_name_list_sorted, \
            temp_setpnt_heat_system, temp_setpnt_cool_system, \
            frac_convective_heat_system, frac_convective_cool_system

    def __get_setpoints_and_convective_fractions(self, h_name_list, c_name_list):
        """ Look up convective fractions and setpoints for heating/cooling """
        # Use default setpoints when there is no heat/cool system or
        # there is no setpoint for the current timestep
        frac_convective_heat = {}
        frac_convective_cool = {}
        temp_setpnt_heat = {}
        temp_setpnt_cool = {}
        for h_name in h_name_list:
            if h_name is not None:
                frac_convective_heat[h_name] = self.__space_heat_systems[h_name].frac_convective()
                temp_setpnt_heat[h_name] = self.__space_heat_systems[h_name].temp_setpnt()
            else:
                frac_convective_heat[h_name] = 1.0
                temp_setpnt_heat[h_name] = None
            if temp_setpnt_heat[h_name] is None:
                temp_setpnt_heat[h_name] = temp_setpnt_heat_none

        for c_name in c_name_list:
            if c_name is not None:
                frac_convective_cool[c_name] = self.__space_cool_systems[c_name].frac_convective()
                temp_setpnt_cool[c_name] = self.__space_cool_systems[c_name].temp_setpnt()
            else:
                frac_convective_cool[c_name] = 1.0
                temp_setpnt_cool[c_name] = None
            if temp_setpnt_cool[c_name] is None:
                temp_setpnt_cool[c_name] = temp_setpnt_cool_none

        return temp_setpnt_heat, temp_setpnt_cool, frac_convective_heat, frac_convective_cool

    def __gains_heat_cool(self, delta_t_h, hc_output_convective, hc_output_radiative):
        gains_heat_cool_convective \
            = sum(hc_output_convective.values()) * units.W_per_kW / delta_t_h
        gains_heat_cool_radiative \
            = sum(hc_output_radiative.values()) * units.W_per_kW / delta_t_h
        return gains_heat_cool_convective, gains_heat_cool_radiative
    
    def __calc_air_changes_per_hour(self, temp_int_air, R_w_arg, intial_p_z_ref_guess, reporting_flag):
        """Calculate the incoming air changes per hour
           intial_p_z_ref_guess is used for calculation in first timestep.
           Later timesteps use the previous timesteps p_z_ref of max and min ACH ,respective to calc.
        
        arg R_w_arg --
        arg intial_p_z_ref_guess -- 
        arg reporting_flag -- 
        """
        if self.__initial_loop:
            self.__internal_pressure_window[reporting_flag] = self.__ventilation.calculate_internal_reference_pressure(
                                                           intial_p_z_ref_guess,
                                                           temp_int_air,
                                                           R_w_arg,
                                                           )
        else:
            self.__internal_pressure_window[reporting_flag] = self.__ventilation.calculate_internal_reference_pressure(
                                                           self.__internal_pressure_window[reporting_flag],
                                                           temp_int_air,
                                                           R_w_arg,
                                                           )    
            
        incoming_air_flow = self.__ventilation.incoming_air_flow(
                                                    self.__internal_pressure_window[reporting_flag],
                                                    temp_int_air,
                                                    R_w_arg,
                                                    reporting_flag,
                                                    report_effective_flow_rate = True,
                                                    )
        air_changes_per_hour = incoming_air_flow / self.__total_volume
        return air_changes_per_hour

    def __calc_space_heating(self, delta_t_h: float, gains_internal_dhw: float):
        """ Calculate space heating demand, heating system output and temperatures

        Arguments:
        delta_t_h -- calculation timestep, in hours
        gains_internal_dhw -- internal gains from hot water system for this timestep, in W
        """
        temp_ext_air = self.__external_conditions.air_temp()
        temp_int_air = self.temp_internal_air()
        # Calculate timestep in seconds
        delta_t = delta_t_h * units.seconds_per_hour

        internal_gains_ductwork = self.__calc_internal_gains_ductwork()
        internal_gains_ductwork_per_m3 = internal_gains_ductwork / self.__total_volume

        internal_gains_buffer_tank = self.__calc_internal_gains_buffer_tank()

        # Windows shut
        ach_windows_shut = self.__calc_air_changes_per_hour(
            temp_int_air,
            R_w_arg = 0,
            intial_p_z_ref_guess = 0,
            reporting_flag = 'min',
            )

        #Windows fully open
        ach_windows_open = self.__calc_air_changes_per_hour(
            temp_int_air,
            R_w_arg = 1,
            intial_p_z_ref_guess = 0,
            reporting_flag = 'max',
            )

        # To indicate the future loop should involve the p_Z_ref from previous calc
        self.__initial_loop = False
        
        if self.__sched_req_vent is not None: 
            ach_target = self.__sched_req_vent[self.__simtime.time_series_idx(
                self.__req_vent_start_day,
                self.__req_vent_time_series_step
                )]

            ach_target = max(ach_windows_shut , min(ach_target, ach_windows_open))
        else:
            ach_target = ach_windows_shut

        gains_internal_zone = {}
        gains_solar_zone = {}
        h_name_list_sorted_zone = {}
        c_name_list_sorted_zone = {}
        temp_setpnt_heat_zone_system = {}
        temp_setpnt_cool_zone_system = {}
        frac_convective_heat_zone_system = {}
        frac_convective_cool_zone_system = {}
        ach_cooling_zone = {}
        ach_to_trigger_heating_zone = {}
        internal_air_temp = {}
        operative_temp = {}
        space_heat_demand_zone = {}
        space_cool_demand_zone = {}
        space_heat_demand_system = {}
        space_cool_demand_system = {}
        space_heat_provided_system = {}
        space_cool_provided_system = {}
        heat_balance_dict = {}

        #Average supply temperature 
        avg_air_supply_temp = self.__ventilation.temp_supply()

        for z_name, zone in self.__zones.items():
            # Calculate internal and solar gains
            gains_internal_zone[z_name] = self.__space_heat_internal_gains_for_zone(
                zone,
                gains_internal_dhw,
                internal_gains_ductwork_per_m3,
                internal_gains_buffer_tank,
                )
            gains_solar_zone[z_name] = zone.gains_solar()

            # Get heating and cooling characteristics for the current zone
            h_name_list_sorted_zone[z_name], c_name_list_sorted_zone[z_name], \
                temp_setpnt_heat_zone_system[z_name], temp_setpnt_cool_zone_system[z_name], \
                frac_convective_heat_zone_system[z_name], frac_convective_cool_zone_system[z_name] \
                = self.__get_heat_cool_systems_for_zone(z_name)

            # Calculate space heating demand based on highest-priority systems,
            # assuming no output from any other systems
            space_heat_demand_zone[z_name], space_cool_demand_zone[z_name], \
                ach_cooling_zone[z_name], ach_to_trigger_heating_zone[z_name] \
                = zone.space_heat_cool_demand(
                    delta_t_h,
                    temp_ext_air,
                    gains_internal_zone[z_name],
                    gains_solar_zone[z_name],
                    frac_convective_heat_zone_system[z_name][h_name_list_sorted_zone[z_name][0]],
                    frac_convective_cool_zone_system[z_name][c_name_list_sorted_zone[z_name][0]],
                    temp_setpnt_heat_zone_system[z_name][h_name_list_sorted_zone[z_name][0]],
                    temp_setpnt_cool_zone_system[z_name][c_name_list_sorted_zone[z_name][0]],
                    avg_air_supply_temp = avg_air_supply_temp,
                    ach_windows_open = ach_windows_open,
                    ach_target = ach_target,
                    )

        # Ventilation required, including for cooling
        is_heating_demand = False
        for demand in space_heat_demand_zone.values():
            if demand > 0.0:
                is_heating_demand = True
        is_cooling_demand = False
        for demand in space_cool_demand_zone.values():
            if demand < 0.0:
                is_cooling_demand = True
        if is_heating_demand:
            # Do not open windows any further than required for ventilation
            # requirement if there is any heating demand in any zone
            ach_cooling = ach_target
        elif is_cooling_demand:
            # Do not open windows any further than required for ventilation
            # requirement if there is any cooling demand in any zone
            ach_cooling = ach_target

            # In this case, will need to recalculate space cooling demand for
            # each zone, this time assuming no window opening for all zones
            # TODO There might be a way to make this more efficient and reduce
            #      the number of times the heat balance solver has to run, but
            #      this would require a wider refactoring of the zone module's
            #      space_heat_cool_demand function
            for z_name, zone in self.__zones.items():
                space_heat_demand_zone[z_name], space_cool_demand_zone[z_name], _, _ \
                    = zone.space_heat_cool_demand(
                        delta_t_h,
                        temp_ext_air,
                        gains_internal_zone[z_name],
                        gains_solar_zone[z_name],
                        frac_convective_heat_zone_system[z_name][h_name_list_sorted_zone[z_name][0]],
                        frac_convective_cool_zone_system[z_name][c_name_list_sorted_zone[z_name][0]],
                        temp_setpnt_heat_zone_system[z_name][h_name_list_sorted_zone[z_name][0]],
                        temp_setpnt_cool_zone_system[z_name][c_name_list_sorted_zone[z_name][0]],
                        avg_air_supply_temp = avg_air_supply_temp,
                        ach_cooling = ach_cooling,
                        )
        else:
            # Subject to the above/below limits, take the maximum required window
            # opening from across all the zones
            ach_cooling = max(ach_cooling_zone.values())

            # Do not open windows to an extent where it would cause any zone
            # temperature to fall below the heating setpoint for that zone
            ach_to_trigger_heating_list = [
                x for x in ach_to_trigger_heating_zone.values() \
                if x is not None
                ]
            if len(ach_to_trigger_heating_list) > 0:
                ach_cooling = min(min(ach_to_trigger_heating_list), ach_cooling)

            # Do not reduce air change rate below ventilation requirement even
            # if it would help with temperature regulation
            ach_cooling = max(ach_cooling, ach_target)

        # Calculate heating/cooling system response and temperature achieved in each zone
        for z_name, zone in self.__zones.items():
            # Check for name clash between heating and cooling systems
            if set(h_name_list_sorted_zone[z_name]).intersection(c_name_list_sorted_zone[z_name]):
                sys.exit('All heating and cooling systems must have unique names')
            # TODO Populate the lists below with minimum output of each heating and cooling system
            hc_output_convective = {
                hc_name: 0.0
                for hc_name in h_name_list_sorted_zone[z_name] + c_name_list_sorted_zone[z_name]
                }
            hc_output_radiative = {
                hc_name: 0.0
                for hc_name in h_name_list_sorted_zone[z_name] + c_name_list_sorted_zone[z_name]
                }

            space_heat_demand_zone_system \
                = {h_name: 0.0 for h_name in h_name_list_sorted_zone[z_name]}
            space_cool_demand_zone_system \
                = {c_name: 0.0 for c_name in c_name_list_sorted_zone[z_name]}
            space_heat_provided_zone_system \
                = {h_name: 0.0 for h_name in h_name_list_sorted_zone[z_name]}
            space_cool_provided_zone_system \
                = {c_name: 0.0 for c_name in c_name_list_sorted_zone[z_name]}
            h_idx = 0
            c_idx = 0
            space_heat_running_time_cumulative = 0.0
            while(    h_idx < len(h_name_list_sorted_zone[z_name])
                  and c_idx < len(c_name_list_sorted_zone[z_name])
                 ):
                h_name = h_name_list_sorted_zone[z_name][h_idx]
                c_name = c_name_list_sorted_zone[z_name][c_idx]
                frac_convective_heat = frac_convective_heat_zone_system[z_name][h_name]
                frac_convective_cool = frac_convective_cool_zone_system[z_name][c_name]
                temp_setpnt_heat = temp_setpnt_heat_zone_system[z_name][h_name]
                temp_setpnt_cool = temp_setpnt_cool_zone_system[z_name][c_name]

                # Calculate space heating/cooling demand, accounting for any
                # output from systems (either output already calculated for
                # higher-priority systems, or min output of current or
                # lower-priority systems).
                # TODO This doesn't yet handle minimum output from current or
                #      lower-priority systems
                gains_heat_cool_convective, gains_heat_cool_radiative \
                    = self.__gains_heat_cool(delta_t_h, hc_output_convective, hc_output_radiative)
                if gains_heat_cool_convective == 0.0 and gains_heat_cool_radiative == 0.0:
                    # If there is no output from any systems, then don't need to
                    # calculate demand again
                    space_heat_demand_zone_system[h_name] = space_heat_demand_zone[z_name]
                    space_cool_demand_zone_system[c_name] = space_cool_demand_zone[z_name]
                else:
                    space_heat_demand_zone_system[h_name], space_cool_demand_zone_system[c_name], \
                        ach_cooling_zone[z_name], _ \
                        = zone.space_heat_cool_demand(
                            delta_t_h,
                            temp_ext_air,
                            gains_internal_zone[z_name],
                            gains_solar_zone[z_name],
                            frac_convective_heat,
                            frac_convective_cool,
                            temp_setpnt_heat,
                            temp_setpnt_cool,
                            avg_air_supply_temp = avg_air_supply_temp,
                            gains_heat_cool_convective = gains_heat_cool_convective,
                            gains_heat_cool_radiative = gains_heat_cool_radiative,
                            ach_cooling = ach_cooling,
                            )

                # If any heating systems potentially require overventilation,
                # calculate running time and throughput factor for current service
                # based on space heating demand assuming only overventilation
                # required for DHW
                # TODO For now, disable overventilation calculation. Only used
                #      for exhaust air heat pump when conditions are out of
                #      range of test data, which is disallowed for now. May also
                #      be needed for combustion appliances in the future but
                #      none of these have been implemented yet. Calculations
                #      using throughput factor have been removed to simplify
                #      ventilation calculations
                '''
                if h_name in self.__heat_system_names_requiring_overvent:
                    time_running_space, throughput_factor_zone \
                        = self.__space_heat_systems[h_name].running_time_throughput_factor(
                            space_heat_demand_zone_system[h_name],
                            space_heat_running_time_cumulative,
                            )
                    # Add running time for the current space heating service to
                    # the cumulative total
                    space_heat_running_time_cumulative += time_running_space

                    # Combine throughput factors for space and water heating.
                    # Note that any additional ventilation due to water heating
                    # is apportioned to the zones in proportion to their volume,
                    # while any additional ventilation due to space heating is
                    # assigned to the zone that the heating demand originates
                    # from, in order to simplify the calculation. If the
                    # additional ventilation for space heating was assigned on a
                    # whole-dwelling basis, then handling any system requiring
                    # overventilation would require recalculation of demand and
                    # system response for the other zone, if it has already been
                    # calculated. This may mean rolling back (or finding a way
                    # to calculate but not commit) system calculations until all
                    # system calculations have been done for both zones.
                    # TODO Make sure this works with more than one system
                    #      requiring overventilation.
                    throughput_factor_zone_overall \
                        = (throughput_factor_zone - 1.0) \
                        + (throughput_factor_dhw - 1.0) \
                        + 1.0

                # If there is overventilation due to heating or hot water system (e.g.
                # exhaust air heat pump) then recalculate space heating/cooling demand
                # with additional ventilation calculated based on throughput factor
                # based on original space heating demand calculation. Note the
                # additional ventilation throughput is the result of the HP running
                # to satisfy both space and water heating demand but will affect
                # space heating demand only
                # TODO The space heating demand is only recalculated once, rather
                #      than feeding back in to the throughput factor calculation
                #      above to get a further-refined space heating demand. This is
                #      consistent with the approach in SAP 10.2 and keeps the
                #      execution time of the calculation bounded. However, the
                #      merits of iterating over this calculation until converging on
                #      a solution should be considered in the future.
                if throughput_factor_zone_overall > 1.0:
                    # Add additional gains from ventilation fans
                    # TODO Remove the branch on the type of ventilation (find a better way)
                    if self.__ventilation is not None \
                    and not isinstance(self.__ventilation, NaturalVentilation):
                        gains_internal_zone[z_name] \
                            += self.__ventilation.fans(
                                zone.volume(),
                                throughput_factor_zone_overall - 1.0,
                                )
                    space_heat_demand_zone_system[h_name], space_cool_demand_zone_system[c_name], _ \
                        = zone.space_heat_cool_demand(
                            delta_t_h,
                            temp_ext_air,
                            gains_internal_zone[z_name],
                            gains_solar_zone[z_name],
                            frac_convective_heat,
                            frac_convective_cool,
                            temp_setpnt_heat,
                            temp_setpnt_cool,
                            gains_heat_cool_convective = sum(hc_output_convective.values()) * units.W_per_kW / delta_t_h,
                            gains_heat_cool_radiative = sum(hc_output_radiative.values()) * units.W_per_kW / delta_t_h,
                            throughput_factor = throughput_factor_zone_overall,
                            )
                    # Need to recalculate space heating demand on zone to account
                    # for extra ventilation, for reporting purposes
                    space_heat_demand_zone[z_name], _, _ \
                        = zone.space_heat_cool_demand(
                            delta_t_h,
                            temp_ext_air,
                            gains_internal_zone[z_name],
                            gains_solar_zone[z_name],
                            frac_convective_heat_zone_system[z_name][h_name_list_sorted_zone[z_name][0]],
                            frac_convective_cool_zone_system[z_name][c_name_list_sorted_zone[z_name][0]],
                            temp_setpnt_heat_zone_system[z_name][h_name_list_sorted_zone[z_name][0]],
                            temp_setpnt_cool_zone_system[z_name][c_name_list_sorted_zone[z_name][0]],
                            throughput_factor = throughput_factor_zone_overall,
                            )
                '''

                # Calculate heating/cooling provided
                if space_heat_demand_zone_system[h_name] > 0.0:
                    space_heat_provided_zone_system[h_name] \
                        = self.__space_heat_systems[h_name].demand_energy(
                            space_heat_demand_zone_system[h_name],
                            )
                    hc_output_convective[h_name] \
                        = space_heat_provided_zone_system[h_name] * frac_convective_heat
                    hc_output_radiative[h_name] \
                        = space_heat_provided_zone_system[h_name] * (1.0 - frac_convective_heat)
                    # If heating has been provided, then next iteration of loop
                    # should use next-priority heating system
                    h_idx += 1
                if space_cool_demand_zone_system[c_name] < 0.0:
                    space_cool_provided_zone_system[c_name] \
                        = self.__space_cool_systems[c_name].demand_energy(
                            space_cool_demand_zone_system[c_name],
                            )
                    hc_output_convective[c_name] \
                        = space_cool_provided_zone_system[c_name] * frac_convective_cool
                    hc_output_radiative[c_name] \
                        = space_cool_provided_zone_system[c_name] * (1.0 - frac_convective_cool)
                    # If cooling has been provided, then next iteration of loop
                    # should use next-priority cooling system
                    c_idx += 1

                # Terminate loop if there is no more demand
                if (space_heat_demand_zone_system[h_name] <= 0.0 \
                and space_cool_demand_zone_system[c_name] >= 0.0) :
                    break
            
            # Call any remaining heating and cooling systems with zero demand
            for h_name in h_name_list_sorted_zone[z_name][h_idx:]:
                if h_name is not None:
                    space_heat_provided_zone_system[h_name] \
                        = self.__space_heat_systems[h_name].demand_energy(0.0)
                else:
                    space_heat_provided_zone_system[h_name] = 0.0
                hc_output_convective[h_name] \
                    = space_heat_provided_zone_system[h_name] * frac_convective_heat
                hc_output_radiative[h_name] \
                    = space_heat_provided_zone_system[h_name] * (1.0 - frac_convective_heat)
            for c_name in c_name_list_sorted_zone[z_name][c_idx:]:
                if c_name is not None:
                    space_cool_provided_zone_system[c_name] \
                        = self.__space_cool_systems[c_name].demand_energy(0.0)
                else:
                    space_cool_provided_zone_system[c_name] = 0.0
                hc_output_convective[c_name] \
                    = space_cool_provided_zone_system[c_name] * frac_convective_cool
                hc_output_radiative[c_name] \
                    = space_cool_provided_zone_system[c_name] * (1.0 - frac_convective_cool)

            # Calculate unmet demand
            self.__unmet_demand(
                delta_t_h,
                temp_ext_air,
                z_name,
                zone,
                gains_internal_zone[z_name],
                gains_solar_zone[z_name],
                temp_setpnt_heat_zone_system[z_name],
                temp_setpnt_cool_zone_system[z_name],
                frac_convective_heat_zone_system[z_name],
                frac_convective_cool_zone_system[z_name],
                h_name_list_sorted_zone[z_name],
                c_name_list_sorted_zone[z_name],
                space_heat_demand_zone[z_name],
                space_cool_demand_zone[z_name],
                hc_output_convective,
                hc_output_radiative,
                ach_windows_open,
                ach_target,
                avg_air_supply_temp,
                )

            # Sum heating gains (+ve) and cooling gains (-ve) and convert from kWh to W
            hc_output_convective_total = sum(hc_output_convective.values())
            hc_output_radiative_total = sum(hc_output_radiative.values())
            gains_heat_cool \
                = (hc_output_convective_total + hc_output_radiative_total) \
                * units.W_per_kW / delta_t_h
            if gains_heat_cool != 0.0:
                frac_convective \
                    = hc_output_convective_total \
                    / (hc_output_convective_total + hc_output_radiative_total)
            else:
                frac_convective = 1.0

            # Calculate final temperatures achieved
            heat_balance_dict[z_name] = zone.update_temperatures(
                delta_t,
                temp_ext_air,
                gains_internal_zone[z_name],
                gains_solar_zone[z_name],
                gains_heat_cool,
                frac_convective,
                ach_cooling,
                avg_air_supply_temp,
                )
            internal_air_temp[z_name] = zone.temp_internal_air()
            operative_temp[z_name] = zone.temp_operative()

            for h_name in h_name_list_sorted_zone[z_name]:
                if h_name not in space_heat_demand_system.keys():
                    space_heat_demand_system[h_name] = 0.0
                    space_heat_provided_system[h_name] = 0.0
                space_heat_demand_system[h_name] += space_heat_demand_zone_system[h_name]
                space_heat_provided_system[h_name] += space_heat_provided_zone_system[h_name]
            for c_name in c_name_list_sorted_zone[z_name]:
                if c_name not in space_cool_demand_system.keys():
                    space_cool_demand_system[c_name] = 0.0
                    space_cool_provided_system[c_name] = 0.0
                space_cool_demand_system[c_name] += space_cool_demand_zone_system[c_name]
                space_cool_provided_system[c_name] += space_cool_provided_zone_system[c_name]

        return \
            gains_internal_zone, gains_solar_zone, \
            operative_temp, internal_air_temp, \
            space_heat_demand_zone, space_cool_demand_zone, \
            space_heat_demand_system, space_cool_demand_system, \
            space_heat_provided_system, space_cool_provided_system, \
            internal_gains_ductwork, heat_balance_dict

    def __get_highest_priority_required_system(self, hc_name_list_sorted, space_heatcool_systems):
        """ Determine highest-priority system that is in its required heating
        or cooling period (and not just the setback period)

        Arguments:
        hc_name_list_sorted -- list of heating or cooling systems (not combined
                               list), sorted in order of priority
        space_heatcool_systems -- dict of space heating or cooling system objects
                                  (not combined list)
        """
        hc_name_highest_req = None
        for hc_name in hc_name_list_sorted:
            if hc_name is not None and space_heatcool_systems[hc_name].in_required_period():
                hc_name_highest_req = hc_name
                break
        return hc_name_highest_req

    def __unmet_demand(
            self,
            delta_t_h,
            temp_ext_air,
            z_name,
            zone,
            gains_internal,
            gains_solar,
            temp_setpnt_heat_system,
            temp_setpnt_cool_system,
            frac_convective_heat_system,
            frac_convective_cool_system,
            h_name_list_sorted,
            c_name_list_sorted,
            space_heat_demand,
            space_cool_demand,
            hc_output_convective,
            hc_output_radiative,
            ach_max,
            ach_target,
            avg_air_supply_temp
            ):
        """ Calculate how much space heating / cooling demand is unmet """
        # Note: Use demand calculated based on highest-priority systems
        # Note: Demand is not considered unmet if it is outside the
        #       required heating/cooling period (which does not include
        #       times when the system is on due to setback or advanced
        #       start). If different systems have different required
        #       heating/cooling periods, unmet demand will be based on the
        #       system with the highest setpoint, ignoring any systems that
        #       are not in required periods (e.g. systems that are in
        #       setback or advanced start periods.
        # Note: Need to check that demand is non-zero, to avoid
        #       reporting unmet demand when heating system is absorbing
        #       energy from zone or cooling system is releasing energy
        #       to zone, which may be the case in some timesteps for
        #       systems with significant thermal mass.

        # Determine highest-priority system that is in its required heating
        # or cooling period (and not just the setback period)
        h_name_highest_req = self.__get_highest_priority_required_system(
            h_name_list_sorted,
            self.__space_heat_systems,
            )
        c_name_highest_req = self.__get_highest_priority_required_system(
            c_name_list_sorted,
            self.__space_cool_systems,
            )

        gains_heat \
            = sum(hc_output_convective[h_name] + hc_output_radiative[h_name] \
                  for h_name in h_name_list_sorted \
                  )
        gains_cool \
            = sum(hc_output_convective[c_name] + hc_output_radiative[c_name] \
                  for c_name in c_name_list_sorted \
                  )
        energy_shortfall_heat = max(0, space_heat_demand - gains_heat)
        energy_shortfall_cool = max(0, -(space_cool_demand - gains_cool))

        if (   h_name_highest_req is not None \
           and space_heat_demand > 0.0 \
           and energy_shortfall_heat > 0.0 \
           ) \
        or (   c_name_highest_req is not None \
           and space_cool_demand < 0.0 \
           and energy_shortfall_cool > 0.0 \
           ):
            if energy_shortfall_heat > 0.0 and h_name_highest_req != h_name_list_sorted[0] \
            or energy_shortfall_cool > 0.0 and c_name_highest_req != c_name_list_sorted[0]:
                # If the highest-priority system is not in required heating
                # period, but a lower-priority system is, calculate demand
                # based on the highest-priority system that is in required
                # heating period

                # Handle case where no heating/cooling system is in required
                # period. In this case, there will be no heat output anyway so
                # the convective fraction doesn't matter
                if h_name_highest_req is None:
                    frac_convective_heat = 1.0
                    temp_setpnt_heat = temp_setpnt_heat_none
                else:
                    frac_convective_heat = frac_convective_heat_system[h_name_highest_req]
                    temp_setpnt_heat = temp_setpnt_heat_system[h_name_highest_req]
                if c_name_highest_req is None:
                    frac_convective_cool = 1.0
                    temp_setpnt_cool = temp_setpnt_cool_none
                else:
                    frac_convective_cool = frac_convective_cool_system[c_name_highest_req]
                    temp_setpnt_cool = temp_setpnt_cool_system[c_name_highest_req]

                space_heat_demand_req, space_cool_demand_req, _, _ = zone.space_heat_cool_demand(
                    delta_t_h, 
                    temp_ext_air, 
                    gains_internal, 
                    gains_solar, 
                    frac_convective_heat, 
                    frac_convective_cool, 
                    temp_setpnt_heat, 
                    temp_setpnt_cool,
                    avg_air_supply_temp = avg_air_supply_temp,
                    ach_windows_open = ach_max,
                    ach_target = ach_target,
                    )
                unmet_demand_heat = max(0, space_heat_demand_req - gains_heat)
                unmet_demand_cool = max(0, -(space_cool_demand_req - gains_cool))
            else:
                # If highest-priority system is in required heating period,
                # use the demand already calculated for the zone
                unmet_demand_heat = energy_shortfall_heat
                unmet_demand_cool = energy_shortfall_cool

            self.__energy_supply_conn_unmet_demand_zone[z_name].demand_energy(
                unmet_demand_heat + unmet_demand_cool,
                )

    def run(self):
        """ Run the simulation """
        timestep_array = []
        gains_internal_dict = {}
        gains_solar_dict = {}
        operative_temp_dict = {}
        internal_air_temp_dict = {}
        space_heat_demand_dict = {}
        space_cool_demand_dict = {}
        space_heat_demand_system_dict = {}
        space_cool_demand_system_dict = {}
        space_heat_provided_dict = {}
        space_cool_provided_dict = {}
        zone_list = []
        hot_water_demand_dict = {}
        hot_water_energy_demand_dict = {}
        hot_water_energy_demand_dict_incl_pipework = {}
        hot_water_energy_output_dict = {}
        hot_water_duration_dict = {}
        hot_water_no_events_dict = {}
        hot_water_pipework_dict = {}
        ductwork_gains_dict = {}
        hot_water_primary_pipework_dict = {}
        hot_water_storage_losses_dict = {}
        heat_balance_all_dict = {'air_node': {}, 'internal_boundary': {},'external_boundary': {}}
        heat_source_wet_results_dict = {}
        heat_source_wet_results_annual_dict = {}
        emitters_output_dict = {}
        vent_output_list = []

        for z_name in self.__zones.keys():
            gains_internal_dict[z_name] = []
            gains_solar_dict[z_name] = []
            operative_temp_dict[z_name] = []
            internal_air_temp_dict[z_name] = []
            space_heat_demand_dict[z_name] = []
            space_cool_demand_dict[z_name] = []
            zone_list.append(z_name)
            for hb_name in heat_balance_all_dict.keys():
                heat_balance_all_dict[hb_name][z_name] = {}

        for z_name, h_name in self.__heat_system_name_for_zone.items():
            space_heat_demand_system_dict[h_name] = []
            space_heat_provided_dict[h_name] = []

        for z_name, c_name in self.__cool_system_name_for_zone.items():
            space_cool_demand_system_dict[c_name] = []
            space_cool_provided_dict[c_name] = []

        hot_water_demand_dict['demand'] = []
        hot_water_energy_demand_dict['energy_demand'] = []
        hot_water_energy_demand_dict_incl_pipework['energy_demand_incl_pipework_loss'] = []
        hot_water_energy_output_dict['energy_output'] = []
        hot_water_duration_dict['duration'] = []
        hot_water_no_events_dict['no_events'] = []
        hot_water_pipework_dict['pw_losses'] = []
        ductwork_gains_dict['ductwork_gains'] = []
        hot_water_primary_pipework_dict['primary_pw_losses'] =[]
        hot_water_storage_losses_dict['storage_losses'] =[]
        self.__initial_loop = True
        self.__internal_pressure_window = {}

        # Loop over each timestep
        for t_idx, t_current, delta_t_h in self.__simtime:
            timestep_array.append(t_current)
            temp_hot_water = self.__hot_water_sources['hw cylinder'].get_temp_hot_water()
            temp_final_drawoff = temp_hot_water
            temp_average_drawoff = temp_hot_water
            hw_demand_vol, hw_demand_vol_target, hw_vol_at_tapping_points, hw_duration, no_events, \
                hw_energy_demand, usage_events, vol_hot_water_equiv_elec_shower \
                = self.__dhw_demand.hot_water_demand(t_idx, temp_hot_water)

            # TODO Remove hard-coding of hot water source name
            # TODO Reporting of the hot water energy output assumes that there
            #      is only one water heating system. If the model changes in
            #      future to allow more than one hot water system, this code may
            #      need to be revised to handle that scenario.

            if isinstance(self.__hot_water_sources['hw cylinder'], StorageTank):
                hw_energy_output, unmet_energy, temp_final_drawoff, temp_average_drawoff, \
                    volume_water_remove_from_tank \
                    = self.__hot_water_sources['hw cylinder'].demand_hot_water(usage_events)
            
                pw_losses_internal, pw_losses_external, gains_internal_dhw_use \
                    = self.__pipework_losses_and_internal_gains_from_hw_StorageTank(
                        delta_t_h,
                        volume_water_remove_from_tank,
                        hw_duration,
                        no_events,
                        temp_final_drawoff, 
                        temp_average_drawoff,
                        temp_hot_water,
                        vol_hot_water_equiv_elec_shower,
                        )
            else:
                hw_energy_output \
                    = self.__hot_water_sources['hw cylinder'].demand_hot_water(hw_demand_vol_target)
                
                pw_losses_internal, pw_losses_external, gains_internal_dhw_use \
                    = self.__pipework_losses_and_internal_gains_from_hw(
                        delta_t_h,
                        hw_vol_at_tapping_points,
                        hw_duration,
                        no_events,
                        temp_hot_water,
                        )

            # Convert from litres to kWh
            cold_water_source = self.__hot_water_sources['hw cylinder'].get_cold_water_source()
            cold_water_temperature = cold_water_source.temperature()
            hw_energy_demand_incl_pipework_loss = misc.water_demand_to_kWh(
                hw_demand_vol,
                temp_hot_water,
                cold_water_temperature,
                )

            gains_internal_dhw \
                = (pw_losses_internal + gains_internal_dhw_use) \
                * units.W_per_kW / self.__simtime.timestep()
            if isinstance(self.__hot_water_sources['hw cylinder'], StorageTank) \
            or isinstance(self.__hot_water_sources['hw cylinder'], BoilerServiceWaterCombi):
                gains_internal_dhw += self.__hot_water_sources['hw cylinder'].internal_gains()

            # loop through on-site energy generation
            for g_name, gen in self.__on_site_generation.items():
                pv = self.__on_site_generation[g_name]
                # Get energy produced for the current timestep
                energy_produced, energy_lost = pv.produce_energy()
                # Add the energy lost figure to the internal gains if it is considered inside the building
                if pv.inverter_is_inside():
                    gains_internal_dhw += energy_lost * units.W_per_kW / self.__simtime.timestep()

            # Addition of primary_pipework_losses_kWh for reporting as part of investigation of issue #31225: FDEV A082
            if isinstance(self.__hot_water_sources['hw cylinder'], StorageTank):
                primary_pw_losses, storage_losses = self.__hot_water_sources['hw cylinder'].toreport()
            else:
                primary_pw_losses = 0.0
                storage_losses = 0.0
            
            gains_internal_zone, gains_solar_zone, \
                operative_temp, internal_air_temp, \
                space_heat_demand_zone, space_cool_demand_zone, \
                space_heat_demand_system, space_cool_demand_system, \
                space_heat_provided, space_cool_provided, \
                ductwork_gains, heat_balance_dict \
                = self.__calc_space_heating(delta_t_h, gains_internal_dhw)

            # Perform calculations that can only be done after all heating
            # services have been calculated.
            for system in self.__timestep_end_calcs:
                system.timestep_end()

            for z_name, gains_internal in gains_internal_zone.items():
                gains_internal_dict[z_name].append(gains_internal)

            for z_name, gains_solar in gains_solar_zone.items():
                gains_solar_dict[z_name].append(gains_solar)

            for z_name, temp in operative_temp.items():
                operative_temp_dict[z_name].append(temp)

            for z_name, temp in internal_air_temp.items():
                internal_air_temp_dict[z_name].append(temp)

            for z_name, demand in space_heat_demand_zone.items():
                space_heat_demand_dict[z_name].append(demand)

            for z_name, demand in space_cool_demand_zone.items():
                space_cool_demand_dict[z_name].append(demand)

            for h_name, demand in space_heat_demand_system.items():
                space_heat_demand_system_dict[h_name].append(demand)

            for c_name, demand in space_cool_demand_system.items():
                space_cool_demand_system_dict[c_name].append(demand)

            for h_name, output in space_heat_provided.items():
                space_heat_provided_dict[h_name].append(output)

            for c_name, output in space_cool_provided.items():
                space_cool_provided_dict[c_name].append(output)

            for z_name, hb_dict in heat_balance_dict.items():
                if hb_dict is not None:
                    for hb_name, gains_losses_dict in hb_dict.items():
                        for heat_gains_losses_name, heat_gains_losses_value in gains_losses_dict.items():
                            if heat_gains_losses_name in heat_balance_all_dict[hb_name][z_name].keys():
                                heat_balance_all_dict[hb_name][z_name][heat_gains_losses_name].append(heat_gains_losses_value)
                            else:
                                heat_balance_all_dict[hb_name][z_name][heat_gains_losses_name] =[heat_gains_losses_value]

            hot_water_demand_dict['demand'].append(hw_demand_vol)
            hot_water_energy_demand_dict['energy_demand'].append(hw_energy_demand)
            hot_water_energy_demand_dict_incl_pipework['energy_demand_incl_pipework_loss'].append(hw_energy_demand_incl_pipework_loss)
            hot_water_energy_output_dict['energy_output'].append(hw_energy_output)
            hot_water_duration_dict['duration'].append(hw_duration)
            hot_water_no_events_dict['no_events'].append(no_events)
            hot_water_pipework_dict['pw_losses'].append(pw_losses_internal + pw_losses_external)
            ductwork_gains_dict['ductwork_gains'].append(ductwork_gains)
            hot_water_primary_pipework_dict['primary_pw_losses'].append(primary_pw_losses)
            hot_water_storage_losses_dict['storage_losses'].append(storage_losses)

            for _, supply in self.__energy_supplies.items():
                supply.calc_energy_import_export_betafactor()

            for diverter in self.__diverters:
                diverter.timestep_end()

        zone_dict = {
            'internal gains': gains_internal_dict,
            'solar gains': gains_solar_dict,
            'operative temp': operative_temp_dict,
            'internal air temp': internal_air_temp_dict,
            'space heat demand': space_heat_demand_dict,
            'space cool demand': space_cool_demand_dict,
            }
        hc_system_dict = {
            'Heating system': space_heat_demand_system_dict,
            'Heating system output': space_heat_provided_dict,
            'Cooling system': space_cool_demand_system_dict,            
            'Cooling system output': space_cool_provided_dict,
            }
        hot_water_dict = {
            'Hot water demand': hot_water_demand_dict,
            'Hot water energy demand incl pipework_loss': hot_water_energy_demand_dict_incl_pipework,
            'Hot water energy demand': hot_water_energy_demand_dict,
            'Hot water duration': hot_water_duration_dict,
            'Hot Water Events': hot_water_no_events_dict,
            'Pipework losses': hot_water_pipework_dict,
            'Primary pipework losses': hot_water_primary_pipework_dict,
            'Storage losses': hot_water_storage_losses_dict,
            }

        # Report detailed outputs from heat source wet objects, if requested and available
        # TODO Note that the below assumes that there is only one water
        #      heating service and therefore that all hot water energy
        #      output is assigned to that service. If the model changes in
        #      future to allow more than one hot water system, this code may
        #      need to be revised to handle that scenario.
        if self.__detailed_output_heating_cooling:
            for name, heat_source_wet in self.__heat_sources_wet.items():
                if  hasattr(heat_source_wet, "output_detailed_results") \
                and callable(heat_source_wet.output_detailed_results):
                    heat_source_wet_results_dict[name], heat_source_wet_results_annual_dict[name] \
                        = heat_source_wet.output_detailed_results(
                            hot_water_energy_output_dict['energy_output']
                            )
            # Emitter detailed output results are stored with respect to heat_system_name
            for heat_system_name, heat_system in self.__space_heat_systems.items():
                if hasattr(heat_system, 'output_emitter_results') \
                and callable(heat_system.output_emitter_results): 
                    emitters_output_dict[heat_system_name] = heat_system.output_emitter_results()

            # Detailed ventilation results collected from ventilation class function      
            if hasattr(self.__ventilation, 'output_vent_results') \
                and callable(self.__ventilation.output_vent_results): 
                    vent_output_list = self.__ventilation.output_vent_results()

        # Return results from all energy supplies
        results_totals = {}
        results_end_user = {}
        energy_import = {}
        energy_export = {}
        energy_generated_consumed = {}
        energy_to_storage = {}
        energy_from_storage = {}
        energy_diverted = {}
        betafactor = {}
        for name, supply in self.__energy_supplies.items():
            results_totals[name] = supply.results_total()
            results_end_user[name] = supply.results_by_end_user()
            energy_import[name] = supply.get_energy_import()
            energy_export[name] = supply.get_energy_export()
            energy_generated_consumed[name] = supply.get_energy_generated_consumed()
            energy_to_storage[name], energy_from_storage[name] = supply.get_energy_to_from_battery()
            energy_diverted[name] = supply.get_energy_diverted()
            betafactor[name] = supply.get_beta_factor()

        hot_water_energy_out = {'hw cylinder': hot_water_energy_output_dict['energy_output']}
        dhw_cop_dict = self.__heat_cool_cop(
            hot_water_energy_out,
            results_end_user,
            self.__energy_supply_conn_names_for_hot_water_source,
            )
        heat_cop_dict = self.__heat_cool_cop(
            space_heat_provided_dict,
            results_end_user,
            self.__energy_supply_conn_name_for_space_heat_system
            )
        cool_cop_dict = self.__heat_cool_cop(
            space_cool_provided_dict,
            results_end_user,
            self.__energy_supply_conn_name_for_space_cool_system
            )

        return \
            timestep_array, results_totals, results_end_user, \
            energy_import, energy_export, energy_generated_consumed, \
            energy_to_storage, energy_from_storage, energy_diverted, betafactor, \
            zone_dict, zone_list, hc_system_dict, hot_water_dict, \
            heat_cop_dict, cool_cop_dict, dhw_cop_dict, \
            ductwork_gains_dict, heat_balance_all_dict, \
            heat_source_wet_results_dict, heat_source_wet_results_annual_dict, \
            emitters_output_dict, vent_output_list

    def __heat_cool_cop(
            self,
            energy_provided_dict,
            results_end_user,
            energy_supply_conn_name_for_space_hc_system,
            ):
        """ Calculate overall CoP over calculation period for each heating and cooling system """
        # Loop over heating systems, get energy output and input, and calculate CoP
        hc_output_overall = {}
        hc_input_overall = {}
        cop_dict = {}
        for hc_name, hc_output in energy_provided_dict.items():
            if hc_name is None:
                continue
            # Take absolute value because cooling system output is reported as a negative value
            hc_output_overall[hc_name] = abs(sum(hc_output))
            hc_input_overall[hc_name] = 0.0
            energy_supply_conn_names = energy_supply_conn_name_for_space_hc_system[hc_name]
            if not isinstance(energy_supply_conn_names, list):
                energy_supply_conn_names = [energy_supply_conn_names]
            for fuel_name, fuel_summary in results_end_user.items():
                if fuel_name == '_unmet_demand':
                    continue
                for conn_name, energy_cons in fuel_summary.items():
                    if conn_name in energy_supply_conn_names:
                        hc_input_overall[hc_name] += sum(energy_cons)

            if hc_input_overall[hc_name] > 0:
                cop_dict[hc_name] = hc_output_overall[hc_name] / hc_input_overall[hc_name]
            else:
                cop_dict[hc_name] = 'DIV/0'

        return cop_dict

