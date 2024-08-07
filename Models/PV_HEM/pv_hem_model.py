#!/usr/bin/env python3

"""
This module contains objects that represent photovoltaic systems.
"""

# Standard library imports
import math
import numpy as np
import core.units as units
from core.space_heat_demand.building_element import projected_height
from core.external_conditions import ExternalConditions
from core.simulation_time import SimulationTime
from core.energy_supply.energy_supply import EnergySupply, EnergySupplyConnection


class PhotovoltaicSystem:
    """ An object to represent a photovoltaic system """

    """ system performance factor lookup
          informative values from table C.4 Annex C BS EN 15316-4-3:2017
          note from 6.2.4.7.2 rear surface free - is if PV system is not integrated.
          Assume this means NOT integrated (BIPV)  or attached (BAPV)
          Increased by 0.85/0.8 based on median quoted performance ratio from
          "Performance of Distributed PV in the UK: A Statistical Analysis of
          Over 7000 systems" conference paper from 31st European Photovoltaic
          Solar Energy Conference and Exhibition, September 2015, assuming this
          applies to moderately ventilated case.
    """
    # Note: BS EN 15316-4-3:2017 section 6.2.4.7.2 states that the performance
    #       factor for "rear surface free" should be 1.0. However, this would
    #       seem to imply that there are no inverter or other system losses,
    #       despite the fact that section 6.2.4.7.5 states that this factor
    #       accounts for these losses. Also, no factor for "rear surface free"
    #       has been given in Table C.4. Therefore, it was decided to use the
    #       same factor for "rear surface free" as for "strongly or forced
    #       ventilated".
    __f_perf_lookup = {
        'unventilated': 0.81,
        'moderately_ventilated': 0.85,
        'strongly_or_forced_ventilated': 0.87,
        'rear_surface_free': 0.87,
    }

    def __init__(self, peak_power, ventilation_strategy, pitch, orientation,
                 base_height, height, width,
                shading, inverter_peak_power, inverter_is_inside):
        """ Construct a PhotovoltaicSystem object

        Arguments:
        peak_power       -- Peak power in kW; represents the electrical power of a photovoltaic
                            system with a given area and a for a solar irradiance of 1 kW/m2
                            on this surface (at 25 degrees)
                            TODO - Could add other options at a later stage.
                            Standard has alternative method when peak power is not available
                            (input type of PV module and Area instead when peak power unknown)
        ventilation_strategy   -- ventilation strategy of the PV system.
                                  This will be used to determine the system performance factor
                                   based on a lookup table
        pitch            -- is the tilt angle (inclination) of the PV panel from horizontal,
                            measured upwards facing, 0 to 90, in degrees.
                            0=horizontal surface, 90=vertical surface.
                            Needed to calculate solar irradiation at the panel surface.
        orientation      -- is the orientation angle of the inclined surface, expressed as the
                            geographical azimuth angle of the horizontal projection of the inclined
                            surface normal, -180 to 180, in degrees;
                            Assumed N 180 or -180, E 90, S 0, W -90
                            TODO - PV standard refers to angle as between 0 to 360?
                            Needed to calculate solar irradiation at the panel surface.
        base_height      -- is the distance between the ground and the lowest edge of the PV panel, in m
        height           -- is the height of the PV panel, in m
        width            -- is the width of the PV panel, in m
        inverter_peak_power -- Peak power in kW; represents the peak electrical power input to the inverter
        inverter_is_inside -- tells us that the inverter is considered inside the building
        overshading      -- TODO could add at a later date. Feed into solar module
        """
        self.__peak_power = peak_power
        self.__f_perf = self.__f_perf_lookup[ventilation_strategy]
        self.__pitch = pitch
        self.__orientation = orientation
        self.__base_height = base_height
        self.__width = width
        self.__projected_height = projected_height(pitch, height)
        self.__shading = shading
        self.__inverter_peak_power = inverter_peak_power
        self.__inverter_is_inside = inverter_is_inside

    def shading_factors_direct_diffuse(self):
        """ return calculated shading factor """
        return self.__external_conditions.shading_reduction_factor_direct_diffuse( \
                self.__base_height, self.__projected_height, self.__width, \
                self.__pitch, self.__orientation, self.__shading)

    def inverter_is_inside(self):
        """ Return whether this unit is considered inside the building or not """
        return self.__inverter_is_inside

    def connect(self, air_temperatures, wind_speeds, wind_directions, diffuse_horizontal_radiation,
                direct_beam_radiation, 
                solar_reflectivity_of_ground):
        """ Produce electrical energy (in kWh) from the PV system
            according to BS EN 15316-4-3:2017 """
        
        simtime = SimulationTime(0, 1, 1)
        proj_dict = {
            "ExternalConditions": {
                "air_temperatures": air_temperatures,
                "wind_speeds": wind_speeds,
                "wind_directions": wind_directions,
                "diffuse_horizontal_radiation": diffuse_horizontal_radiation,
                "direct_beam_radiation": direct_beam_radiation,
                "solar_reflectivity_of_ground": solar_reflectivity_of_ground,
                "latitude": 51.42,
                "longitude": -0.75,
                "timezone": 0,
                "start_day": 0,
                "end_day": 0,
                "time_series_step": 1,
                "january_first": 1,
                "daylight_savings": "not applicable",
                "leap_day_included": False,
                "direct_beam_conversion_needed": False,
                "shading_segments":[{"number": 1, "start": 180, "end": 135},
                                    {"number": 2, "start": 135, "end": 90,
                                     "shading": [
                                         {"type": "overhang", "height": 2.2, "distance": 6}
                                         ]
                                     },
                                    {"number": 3, "start": 90, "end": 45},
                                    {"number": 4, "start": 45, "end": 0, 
                                     "shading": [
                                         {"type": "obstacle", "height": 40, "distance": 4},
                                         {"type": "overhang", "height": 3, "distance": 7}
                                         ]
                                     },
                                    {"number": 5, "start": 0, "end": -45,
                                     "shading": [
                                         {"type": "obstacle", "height": 3, "distance": 8},
                                         ]
                                     },
                                    {"number": 6, "start": -45, "end": -90},
                                    {"number": 7, "start": -90, "end": -135},
                                    {"number": 8, "start": -135, "end": -180}],
            }
        }
        external_conditions = ExternalConditions(
            simtime,
            proj_dict['ExternalConditions']['air_temperatures'],
            proj_dict['ExternalConditions']['wind_speeds'],
            proj_dict['ExternalConditions']['diffuse_horizontal_radiation'],
            proj_dict['ExternalConditions']['direct_beam_radiation'],
            proj_dict['ExternalConditions']['solar_reflectivity_of_ground'],
            proj_dict['ExternalConditions']['latitude'],
            proj_dict['ExternalConditions']['longitude'],
            proj_dict['ExternalConditions']['timezone'],
            proj_dict['ExternalConditions']['start_day'],
            proj_dict['ExternalConditions']['end_day'],
            proj_dict['ExternalConditions']["time_series_step"],
            proj_dict['ExternalConditions']['january_first'],
            proj_dict['ExternalConditions']['daylight_savings'],
            proj_dict['ExternalConditions']['leap_day_included'],
            proj_dict['ExternalConditions']['direct_beam_conversion_needed'],
            proj_dict['ExternalConditions']['shading_segments']
            )
        energysupply = EnergySupply("electricity", simtime)
        energysupplyconn = energysupply.connection("pv generation without shading")

        self.__simulation_time = simtime
        self.__external_conditions = external_conditions
        self.__energy_supply_conn = energysupplyconn

        #solar_irradiance in W/m2
        i_sol_dir, i_sol_dif, _ = self.__external_conditions.calculated_direct_diffuse_total_irradiance(
            self.__pitch,
            self.__orientation
            )
        #shading factors
        f_sh_dir, f_sh_dif = self.shading_factors_direct_diffuse()
        #solar_irradiation in kWh/m2
        solar_irradiation = (i_sol_dir * f_sh_dir + i_sol_dif * f_sh_dif) \
                            * self.__simulation_time.timestep() / units.W_per_kW

        #reference_solar_irradiance kW/m2
        ref_solar_irradiance = 1

        #CALCULATION
        #E.el.pv.out.h = E.sol.pv.h * P.pk * f.perf / I.ref
        #energy_produced = solar_irradiation * peak_power * system_performance_factor
        #                    / reference_solar_irradiance
        # energy input in kWh; now need to calculate total energy produce taking into account inverter efficiency
        energy_input \
            = solar_irradiation * self.__peak_power * self.__f_perf / 0.92 / ref_solar_irradiance 
            # f_perf is divided by 0.92 to avoid double-applying the inverter efficiency, 
            # which is applied separately below via 'inverter_dc_ac_efficiency', since 
            # inverter efficiency was inherently included in the factors taken from
            # from BS EN 15316-4-3:2017.
        
        # power output from PV panel in kW used to calculate ratio for efficiency loss of inverters from DC to AC
        power_input_inverter = energy_input / self.__simulation_time.timestep()

        # Calculate Ratio of Rated Power
        ratio_of_rated_output = min(power_input_inverter, self.__inverter_peak_power) / self.__inverter_peak_power

        # Using Ratio of Rated Power, calculate Inverter DC to AC efficiency 
        # equation was estimated based on graph from 
        # https://www.researchgate.net/publication/260286647_Performance_of_PV_inverters figure 9
        if ratio_of_rated_output == 0:
            inverter_dc_ac_efficiency = 0
        else:
            inverter_dc_ac_efficiency = 0.92 * (math.tanh(4.67375 * ratio_of_rated_output) ** 0.137951)

        # Calculate energy produced output taking into account peak power of inverter + array 
        # and inverter DC to AC efficiency
        energy_produced \
            = min (energy_input, self.__inverter_peak_power * self.__simulation_time.timestep()) \
                 * inverter_dc_ac_efficiency
        
        # Add energy produced to the applicable energy supply connection (this will reduce demand)
        self.__energy_supply_conn.supply_energy(energy_produced)

        # energy_lost = energy_input - energy_produced
    
        return {'pv_gen': energy_produced, 'total_irr': diffuse_horizontal_radiation + direct_beam_radiation}


if __name__ == "__main__":

    pv_system = PhotovoltaicSystem(
                2.5,                        # Peak power
                "moderately_ventilated",    # Ventilation strat
                30,                         # Pitch
                0,                          # Orientation
                10,                         # Base height
                2,                          # Height
                3,                          # Width
                [],                         # Shading
                2.5,                        # Peak power
                inverter_is_inside=False,
                )

    result = pv_system.connect(
         0.0,
         3.9,
         220,
         11,
         11,
         0.2
    )


