#!/usr/bin/env python3

"""
This module provides objects to represent ventilation elements.
"""

# Standard library imports
import sys
from enum import IntEnum
from math import sqrt

# Third-party imports
import numpy as np

# Local imports
from core.units import seconds_per_hour, litres_per_cubic_metre, W_per_kW,\
    Celcius2Kelvin

# Define constants
p_a = 1.204 # Air density at 20 degrees C, in kg/m^3 , BS EN ISO 52016-1:2017, Section 6.3.6
c_a = 1006.0 # Specific heat of air at constant pressure, in J/(kg K), BS EN ISO 52016-1:2017, Section 6.3.6


def air_change_rate_to_flow_rate(air_change_rate, zone_volume):
    """ Convert infiltration rate from ach to m^3/s """
    return air_change_rate * zone_volume / seconds_per_hour


# TODO Throughput factor only applies to MVHR and WHEV, therefore only these
#      systems accept throughput_factor as an argument to the h_ve function.
#      This means that the MVHR and WHEV classes no longer have the same
#      interface as other ventilation element classes, which could make future
#      development more difficult. Ideally, we would find a cleaner way to
#      implement this difference.

class VentilationElementInfiltration:
    """ A class to represent infiltration ventilation elements """

    # Infiltration rates for openings (m3 per hour)
    # TODO Reference these
    __INF_RATE_CHIMNEY_OPEN = 80.0
    __INF_RATE_CHIMNEY_BLOCKED = 20.0
    __INF_RATE_FLUE_OPEN = 20.0
    __INF_RATE_FLUE_SOLID_FUEL_BOILER = 20.0
    __INF_RATE_FLUE_OTHER_HEATER = 35.0
    __INF_RATE_FIRE_CLOSED = 10.0
    __INF_RATE_FIRE_GAS = 40.0
    __INF_RATE_EXTRACT_FAN = 10.0
    __INF_RATE_PASSIVE_STACK_VENT = 10.0

    class ShelterType(IntEnum):
        # Values match column indices in table of divisors
        VERY_SHELTERED = 0
        SHELTERED = 1
        NORMAL = 2
        EXPOSED = 3

        @classmethod
        def from_string(cls, strval):
            if strval == "very sheltered":
                return cls.VERY_SHELTERED
            elif strval == "sheltered":
                return cls.SHELTERED
            elif strval == "normal":
                return cls.NORMAL
            elif strval == "exposed":
                return cls.EXPOSED
            else:
                sys.exit('ShelterType (' + str(strval) + ') not valid.')
                # TODO Exit just the current case instead of whole program entirely?

    class DwellingType(IntEnum):
        # Values match row indices in table of divisors
        HOUSE_1_STOREY = 0
        HOUSE_2_STOREY = 1
        FLAT_STOREY_1_TO_5 = 2
        FLAT_STOREY_6_TO_10 = 3
        FLAT_STOREY_11_PLUS = 4

        @classmethod
        def from_string(cls, strval_type, strval_storey):
            if strval_type == "house":
                if strval_storey == 1:
                    return cls.HOUSE_1_STOREY
                elif strval_storey >= 2:
                    return cls.HOUSE_2_STOREY
            elif strval_type == "flat":
                if 0 < strval_storey <= 5:
                    return cls.FLAT_STOREY_1_TO_5
                elif 5 < strval_storey <= 10:
                    return cls.FLAT_STOREY_6_TO_10
                elif strval_storey > 10:
                    return cls.FLAT_STOREY_11_PLUS
            else:
                sys.exit('DwellingType (' + str(strval_type) + ') not valid.')
                # TODO Exit just the current case instead of whole program entirely?

    # Divisors to convert air change rate at 50 Pa to infiltration
    # Values for "Normal" House 1-2 storey and Flat storeys 1-10 are from CIBSE Guide A
    # Values for "Exposed" based on CIBSE Guida A: "on severely exposed sites, a
    # 50% increase to the tabulated values should be allowed.
    # Values for "Sheltered" based on CIBSE Guide A: "on sheltered sites, the
    # infiltration rate may be reduced by 33%".
    # Values for "Very sheltered" assume a reduction of 50%.
    # Values for Flat storeys 11+ are extrapolated based on profiles of wind
    # speed vs. height, assuming storey height of 3.5 metres.
    __DIVISORS = [
        # Very sheltered | Sheltered | Normal | Exposed |
        [           41.2 ,      30.7 ,   20.6 ,    13.7 ], # House 1-storey
        [           34.0 ,      25.4 ,   17.0 ,    11.3 ], # House 2-storey
        [           34.6 ,      25.8 ,   17.3 ,    11.5 ], # Flat storeys 1-5
        [           30.2 ,      22.5 ,   15.1 ,    10.1 ], # Flat storeys 6-10
        [           29.3 ,      19.9 ,   13.7 ,     9.3 ], # Flat storeys 11+
        ]

    def __init__(self,
            storeys_in_building,
            shelter,
            build_type,
            pressure_test_result_ach,
            pressure_test_type,
            env_area,
            volume,
            sheltered_sides,
            open_chimneys,
            open_flues,
            closed_fire,
            flues_d,
            flues_e,
            blocked_chimneys,
            extract_fans,
            passive_vents,
            gas_fires,
            ext_cond,
            storey_of_dwelling = None,
            ):
        """ Construct a VentilationElementInfiltration object """

        """Arguments:
        storeys_in_building   -- total number of storeys in building
        shelter               -- exposure level of the building i.e. very sheltered, sheltered, normal, or exposed
        build_type            -- type of building e.g. house, flat, etc.
        pressure_test_result_ach -- result of pressure test, in ach
        pressure_test_type    -- measurement used for pressure test i.e. based on air change rate value at 50 Pa (50Pa) or 4 Pa (4Pa)
        env_area              -- total envelope area of the building including party walls and floors, in m^2
        volume                -- total volume of dwelling, m^3
        sheltered_sides       -- number of sides of the building which are sheltered
        open_chimneys         -- number of open chimneys
        open_flues            -- number of open flues
        closed_fire           -- number of chimneys / flues attched to closed fire
        flues_d               -- number of flues attached to soild fuel boiler
        flues_e               -- number of flues attached to other heater
        blocked_chimneys      -- number of blocked chimneys
        extract_fans          -- number of intermittent extract fans
        passive_vents         -- number of passive vents
        gas_fires             -- number of flueless gas fires
        ext_cond              -- reference to ExternalConditions object
        storey_of_dwelling    -- for flats only, storey number within building
        """

        self.__external_conditions = ext_cond
        self.__volume = volume
        sheltered_sides = int(sheltered_sides)

        # Calculate the infiltration rate from openings (chimneys, flues, fans, PSVs, etc.)
        # in air-changes per hour
        self.__inf_openings \
            = ( (open_chimneys    * self.__INF_RATE_CHIMNEY_OPEN) \
              + (open_flues       * self.__INF_RATE_FLUE_OPEN) \
              + (closed_fire      * self.__INF_RATE_FIRE_CLOSED) \
              + (flues_d          * self.__INF_RATE_FLUE_SOLID_FUEL_BOILER) \
              + (flues_e          * self.__INF_RATE_FLUE_OTHER_HEATER) \
              + (blocked_chimneys * self.__INF_RATE_CHIMNEY_BLOCKED) \
              + (extract_fans     * self.__INF_RATE_EXTRACT_FAN) \
              + (passive_vents    * self.__INF_RATE_PASSIVE_STACK_VENT)\
              + (gas_fires        * self.__INF_RATE_FIRE_GAS)
              ) \
            / volume

        if build_type == 'flat':
            storey = storey_of_dwelling
        elif build_type == 'house':
            storey = storeys_in_building
        else:
            sys.exit('Error: build type not applicable')

        # Choose correct divisor to apply to Q50:
        # TODO add options for bungalow and maisonette
        def init_divisor():
            dwelling_type = self.DwellingType.from_string(build_type, storey)
            shelter_type = self.ShelterType.from_string(shelter)
            return self.__DIVISORS[dwelling_type][shelter_type]

        self.__divisor = init_divisor()

        # Calculate shelter factor
        def init_shelter_factor():
            # TODO check shelter correction - option 1
            if sheltered_sides < 0 or sheltered_sides > 4:
                sys.exit( ' Number of sheltered sides not recognised.' )
                # TODO Exit just the current case instead of whole program entirely?
            # Calculate shelter factor based on formula from SAP 10.2
            # TODO Reference the origin of this formula.
            return 1.0 - (0.075 * sheltered_sides)

        self.__shelter_factor = init_shelter_factor()

        # Calculate infiltration rate
        def init_infiltration():
            if pressure_test_type == "4Pa":
                # If test results are at 4 Pa, convert to equivalent 50 Pa result
                # before applying divisor.
                # SAP 10 Technical Paper S10TP-19 "Use of low pressure pulse
                # test data in SAP" gives the relationship between air
                # permeability measured at 50 Pa and 4 Pa. The equation below is
                # based on this but has been converted to work with test results
                # expressed in ach rather than m3/m2/h.
                test_result_ach_50Pa \
                    = 5.254 * (pressure_test_result_ach**0.9241) * ((env_area / volume)**(1-0.9241))
            elif pressure_test_type == "50Pa":
                test_result_ach_50Pa = pressure_test_result_ach
            else:
                sys.exit( ' Pressure test result type not recognised.' )
                # TODO Exit just the current case instead of whole program entirely?

            return (test_result_ach_50Pa / self.__divisor) \
                 + (self.__inf_openings * self.__shelter_factor)

        self.__infiltration = init_infiltration()

    def infiltration(self):
        return self.__infiltration

    def h_ve(self, zone_volume):
        """ Calculate the heat transfer coefficient (h_ve), in W/K,
        according to ISO 52016-1:2017, Section 6.5.10.1
        
        Arguments:
        zone_volume -- volume of zone, in m3
        """

        # Apply wind speed correction factor
        wind_factor = self.__external_conditions.wind_speed() / 4.0 # 4.0 m/s represents the average wind speed
        inf_rate = self.__infiltration * wind_factor

        # Convert infiltration rate from ach to m^3/s
        q_v = inf_rate * zone_volume / seconds_per_hour

        # Calculate h_ve according to BS EN ISO 52016-1:2017 section 6.5.10 equation 61
        h_ve = p_a * c_a * q_v
        return h_ve
        # TODO b_ztu needs to be applied in the case if ventilation element
        #      is adjacent to a thermally unconditioned zone.

    def h_ve_average(self, zone_volume):
        """ Calculate the heat transfer coefficient (h_ve), in W/K,
        according to ISO 52016-1:2017, Section 6.5.10.1, for a constant average windspeed
        
        Arguments:
        zone_volume -- volume of zone, in m3
        """

        # Apply wind speed correction factor
        wind_factor = self.__external_conditions.wind_speed_annual() / 4.0  # using average annual wind speed
        inf_rate = self.__infiltration * wind_factor

        # Convert infiltration rate from ach to m^3/s
        q_v = inf_rate * zone_volume / seconds_per_hour

        # Calculate h_ve according to BS EN ISO 52016-1:2017 section 6.5.10 equation 61
        h_ve_average = p_a * c_a * q_v
        return h_ve_average
        # TODO b_ztu needs to be applied in the case if ventilation element
        #      is adjacent to a thermally unconditioned zone.

    def temp_supply(self):
        """ Calculate the supply temperature of the air flow element
        according to ISO 52016-1:2017, Section 6.5.10.2 """
        return self.__external_conditions.air_temp()
        # TODO For now, this only handles ventilation elements to the outdoor
        #      environment, not e.g. elements to adjacent zones.


class MechnicalVentilationHeatRecovery:
    """ A class to represent ventilation with heat recovery (MVHR) elements """

    def __init__(
            self,
            required_air_change_rate,
            specific_fan_power,
            efficiency_hr,
            energy_supply_conn,
            ext_con,
            simulation_time,
            ):
        """ Construct a MechnicalVentilationHeatRecovery object

        Arguments:
        required_air_change_rate -- ach (l/s)m-3 calculated according to Part F
        efficiency_hr -- heat recovery efficiency (0 to 1) allowing for in-use factor
        ext_con -- reference to ExternalConditions object
        """
        self.__air_change_rate = required_air_change_rate
        self.__sfp = specific_fan_power
        self.__efficiency = efficiency_hr
        self.__energy_supply_conn = energy_supply_conn
        self.__external_conditions = ext_con
        self.__simtime = simulation_time

    def h_ve(self, zone_volume, throughput_factor=1.0):
        """ Calculate the heat transfer coefficient (h_ve), in W/K,
        according to ISO 52016-1:2017, Section 6.5.10.1

        Arguments:
        zone_volume -- volume of zone, in m3
        throughput_factor -- proportional increase in ventilation rate due to
                             overventilation requirement
        """

        q_v = air_change_rate_to_flow_rate(self.__air_change_rate, zone_volume) * throughput_factor

        # Calculate effective flow rate of external air
        # NOTE: Technically, the MVHR system supplies air at a higher temperature
        # than the outside air. However, it is simpler to adjust the heat
        # transfer coefficient h_ve to account for the heat recovery effect
        # using an "equivalent" or "effective" flow rate of external air.
        q_v_effective = q_v * (1 - self.__efficiency)

        # Calculate h_ve according to BS EN ISO 52016-1:2017 section 6.5.10 equation 61
        h_ve = p_a * c_a * q_v_effective
        return h_ve
        # TODO b_ztu needs to be applied in the case if ventilation element
        #      is adjacent to a thermally unconditioned zone.

    def h_ve_average(self, zone_volume):
        return self.h_ve(zone_volume)
        # TODO b_ztu needs to be applied in the case if ventilation element
        #      is adjacent to a thermally unconditioned zone.

    def fans(self, zone_volume, throughput_factor=1.0):
        """ Calculate gains and energy use due to fans"""
        # Calculate energy use by fans (only fans on intake/supply side
        # contribute to internal gains - assume that this is half of the fan
        # power)
        q_v = air_change_rate_to_flow_rate(self.__air_change_rate, zone_volume) \
            * throughput_factor
        fan_power_W = self.__sfp * (q_v * litres_per_cubic_metre)
        fan_energy_use_kWh = (fan_power_W  / W_per_kW) * self.__simtime.timestep()

        self.__energy_supply_conn.demand_energy(fan_energy_use_kWh)
        return fan_energy_use_kWh / 2.0

    def temp_supply(self):
        """ Calculate the supply temperature of the air flow element
        according to ISO 52016-1:2017, Section 6.5.10.2 """
        # NOTE: Technically, the MVHR system supplies air at a higher temperature
        # than the outside air, i.e.:
        #     temp_supply = self.__efficiency * temp_int_air \
        #                 + (1 - self.__efficiency) * self.__external_conditions.air_temp()
        # However, calculating this requires the internal air temperature, which
        # has not been calculated yet. Calculating this properly would require
        # the equation above to be added to the heat balance solver. Therefore,
        # it is simpler to adjust the heat transfer coefficient h_ve to account
        # for the heat recovery effect using an "equivalent" flow rate of
        # external air.
        return self.__external_conditions.air_temp()
        
    def efficiency(self):
        return self.__efficiency


class WholeHouseExtractVentilation:
    """ A class to represent whole house extract ventilation elements """

    def __init__(
            self,
            required_air_change_rate,
            specific_fan_power,
            infiltration_rate,
            energy_supply_conn,
            ext_con,
            simulation_time
            ):
        """ Construct a WholeHouseExtractVentilation object

        Arguments:
        required_air_change_rate -- in ach
        specific_fan_power -- in W / (litre / second), inclusive of any in-use factors
        infiltration_rate -- in ach, not adjusted for wind speed
        energy_supply_conn -- reference to EnergySupplyConnection object
        ext_con -- reference to ExternalConditions object
        """

        self.__air_change_rate_req = required_air_change_rate
        self.__infiltration_rate = infiltration_rate
        self.__sfp = specific_fan_power
        self.__energy_supply_conn = energy_supply_conn
        self.__external_conditions = ext_con
        self.__simtime = simulation_time

    def air_change_rate(self, infiltration_rate):
        """ Calculate air change rate for the system

        Arguments:
        infiltration_rate -- in ach, adjusted for wind speed
        """
        # The calculation below is based on SAP 10.2 equations, but with the
        # sharp "elbow" at infiltration_rate == 0.5 * air_change_rate_req
        # replaced by a smooth curve between infiltration_rate == 0 and
        # infiltration_rate = air_change_rate_req.
        # As we are already accounting for infiltration separately, it is
        # subtracted from the totals here, compared to the equations in SAP 10.2
        if infiltration_rate < self.__air_change_rate_req:
            ach = self.__air_change_rate_req - infiltration_rate \
                + (infiltration_rate ** 2 * 0.5 / self.__air_change_rate_req)
        else:
            ach = 0.5 * self.__air_change_rate_req
        return ach

    def h_ve(self, zone_volume, throughput_factor=1.0):
        """ Calculate the heat transfer coefficient (h_ve), in W/K,
        according to ISO 52016-1:2017, Section 6.5.10.1

        Arguments:
        zone_volume -- volume of zone, in m3
        throughput_factor -- proportional increase in ventilation rate due to
                             overventilation requirement
        """

        infiltration_rate_adj \
            = self.__infiltration_rate * self.__external_conditions.wind_speed() / 4.0
        ach = self.air_change_rate(infiltration_rate_adj)
        q_v = air_change_rate_to_flow_rate(ach, zone_volume) * throughput_factor

        # Calculate h_ve according to BS EN ISO 52016-1:2017 section 6.5.10 equation 61
        h_ve = p_a * c_a * q_v
        return h_ve
        # TODO b_ztu needs to be applied in the case if ventilation element
        #      is adjacent to a thermally unconditioned zone.

    def h_ve_average(self, zone_volume):
        """ Calculate the heat transfer coefficient (h_ve), in W/K,
        according to ISO 52016-1:2017, Section 6.5.10.1 using an average windspeed of 4 m/s

        Arguments:
        zone_volume -- volume of zone, in m3
        """

        infiltration_rate_adj \
            = self.__infiltration_rate * self.__external_conditions.wind_speed_annual() / 4.0  # using average annual wind speed
        ach = self.air_change_rate(infiltration_rate_adj)
        q_v = air_change_rate_to_flow_rate(ach, zone_volume)

        # Calculate h_ve according to BS EN ISO 52016-1:2017 section 6.5.10 equation 61
        h_ve_average = p_a * c_a * q_v
        return h_ve_average
        # TODO b_ztu needs to be applied in the case if ventilation element
        #      is adjacent to a thermally unconditioned zone.

    def fans(self, zone_volume, throughput_factor=1.0):
        """ Calculate gains and energy use due to fans """
        # Calculate energy use by fans (does not contribute to internal gains as
        # this is extract-only ventilation)
        q_v = air_change_rate_to_flow_rate(self.__air_change_rate_req, zone_volume) \
            * throughput_factor
        fan_power_W = self.__sfp * (q_v * litres_per_cubic_metre)
        fan_energy_use_kWh = (fan_power_W  / W_per_kW) * self.__simtime.timestep()

        self.__energy_supply_conn.demand_energy(fan_energy_use_kWh)
        return 0.0

    def temp_supply(self):
        """ Calculate the supply temperature of the air flow element
        according to ISO 52016-1:2017, Section 6.5.10.2 """
        return self.__external_conditions.air_temp()
        # TODO For now, this only handles ventilation elements to the outdoor
        #      environment, not e.g. elements to adjacent zones.


class NaturalVentilation:
    """ A class to represent natural ventilation """

    def __init__(
            self,
            required_air_change_rate,
            infiltration_rate,
            ext_con,
            ):
        """ Construct a NaturalVentilation object

        Arguments:
        required_air_change_rate -- in ach
        infiltration_rate -- in ach, not adjusted for wind speed
        ext_con -- reference to ExternalConditions object
        """
        self.__air_change_rate_req = required_air_change_rate
        self.__infiltration_rate = infiltration_rate
        self.__external_conditions = ext_con

    def air_change_rate(self, infiltration_rate):
        """ Calculate air change rate for the system

        Arguments:
        infiltration_rate -- in ach, adjusted for wind speed
        """
        # The calculation below is based on SAP 10.2 equations, but with the curve between
        # infiltration_rate == 0 and infiltration_rate == 2 * air_change_rate_req
        # adjusted to handle values of air_change_rate_req other than 0.5.
        # As we are already accounting for infiltration separately, it is
        # subtracted from the totals here, compared to the equations in SAP 10.2
        if infiltration_rate < 2.0 * self.__air_change_rate_req:
            ach = self.__air_change_rate_req - infiltration_rate \
                + (infiltration_rate ** 2 * 0.25 / self.__air_change_rate_req)
        else:
            # No additional ventilation (infiltration already accounted for separately)
            ach = 0.0
        return ach

    def h_ve(self, zone_volume):
        """ Calculate the heat transfer coefficient (h_ve), in W/K,
        according to ISO 52016-1:2017, Section 6.5.10.1

        Arguments:
        zone_volume -- volume of zone, in m3
        inf_rate -- air change rate of ventilation system
        """

        infiltration_rate_adj \
            = self.__infiltration_rate * self.__external_conditions.wind_speed() / 4.0
        ach = self.air_change_rate(infiltration_rate_adj)
        q_v = air_change_rate_to_flow_rate(ach, zone_volume)

        # Calculate h_ve according to BS EN ISO 52016-1:2017 section 6.5.10 equation 61
        h_ve = p_a * c_a * q_v

        return h_ve
        # TODO b_ztu needs to be applied in the case if ventilation element
        #      is adjacent to a thermally unconditioned zone.

    def h_ve_average(self, zone_volume):
        """ Calculate the heat transfer coefficient (h_ve), in W/K,
        according to ISO 52016-1:2017, Section 6.5.10.1

        Arguments:
        zone_volume -- volume of zone, in m3
        """

        infiltration_rate_adj \
            = self.__infiltration_rate * self.__external_conditions.wind_speed_annual() / 4.0  # using average annual wind speed
        ach = self.air_change_rate(infiltration_rate_adj)
        q_v = air_change_rate_to_flow_rate(ach, zone_volume)

        # Calculate h_ve according to BS EN ISO 52016-1:2017 section 6.5.10 equation 61
        h_ve_average = p_a * c_a * q_v
        return h_ve_average
        # TODO b_ztu needs to be applied in the case if ventilation element
        #      is adjacent to a thermally unconditioned zone.

    def temp_supply(self):
        """ Calculate the supply temperature of the air flow element
        according to ISO 52016-1:2017, Section 6.5.10.2 """
        return self.__external_conditions.air_temp()
        # TODO For now, this only handles ventilation elements to the outdoor
        #      environment, not e.g. elements to adjacent zones.


class WindowOpeningForCooling:
    """
    Object to represent window opening for providing additional ventilation
    in response to high internal temperature.

    This is based on the simple building layouts from CIBSE Guide A section 4.6.2
    (Method 2) and Tables 4.25 and 4.26. First, the layout that best fits the building
    is determined, then the characteristics of the openings are aggregated into
    the parameters needed by the equations.
    """

    def __init__(
            self,
            window_area_equivalent,
            external_conditions,
            openings,
            control,
            natvent = None,
            ):
        """ Construct a WindowOpeningForCooling object

        Arguments:
        window_area_equivalent -- maximum equivalent area of all openings in the relevant zone
        external_conditions -- reference to ExternalConditions object
        openings -- list of openings to be considered
        control -- reference to control object (must implement setpnt function)
        natvent -- reference to NaturalVentilation object, if building is naturally ventilated
        """
        self.__window_area_equiv = window_area_equivalent
        self.__external_conditions = external_conditions
        self.__control = control
        self.__natvent = natvent

        # Assign equivalent areas to each window/group in proportion to actual area
        opening_area_total = sum(op.area for op in openings)
        opening_area_equiv_total_ratio = window_area_equivalent / opening_area_total

        # Find orientation of largest window
        # Find height of highest and lowest windows
        largest_op = openings[0]
        highest_op = openings[0]
        lowest_op = openings[0]
        for op in openings[1:]:
            # TODO What do we do if two openings are joint-largest? Calculate for
            #      both (or all windows) and take whichever gives largest areas?
            if op.area > largest_op.area:
                largest_op = op
            if op.mid_height() > highest_op.mid_height():
                highest_op = op
            if op.mid_height() < lowest_op.mid_height():
                lowest_op = op
        largest_op_orientation = largest_op.orientation()
        op_height_threshold = (highest_op.mid_height() - lowest_op.mid_height()) / 2.0

        openings_same_side = []
        openings_opp_side = []
        openings_low = []
        openings_high = []
        for op in openings:
            # Determine orientation of other windows relative to largest
            op_rel_orientation = abs(op.orientation() - largest_op_orientation)
            if op_rel_orientation > 360:
                op_rel_orientation -= 360
            # Group windows into same, opposite and adjacent sides
            if op_rel_orientation <= 45:
                # Opening is on same side
                openings_same_side.append(op)
            elif op_rel_orientation >= 135:
                # Opening is on opp side
                openings_opp_side.append(op)
            # Else opening is on adjacent side, so ignore

            # Assign windows to high and low groups based on which they are closest to
            if op.mid_height() < op_height_threshold:
                openings_low.append(op)
            else:
                openings_high.append(op)

        # Determine whether cross ventilation is possible
        self.__cross_vent = len(openings_opp_side) > 0

        if self.__cross_vent:
            openings_same_side_high = list(set(openings_same_side).intersection(openings_high))
            openings_same_side_low = list(set(openings_same_side).intersection(openings_low))
            openings_opp_side_high = list(set(openings_opp_side).intersection(openings_high))
            openings_opp_side_low = list(set(openings_opp_side).intersection(openings_low))

            # Calculate high and low opening areas on same and opposite sides of building
            A1 = opening_area_equiv_total_ratio * sum(op.area for op in openings_same_side_high)
            A2 = opening_area_equiv_total_ratio * sum(op.area for op in openings_same_side_low)
            A3 = opening_area_equiv_total_ratio * sum(op.area for op in openings_opp_side_high)
            A4 = opening_area_equiv_total_ratio * sum(op.area for op in openings_opp_side_low)
            self.__A_w = sqrt(1.0 / ( (1.0 / ((A1 + A2) ** 2)) + 1.0 / ((A3 + A4) ** 2) ))
            if A2 + A4 == 0.0:
                self.__A_b = 0
                self.__opening_height_diff = 0.0
            else:
                self.__A_b = sqrt(1.0 / ( (1.0 / ((A1 + A3) ** 2)) + 1.0 / ((A2 + A4) ** 2) ))

                # Calculate area-weighted average height of windows in high and low groups
                opening_mid_height_ave_upper \
                    = ( sum(op.mid_height() * op.area for op in openings_same_side_high)
                      + sum(op.mid_height() * op.area for op in openings_opp_side_high)
                      ) \
                    / ( sum(op.area for op in openings_same_side_high)
                      + sum(op.area for op in openings_opp_side_high)
                      )
                opening_mid_height_ave_lower \
                    = ( sum(op.mid_height() * op.area for op in openings_same_side_low)
                      + sum(op.mid_height() * op.area for op in openings_opp_side_low)
                      ) \
                    / ( sum(op.area for op in openings_same_side_low)
                      + sum(op.area for op in openings_opp_side_low)
                      )
                # Calculate opening height difference
                self.__opening_height_diff \
                    = opening_mid_height_ave_upper - opening_mid_height_ave_lower
        else:
            # Determine whether stack ventilation is possible
            if len(openings_high) > 1 and len(openings_low) > 1:
                self.__stack_vent = True
                opening_area_upper = sum(op.area for op in openings_high)
                opening_area_lower = sum(op.area for op in openings_low)

                # Calculate opening area ratio
                self.__opening_area_ratio = opening_area_upper / opening_area_lower

                # Calculate area-weighted average height of windows in high and low groups
                opening_mid_height_ave_upper \
                    = sum(op.mid_height() * op.area for op in openings_high) \
                    / opening_area_upper
                opening_mid_height_ave_lower \
                    = sum(op.mid_height() * op.area for op in openings_low) \
                    / opening_area_lower
                # Calculate opening height difference
                self.__opening_height_diff \
                    = opening_mid_height_ave_upper - opening_mid_height_ave_lower
            else:
                self.__stack_vent = False
                self.__openings = openings

    def temp_setpnt(self):
        return self.__control.setpnt()

    def h_ve_max(self, zone_volume, temp_int):
        g = 9.81 # m/s
        C_d = 0.62 # Discharge coeff
        # Difference in wind pressure coeff from CIBSE Guide A Table 4.12 for urban environment
        # TODO Should we differentiate between urban and rural, perhaps based on
        #      the shelter input to the infiltration calculation?
        dC_p = 0.2 - (-0.25)
        wind_speed = self.__external_conditions.wind_speed()
        temp_ext = self.__external_conditions.air_temp()
        temp_diff = abs(temp_int - temp_ext)
        temp_average_C = (temp_int + temp_ext) / 2.0
        temp_average_K = Celcius2Kelvin(temp_average_C)

        if self.__cross_vent:
            q_v_wind = C_d * self.__A_w * wind_speed * dC_p ** 0.5
            q_v_stack \
                = C_d * self.__A_b \
                * ((2.0 * temp_diff * self.__opening_height_diff * g) / temp_average_K) ** 0.5
            if wind_speed / sqrt(temp_diff) \
             < 0.26 * (self.__A_b / self.__A_w) * (self.__opening_height_diff * dC_p) ** 0.5:
                q_v = q_v_stack
            else:
                q_v = q_v_wind
        else:
            q_v_wind = 0.025 * self.__window_area_equiv * wind_speed
            if self.__stack_vent:
                q_v_stack = C_d * self.__window_area_equiv \
                          * ( self.__opening_area_ratio * sqrt(2)
                            / ( (1 + self.__opening_area_ratio)
                              * (1 + self.__opening_area_ratio ** 2) ** 0.5
                              )
                            ) \
                          * ( (temp_diff * self.__opening_height_diff * g)
                              / temp_average_K
                            ) \
                            ** 0.5
            else: # Stack effect is only between top and bottom of each opening
                q_v_stack = 0.0
                for op in self.__openings:
                    q_v_stack += C_d * self.__window_area_equiv / 3.0 \
                               * ( (temp_diff * op.projected_height() * g)
                                   / temp_average_K
                                 ) \
                                 ** 0.5
            q_v = max(q_v_wind, q_v_stack)

        # Calculate h_ve according to BS EN ISO 52016-1:2017 section 6.5.10 equation 61
        h_ve = p_a * c_a * q_v

        # Calculate max h_ve achievable for window opening
        if self.__natvent is not None:
            h_ve_nat_vent = self.__natvent.h_ve(zone_volume)
        else:
            h_ve_nat_vent = 0.0

        # Avoid double-counting window opening for ventilation and cooling purposes
        return max(0.0, h_ve - h_ve_nat_vent)

    def temp_supply(self):
        """ Calculate the supply temperature of the air flow element
        according to ISO 52016-1:2017, Section 6.5.10.2 """
        return self.__external_conditions.air_temp()
