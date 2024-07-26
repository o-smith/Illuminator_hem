#!/usr/bin/env python3

"""
This module provides objects to represent Infiltration and Ventilation.
The calculations are based on Method 1 of BS EN 16798-7.
"""

# Standard library imports
import sys
from scipy.optimize import fsolve, root_scalar
import warnings

# Third-party imports
import numpy as np

# Local imports
from core.units import seconds_per_hour, Celcius2Kelvin, litres_per_cubic_metre, \
    W_per_kW
from core.simulation_time import SimulationTime
from core.external_conditions import ExternalConditions
from core.material_properties import AIR
from cmath import pi

# Define constants
p_a_ref = AIR.density_kg_per_m3()
c_a = AIR.specific_heat_capacity_kWh()

# (Default values from BS EN 16798-7, Table 11)
# Coefficient to take into account stack effect in airing calculation in (m/s)/(m*K)
C_stack = 0.0035
# Coefficient to take into account wind speed in airing calculation in 1/(m/s)
C_wnd = 0.001
# Gravitational constant in m/s2
g = 9.81
#Room temperature in degrees K
T_e_ref = 293.15
#Absolute zero in degrees K
T_0_abs = 273.15


def calculate_pressure_difference_at_an_airflow_path(
        h_path,
        C_p_path,
        u_site,
        T_e,
        T_z,
        p_z_ref
        ):
    """Calculate pressure difference between the exterior and the interior of the dwelling
    for a flow path (at it's elevavation above the vent zone floor)

    Arguments:
    h_path -- height of air flow path (m)
    C_p_path -- wind pressure coefficient
    u_site -- wind velocity at zone level (m/s)
    T_e -- external air temperature (K)
    T_z -- thermal zone air temperature (K)
    p_z_ref -- internal reference pressure (Pa)
    """
    p_e_path = p_a_ref * T_e_ref/T_e * (0.5 * C_p_path * u_site**2 - h_path * g) #(5)
    p_z_path = p_z_ref - p_a_ref * h_path * g * T_e_ref/T_z #(6)
    delta_p_path = p_e_path - p_z_path #(4)
    return delta_p_path

def air_change_rate_to_flow_rate(air_change_rate, zone_volume):
    """ Convert infiltration rate from ach to m^3/s """
    return air_change_rate * zone_volume / seconds_per_hour

def get_fuel_flow_factor(fuel_type, appliance_type):
    """ Table B.3 fuel flow factors

    Arguments:
    fuel_type -- options are 'wood', 'gas', 'oil' or 'coal'.
    appliance_type -- options are 'open_fireplace', 'closed_with_fan',
        open_gas_flue_balancer', 'open_gas_kitchen_stove', 'open_gas_fire' or 'closed_fire'
    """
    if fuel_type == 'wood':
        if appliance_type == 'open_fireplace':
            f_ff = 2.8
        else:
            sys.exit("appliance_type: "+ str(appliance_type) +" not applicable for the fuel_type: "+ str(fuel_type) +" selected.")
    elif fuel_type == 'gas':
        if appliance_type == 'closed_with_fan':
            f_ff = 0.38
        elif appliance_type == 'open_gas_flue_balancer':
            f_ff = 0.78
        elif appliance_type == 'open_gas_kitchen_stove':
            f_ff = 3.35
        elif appliance_type == 'open_gas_fire':
            f_ff = 3.35
        else:
            sys.exit("appliance_type: "+ str(appliance_type) +" not applicable for the fuel_type: "+ str(fuel_type) +" selected.")
    elif fuel_type == 'oil':
        if appliance_type == 'closed_fire':
            f_ff = 0.32
        else:
            sys.exit("appliance_type: "+ str(appliance_type) +" not applicable for the fuel_type: "+ str(fuel_type) +" selected.")
    elif fuel_type == 'coal':
        if appliance_type == 'closed_fire':
            f_ff = 0.52
        else:
            sys.exit("appliance_type: "+ str(appliance_type) +" not applicable for the fuel_type: "+ str(fuel_type) +" selected.")
    else:
        sys.exit(" fuel_type: "+ str(fuel_type) +" not found.")
    return f_ff

def get_appliance_system_factor(supply_situation, exhaust_situation):
    """Interpreted from Table B.2 from BS EN 16798-7, get the appliance system factor
    for a combustion appliance.

    Arguments:
    supply_situation -- Combustion air supply situation: 'room_air' or 'outside'
    exhaust_situation -- flue gas exhaust situation: 'into_room', 'into_separate_duct' or 'into_mech_vent'
    """
    if supply_situation == "outside":
        f_as = 0
    elif supply_situation == "room_air":
        if exhaust_situation == "into_room":
            f_as = 0
        elif exhaust_situation == "into_separate_duct":
            f_as = 1
        elif exhaust_situation == "into_mech_vent":
            sys.exit("Cannot currently handle 'exhaust_situation': "+ str(exhaust_situation) +" for combustion appliance.")
        else:
            sys.exit("'exhaust_situation': "+ str(exhaust_situation) +" not recognised for combustion appliance")
    else:
        sys.exit("'supply_situation': "+ str(supply_situation) +" not recognised for combustion appliance")
    return f_as

def adjust_air_density_for_altitude(h_alt):
    """Adjust air density for altitude above sea level.

    Arguments:
    h_alt -- altitude above sea level (m)
    """
    p_a_alt = p_a_ref * (1 - ((0.00651 * h_alt) / 293))** 4.255
    return p_a_alt

def air_density_at_temp(temperature, air_density_adjusted_for_alt):
    """Recalculate air density based on the current temperature

    Arguments:
    temperature -- temperature to adjust (K)
    air_density_adjusted_for_alt - The air density after adjusting for altitude (Kg/m3)
    """
    return T_e_ref / temperature * air_density_adjusted_for_alt

def convert_to_mass_air_flow_rate(qv_in, qv_out, T_e, T_z, p_a_alt):
    """Converts volume air flow rate (qv) to mass air flow rate (qm).
    (Equations 65 & 66 from BS EN 16798-7)

    Arguments:
    qv_in -- volume flow rate of air entering the dwelling
    qv_out -- volume flow rate of air leaving the dwelling
    T_e -- External air temperature (K)
    T_e -- Thermal zone air temperature (K)
    p_a_alt -- The air density after adjusting for altitude (Kg/m3)
    """
    qm_in = convert_volume_flow_rate_to_mass_flow_rate(qv_in, T_e, p_a_alt)
    qm_out = convert_volume_flow_rate_to_mass_flow_rate(qv_out, T_z, p_a_alt)
    return qm_in, qm_out

def convert_volume_flow_rate_to_mass_flow_rate(qv, temperature, p_a_alt):
    """ Convert volume flow rate in m3/hr to mass flow rate in kg/hr, at temperature in Kelvin 

    Arguments:
    qv -- volume flow rate (m3/h)
    temperature -- air temperature (K)
    p_a_alt -- The air density after adjusting for altitude (Kg/m3)
    """
    return qv * air_density_at_temp(temperature, p_a_alt)

def convert_mass_flow_rate_to_volume_flow_rate(qm, temperature, p_a_alt):
    """ Convert mass flow rate in kg/hr to volume flow rate in m3/hr, at temperature in Kelvin

    Arguments:
    qm -- mass flow rate (Kg/h)
    temperature -- air temperature (K)
    p_a_alt -- The air density after adjusting for altitude (Kg/m3)
    """
    return qm / air_density_at_temp(temperature, p_a_alt)

def ter_class_to_roughness_coeff(TER_CLASS):
    """Interpreted from Table B.13 in BS EN 16798-7.
    Terrain Class input to roughness coefficient at building site at 10m
    Arguments:
    TER_CLASS -- Terrain class, one of 'Open terrain', 'Country' or 'Urban'
    """
    if TER_CLASS == "Open terrain":
        C_rgh_10_site = 1.0
    elif TER_CLASS == "Country":
        C_rgh_10_site = 0.9
    elif TER_CLASS == "Urban":
        C_rgh_10_site = 0.8
    return C_rgh_10_site

def wind_speed_at_zone_level(C_rgh_site, u_10, C_top_site=1, C_rgh_met=1, C_top_met=1):
    """Meteorological wind speed at 10 m corrected to reference wind speed at zone level of the dwelling

    Arguments:
    C_rgh_site -- roughness coefficient at building site
    u_10 -- wind velocity at 10m (m/s)
    C_top_site -- topography coefficient at building site
    C_rgh_met -- roughness coefficient at 10m depending on meteorological station
    C_top_met -- topography coefficient at building height depending on meteorological station
    """
    u_site = ((C_rgh_site * C_top_site) / (C_rgh_met * C_top_met)) * u_10
    return u_site

#Not needed as using internal pressure for windows, not cross ventilation
# def wind_speed_at_10m_height(C_rgh_10_site, u_10, C_top_10_site=1, C_rgh_met=1, C_top_met=1):
#     """The meteorological wind speed at 10 m is corrected as follows to obtain the reference wind speed at
#     site at 10m height:
#     """
#     u_10_site = ((C_rgh_10_site * C_top_10_site) / (C_rgh_met * C_top_met)) * u_10
#     return u_10_site

def orientation_difference(orientation1, orientation2):
    """Determine orientation of other windows relative to largest"""
    op_rel_orientation = abs(orientation1 - orientation2)
    if op_rel_orientation > 360:
        op_rel_orientation -= 360
    return op_rel_orientation

def get_facade_direction(f_cross, orientation, pitch, wind_direction):
    """Gets direction of the facade from pitch and orientation

    Arguments:
    f_cross -- boolean, dependant on if cross ventilation is possible or not
    orientation -- orientation of the facade (degrees)
    pitch -- pitch of the facade (degrees)
    wind_direction -- direction the wind is blowing (degrees)
    """
    if f_cross == True:
        if pitch < 10:
            facade_direction = "Roof10"
        elif pitch <= 30:
            facade_direction = "Roof10_30"
        elif pitch < 60:
            facade_direction = "Roof30"
        else:
            orientation_diff = orientation_difference(orientation, wind_direction)
            if orientation_diff < 90:
                facade_direction = "Windward"
            else:
                facade_direction = "Leeward"
    elif f_cross == False:
        if pitch < 60:
            facade_direction = "Roof"
        else:
            orientation_diff = orientation_difference(orientation, wind_direction)
            if orientation_diff < 90:
                facade_direction = "Windward"
            else:
                facade_direction = "Leeward"
    else:
        sys.exit('f_cross should be set to True or False')
    return facade_direction

def get_C_p_path(f_cross,
                 shield_class,
                 h_path,
                 wind_direction,
                 orientation=None,
                 pitch=None,
                 facade_direction=None,
                 ):
    """Interpreted from Table B.7 for determining dimensionless wind pressure coefficients

    Arguments:
    f_cross -- boolean, dependant on if cross ventilation is possible or not
    shield_class -- indicates exposure to wind
    h_path - height of flow path (m)
    wind_direction -- direction the wind is blowing (degrees)
    orientation -- orientation of the facade (degrees)
    pitch -- pitch of the facade (degrees)
    facade_direction -- direction of the facade (from get_facade_direction or manual entry)
    """
    if facade_direction == None:
        if orientation is None or pitch is None:
            sys.exit('You have not entered a pitch or orientation for an opaque building element.')
        else:
            facade_direction = get_facade_direction(f_cross, orientation, pitch, wind_direction)
    if f_cross == True:
        if h_path < 15:
            if shield_class == "Open":
                C_p = {"Windward": 0.50, "Leeward": -0.70, "Roof10": -0.70, "Roof10_30": -0.60, "Roof30": -0.20}
            elif shield_class == "Normal":
                C_p = {"Windward": 0.25, "Leeward": -0.50, "Roof10": -0.60, "Roof10_30": -0.50, "Roof30": -0.20}
            elif shield_class == "Shielded":
                C_p = {"Windward": 0.05, "Leeward": -0.30, "Roof10": -0.50, "Roof10_30": -0.40, "Roof30": -0.20}
        elif 15 <= h_path < 50:
            if shield_class == "Open":
                C_p = {"Windward": 0.65, "Leeward": -0.70, "Roof10": -0.70, "Roof10_30": -0.60, "Roof30": -0.20}
            elif shield_class == "Normal":
                C_p = {"Windward": 0.45, "Leeward": -0.50, "Roof10": -0.60, "Roof10_30": -0.50, "Roof30": -0.20}
            elif shield_class == "Shielded":
                C_p = {"Windward": 0.25, "Leeward": -0.30, "Roof10": -0.50, "Roof10_30": -0.40, "Roof30": -0.20}
        elif h_path >= 50:
            if shield_class == "Open":
                C_p = {"Windward": 0.80, "Leeward": -0.70, "Roof10": -0.70, "Roof10_30": -0.60, "Roof30": -0.20}
    elif f_cross == False:
        C_p = {"Windward": 0.05, "Leeward": -0.05, "Roof": 0}
    else:
        sys.exit('f_cross should be set to True or False')
    return C_p[facade_direction]


class Window:
    """ An object to represent Windows """
    def __init__(
            self,
            extcond,
            h_w_fa,
            h_w_path,
            A_w_max,
            window_part_list,
            orientation,
            pitch,
            altitude,
            on_off_ctrl_obj,
            ):
        """Construct a Window object

        Arguments:
            extcond -- The external conditions.
            h_w_fa -- The free area height of the window
            h_w_path -- The midheight of the window.
            A_w_max -- The maximum window opening area.
            window_part_list -- The list of window parts.
            orientation -- The orientation of the window.
            pitch -- The pitch of the window.
            altitude -- altitude of dwelling above sea level (m)
            on_off_ctrl_obj - 

        Method
            - Based on Section 6.4.3.5 Airflow due to windows opening section.
        """
        self.__h_w_fa = h_w_fa
        self.__h_w_path = h_w_path
        self.__A_w_max = A_w_max
        self.__C_D_w = 0.67 # Discharge coefficient for windows based on Section B.3.2.1 of BS EN 16798-7:2017
        self.__n_w = 0.5 # Flow exponent of window based on Section B.3.2.2 of BS EN 16798-7:2017
        self.__orientation = orientation 
        self.__pitch = pitch
        self.__external_conditions = extcond
        self.__N_w_div = max(len(window_part_list) - 1 , 0) 
        self.__on_off_ctrl_obj = on_off_ctrl_obj
        self.__altitude = altitude
        self.__p_a_alt = adjust_air_density_for_altitude(altitude)

        self.__window_parts = []
        for window_part_number, window_part in enumerate(window_part_list):
            self.__window_parts.append(
                WindowPart(
                    window_part["mid_height_air_flow_path"], #h_w_path
                    self.__h_w_fa,
                    self.__N_w_div,
                    window_part_number + 1,
                    )
                )

    def calculate_window_opening_free_area(self, R_w_arg):
        """The window opening free area A_w for a window 
        Equation 40 in BS EN 16798-7.
        Arguments:
        R_w_arg -- ratio of window opening (0-1)
        """
        # Assume windows are shut if the control object is empty      
        if self.__on_off_ctrl_obj is None:
            R_w_arg = 0
        if self.__on_off_ctrl_obj is not None and self.__on_off_ctrl_obj.is_on() == False:
            R_w_arg = 0        
        return R_w_arg * self.__A_w_max #self.__A_w

    def calculate_flow_coeff_for_window(self, R_w_arg):
        """The C_w_path flow coefficient for a window
        Equation 54 from BS EN 16798-7
        Arguments:
        R_w_arg -- ratio of window opening (0-1)
        """
        # Assume windows are shut if the control object is empty      
        if self.__on_off_ctrl_obj is None:
            R_w_arg = 0
        if self.__on_off_ctrl_obj is not None and self.__on_off_ctrl_obj.is_on() == False:
            R_w_arg = 0
        A_w = self.calculate_window_opening_free_area(R_w_arg)
        return (3600 * self.__C_D_w * A_w * (2/p_a_ref)** self.__n_w)

    def calculate_flow_from_internal_p(self, u_site, T_z, p_z_ref, f_cross, shield_class, R_w_arg = None):
        """Calculate the airflow through window opening based on the how open the window is and internal pressure

        Arguments:
        u_site -- wind velocity at zone level (m/s)
        T_z -- thermal zone air temperature (K)
        p_z_ref -- internal reference pressure (Pa)
        f_cross -- boolean, dependant on if cross ventilation is possible or not
        shield_class -- indicates exposure to wind
        R_w_arg -- ratio of window opening (0-1)
        """
        wind_direction = self.__external_conditions.wind_direction()
        T_e = Celcius2Kelvin(self.__external_conditions.air_temp())   
        # Assume windows are shut if the control object is empty      
        if self.__on_off_ctrl_obj is None:
            R_w_arg = 0
        if self.__on_off_ctrl_obj is not None and self.__on_off_ctrl_obj.is_on() == False:
            R_w_arg = 0

        # Wind pressure coefficient for the window
        C_p_path = get_C_p_path(
            f_cross,
            shield_class,
            self.__h_w_path,
            wind_direction,
            self.__orientation,
            self.__pitch,
            )

        # Airflow coefficient of the window
        C_w_path = self.calculate_flow_coeff_for_window(R_w_arg)

        # Sum airflow through each window part entering and leaving - based on Equation 56 and 57
        qv_in_through_window_opening = 0
        qv_out_through_window_opening = 0
        for window_part in self.__window_parts:
            air_flow = window_part.calculate_ventilation_through_windows_using_internal_p(
                u_site,
                T_e,
                T_z,
                C_w_path,
                p_z_ref,
                C_p_path,
                )
            if air_flow >= 0:
                qv_in_through_window_opening += air_flow
            else:
                qv_out_through_window_opening += air_flow
        
        # Convert volume air flow rate to mass air flow rate
        qm_in_through_window_opening, qm_out_through_window_opening = convert_to_mass_air_flow_rate(
            qv_in_through_window_opening,
            qv_out_through_window_opening,
            T_e,
            T_z,
            self.__p_a_alt
        )
        return qm_in_through_window_opening, qm_out_through_window_opening


class WindowPart:
    def __init__(
            self,
            h_w_path,
            h_w_fa,
            N_w_div,
            window_part_number,
            ):
        """
        Construct a WindowPart object

        Argument:
            h_w_path -- Mid-height of the window
            h_w_fa -- free area height of the window
            N_w_div -- number of window divisions
            window_part_number -- The identifying number of the window part
        """
        self.__h_w_path = h_w_path #Find from inp file - Midheight of window part
        self.__h_w_fa = h_w_fa #Find from inp file
        self.__N_w_div = N_w_div
        self.__h_w_div_path = self.calculate_height_for_delta_p_w_div_path(window_part_number)
        self.__n_w = 0.5 # Flow exponent of window

    def calculate_ventilation_through_windows_using_internal_p(
            self,
            u_site,
            T_e,
            T_z,
            C_w_path,
            p_z_ref,
            C_p_path,
            ):
        """Calculate the airflow through window parts from internal pressure

        Arguments:
        u_site -- wind velocity at zone level (m/s)
        T_e -- external air temperature (K)    
        T_z -- thermal zone air temperature (K)
        C_w_path -- wind pressure coefficient at height of the window
        p_z_ref -- internal reference pressure (Pa)
        C_p_path -- wind pressure coefficient at the height of the window part
        """
        delta_p_path = calculate_pressure_difference_at_an_airflow_path(
            self.__h_w_div_path,
            C_p_path,
            u_site,
            T_e,
            T_z,
            p_z_ref
        )
        # Based on Equation 53
        qv_w_div_path = C_w_path / (self.__N_w_div + 1) \
                      * np.sign(delta_p_path) \
                      * abs(delta_p_path) ** self.__n_w
        return qv_w_div_path

    def calculate_height_for_delta_p_w_div_path(self, j):
        """ The height to be considered for delta_p_w_div_path
        Equation 55 from BS EN 16798-7"""
        h_w_div_path = self.__h_w_path - self.__h_w_fa/2 + self.__h_w_fa/(2 * (self.__N_w_div + 1)) +\
                      (self.__h_w_fa/(self.__N_w_div + 1)) * (j-1)
        return h_w_div_path


class Vent:
    """ An object to represent Vents """
    def __init__(
            self,
            extcond,
            h_path,
            A_vent,
            delta_p_vent_ref,
            orientation,
            pitch,
            altitude,
            ):
        """ Construct a Vent object

        Arguments:
            external_conditions -- reference to ExternalConditions object
            h_path -- mid height of air flow path relative to ventilation zone (m)
            A_vent - Equivalent area of a vent (m2)
            delta_p_vent_ref -- reference pressure difference for vent (Pa)
            orientation -- The orientation of the vent (degrees)
            pitch -- The pitch of the vent (degrees)
            altitude -- altitude of dwelling above sea level (m)
        
        Method:
            - Based on Section 6.4.3.6 Airflow through vents from BS EN 16798-7
        """
        self.__h_path = h_path
        self.__A_vent = A_vent 
        self.__delta_p_vent_ref = delta_p_vent_ref
        self.__orientation = orientation
        self.__pitch = pitch
        self.__altitude = altitude
        self.__external_conditions = extcond
        self.__n_vent = 0.5 # Flow exponent for vents based on Section B.3.2.2 from BS EN 16798-7
        self.__C_D_vent = 0.6 # Discharge coefficient of vents based on B.3.2.1 from BS EN 16798-7
        self.__C_vent_path = self.calculate_flow_coeff_for_vent()
        self.__p_a_alt = adjust_air_density_for_altitude(altitude)

    def calculate_flow_coeff_for_vent(self):
        """The airflow coefficient of the vent calculated from equivalent area A_vent_i
        according to EN 13141-1 and EN 13141-2. 
        Based on Equation 59 from BS EN 16798-7.
        """
        #NOTE: The standard does not define what the below 3600 and 10000 are.
        C_vent_path = (3600/10000) \
                    * self.__C_D_vent * self.__A_vent \
                    * (2/p_a_ref)**0.5 * (1/self.__delta_p_vent_ref)**(self.__n_vent - 0.5)
        return C_vent_path

    def calculate_ventilation_through_vents_using_internal_p(
            self,
            u_site,
            T_e,
            T_z,
            C_vent_path,
            C_p_path,
            p_z_ref,
            ):
        """Calculate the airflow through vents from internal pressure

        Arguments:
        u_site -- wind velocity at zone level (m/s)
        T_e -- external air temperature (K)
        T_z -- thermal zone air temperature (K)
        C_vent_path -- wind pressure coefficient at height of the vent
        C_p_path -- wind pressure coefficient at the height of the window part
        p_z_ref -- internal reference pressure (Pa)
        """
        # Pressure_difference at the vent level
        delta_p_path = calculate_pressure_difference_at_an_airflow_path(
            self.__h_path,
            C_p_path,
            u_site,
            T_e,
            T_z,
            p_z_ref
        )

        # Air flow rate for each couple of height and wind pressure coeficient associated with vents.
        # Based on Equation 58
        qv_vent_path = C_vent_path * np.sign(delta_p_path) \
                     * abs(delta_p_path) \
                     ** self.__n_vent
        return qv_vent_path

    def calculate_flow_from_internal_p(self, u_site, T_z, p_z_ref, f_cross, shield_class):
        """Calculate the airflow through vents from internal pressure

        Arguments:
        u_site -- wind velocity at zone level (m/s)
        T_z -- thermal zone air temperature (K)
        p_z_ref -- internal reference pressure (Pa)
        f_cross -- boolean, dependant on if cross ventilation is possible or not
        shield_class -- indicates exposure to wind
        """
        wind_direction = self.__external_conditions.wind_direction()
        T_e = Celcius2Kelvin(self.__external_conditions.air_temp())

        # Wind pressure coefficient for the air flow path
        C_p_path = get_C_p_path(
            f_cross,
            shield_class,
            self.__h_path,
            wind_direction,
            self.__orientation,
            self.__pitch,
            )

        # Calculate airflow through each vent
        air_flow = self.calculate_ventilation_through_vents_using_internal_p(
            u_site,
            T_e,
            T_z,
            self.__C_vent_path ,
            C_p_path,
            p_z_ref,
            )

        # Sum airflows entering and leaving - based on Equation 60 and 61
        qv_in_through_vent = 0
        qv_out_through_vent = 0
        if air_flow >= 0:
            qv_in_through_vent += air_flow
        else:
            qv_out_through_vent += air_flow

        # Convert volume air flow rate to mass air flow rate
        qm_in_through_vent, qm_out_through_vent = convert_to_mass_air_flow_rate(
            qv_in_through_vent,
            qv_out_through_vent,
            T_e,
            T_z,
            self.__p_a_alt
        )
        return qm_in_through_vent, qm_out_through_vent


class Leaks:
    """ An object to represent Leaks """
    def __init__(
            self,
            extcond,
            h_path,
            delta_p_leak_ref,
            qv_delta_p_leak_ref,
            facade_direction,
            A_roof,
            A_facades,
            A_leak,
            altitude,
            ):
        """ Construct a Leaks object

        Arguments:
            extcond -- reference to ExternalConditions object
            h_path -- mid height of the air flow path relative to ventlation zone floor level
            delta_p_leak_ref -- Reference pressure difference (From pressure test e.g blower door = 50Pa)
            C_leak_path -- Coefficient of leakage in the external envelope
            qv_delta_p_leak_ref -- flow rate through
            facade_direction -- The direction of the facade the leak is on.
            A_roof -- Surface area of the roof of the ventilation zone (m2)
            A_facades -- Surface area of facades (m2)
            A_leak - Reference area of the envelope airtightness index qv_delta_p_leak_ref (depends on national context)
            altitude -- altitude of dwelling above sea level (m)
        
        Method:
            - - Based on Section 6.4.3.6 Airflow through leaks from BS EN 16798-7.
        
        """
        self.__h_path = h_path
        self.__delta_p_leak_ref = delta_p_leak_ref
        self.__A_roof = A_roof
        self.__A_facades = A_facades
        self.__A_leak = A_leak
        self.__qv_delta_p_leak_ref = qv_delta_p_leak_ref
        self.__facade_direction = facade_direction
        self.__altitude = altitude
        self.__external_conditions = extcond
        self.__n_leak = 0.667 # Flow exponent through leaks based on value in B.3.3.14
        self.__C_leak_path = self.calculate_flow_coeff_for_leak()
        self.__p_a_alt = adjust_air_density_for_altitude(altitude)

    def calculate_flow_coeff_for_leak(self):
        """The airflow coefficient of the leak """
        # C_leak - Leakage coefficient of ventilation zone
        C_leak = self.__qv_delta_p_leak_ref * self.__A_leak / (self.__delta_p_leak_ref) ** self.__n_leak

        # Leakage coefficient of roof, estimated to be proportional to ratio 
        # of surface area of the facades to that of the facades plus the roof.
        if self.__facade_direction != 'Windward' and self.__facade_direction != 'Leeward': #leak in roof
            C_leak_roof = C_leak * self.__A_roof / (self.__A_facades + self.__A_roof)
            C_leak_path = C_leak_roof #Table B.12
        # Leakage coefficient of facades, estimated to be proportional to ratio 
        # of surface area of the roof to that of the facades plus the roof.
        else: #leak in facades
            C_leak_facades = C_leak * self.__A_facades / (self.__A_facades + self.__A_roof)
            C_leak_path = 0.25 * C_leak_facades #Table B.12
        return C_leak_path

    def calculate_ventilation_through_leaks_using_internal_p(
            self,
            u_site,
            T_e,
            T_z,
            C_p_path,
            p_z_ref,
            ):
        """Calculate the airflow through leaks from internal pressure

        Arguments:
        u_site -- wind velocity at zone level (m/s)
        T_e -- external air temperature (K)
        T_z -- thermal zone air temperature (K)
        C_p_path -- wind pressure coefficient at the height of the window part
        p_z_ref -- internal reference pressure (Pa)
        """
        # For each couple of height and wind pressure coeficient associated with vents,
        # the air flow rate.
        delta_p_path = calculate_pressure_difference_at_an_airflow_path(
            self.__h_path,
            C_p_path,
            u_site,
            T_e,
            T_z,
            p_z_ref,
            )

        # Airflow through leaks based on Equation 62
        qv_leak_path = self.__C_leak_path * np.sign(delta_p_path) \
                     * abs(delta_p_path) \
                     ** self.__n_leak

        return qv_leak_path

    def calculate_flow_from_internal_p(self, u_site, T_z, p_z_ref, f_cross, shield_class):
        wind_direction = self.__external_conditions.wind_direction()
        T_e = Celcius2Kelvin(self.__external_conditions.air_temp())

        # Wind pressure coefficient for the air flow path
        C_p_path = get_C_p_path(f_cross,
                                shield_class,
                                self.__h_path,
                                wind_direction,
                                None,
                                None,
                                self.__facade_direction,
                                ) #TABLE from annex B
        # Calculate airflow through each leak
        qv_in_through_leak = 0
        qv_out_through_leak = 0
        air_flow = self.calculate_ventilation_through_leaks_using_internal_p(
            u_site,
            T_e,
            T_z,
            C_p_path,
            p_z_ref,
            )

        # Add airflow entering and leaving through leak
        if air_flow >= 0:
            qv_in_through_leak += air_flow
        else:
            qv_out_through_leak += air_flow

        # Convert volume air flow rate to mass air flow rate
        qm_in_through_leak, qm_out_through_leak = convert_to_mass_air_flow_rate(
            qv_in_through_leak,
            qv_out_through_leak,
            T_e,
            T_z,
            self.__p_a_alt
        )
        return qm_in_through_leak, qm_out_through_leak


class AirTerminalDevices:
    """An object to represent AirTerminalDevices"""
    def __init__(
            self,
            A_ATD,
            delta_p_ATD_ref,
            ):
        """
        Construct a AirTerminalDevices object

        Arguments:
            A_ATD -- equivalent area of the air terminal device (m2)
            delta_p_ATD_ref -- Reference pressure difference for an air terminal device (Pa)
        
        Method:
            Based on Section 6.4.3.2.2 from BS EN 16798-7
        """
        self.__C_D_ATD = 0.6 # Discharge coefficient for air terminal devices based on B.3.2.1
        self.__n_ATD = 0.5 # Flow exponent of air terminal devices based on B.3.2.2
        self.__A_ATD = A_ATD 
        self.__delta_p_ATD_ref = delta_p_ATD_ref
        self.__C_ATD_path = self.calculate_flow_coeff_for_ATD()

    def calculate_flow_coeff_for_ATD(self):
        """The airflow coefficient of the ATD is calculated 
        from the equivalent area A_vent value, according to
        EN 13141-1 and EN 13141-2.
        Equation 26 from BS EN 16798-7."""
        #NOTE: The standard does not define what the below 3600 and 10000 are.
        C_ATD_path = (3600/10000) * self.__C_D_ATD \
                                  * self.__A_ATD \
                                  * (2 / p_a_ref)**0.5 \
                                  * (1 / self.__delta_p_ATD_ref) \
                                  **(self.__n_ATD - 0.5)
        return C_ATD_path

    def calculate_pressure_difference_atd(self, qv_pdu):
        """The pressure loss at internal air terminal devices is calculated from
        the total air flow rate passing through the device.
        Equation 25 from BS EN 16798-7.
        Solving for qv_pdu.

        Arguments:
        qv_pdu - volume flow rate through passive and hybrid ducts.
        """
        delta_p_ATD = -(np.sign(qv_pdu)) \
                    * (abs(qv_pdu)/self.__C_ATD_path) \
                    ** (1/self.__n_ATD)
        return delta_p_ATD

class Cowls:
    """An object to represent Cowls"""
    def __init__(self,
                 height,
                 ):
        """
        Construct a Cowls object

        Arguments:
        height - height Between the top of the roof and the roof outlet in m (m)
        """
        self.__C_p_cowl_roof = 0 #Default B.3.3.5
        self.__delta_cowl_height = self.get_delta_cowl_height(height)

    def get_delta_cowl_height(self, height):
        """Interpreted Table B.9 from BS EN 16798-7
        Get values for delta_C_cowl_height.

        Arguments:
        height - height Between the top of the roof and the roof outlet in m (m)
        """
        if height < 0.5:
            delta_p_cowl_height = -0.0
        elif height >= 0.5 and height <= 1.0:
            delta_p_cowl_height = -0.1
        else:
            delta_p_cowl_height = -0.2
        return delta_p_cowl_height

class CombustionAppliances:
    """ An object to represent CombustionAppliances """
    def __init__(self,
                 supply_situation,
                 exhaust_situation,
                 fuel_type,
                 appliance_type,
                 ):
        """Construct a CombustionAppliances object
        
        Arguments:
        f_op_comp -- Operation requirement signal (combustion appliance) (0 =  OFF ; 1 = ON)
        supply_situation - Combustion air supply situation: 'room_air' or 'outside'
        exhaust_situation - flue gas exhaust situation: 'into_room', 'into_separate_duct' or 'into_mech_vent'
        f_ff -- combustion air flow factor
        P_h_fi - Combustion appliance heating fuel input power (kW)
        """
        self.__f_as = get_appliance_system_factor(supply_situation, exhaust_situation) #Combustion appliance system factor (0 or 1)
        self.__f_ff = get_fuel_flow_factor(fuel_type, appliance_type) #Fuel flow factor (0-5)

    def calculate_air_flow_req_for_comb_appliance(self, f_op_comp, P_h_fi):
        """Calculate additional air flow rate required for the operation of
        combustion appliance q_v_comb.

        Arguments:
        f_op_comp --  Operation requirement signal (combustion appliance) (0 =  OFF ; 1 = ON)
        P_h_fi -- Combustion appliance heating fuel input power (kW)
        """
        if f_op_comp == 1:
            q_v_comb = 3.6 * f_op_comp * self.__f_as * self.__f_ff * P_h_fi #(35)
            q_v_in_through_comb = 0 #(37)
            q_v_out_through_comb = -q_v_comb # (38)
            #temp associated with q_v_out_through_comb is ventilation zone temperature Tz
        else:
            q_v_in_through_comb = 0
            q_v_out_through_comb = 0
            #TODO flue is considered as vertical passive duct, standard formulas in CIBSE guide.
        return q_v_in_through_comb, q_v_out_through_comb
    


class MechanicalVentilation:
    """ An object to represent Mechanical Ventilation """
    def __init__(
            self,
            extcond,
            sup_air_flw_ctrl,
            sup_air_temp_ctrl,
            Q_H_des,
            Q_C_des,
            vent_type,
            specific_fan_power,
            design_outdoor_air_flow_rate,
            simulation_time,
            energy_supply_conn,
            total_volume,
            altitude,
            ctrl_intermittent_MEV = None,
            mvhr_eff = 0.0,
            theta_ctrl_sys = None, #Only required if sup_air_temp_ctrl = LOAD_COM
            ):
        """Construct a Mechanical Ventilation object
        
        Arguments:
            external_conditions -- reference to ExternalConditions object
            sup_air_flw_ctrl -- supply air flow rate control
            sup_air_temp_ctrl --supply air temperature control
            Q_H_des -- deisgn zone heating need to be covered by the mechanical ventilation system
            Q_C_des -- deisgn zone cooling need to be covered by the mechanical ventilation system
            vent_type -- ventilation system type
            specific_fan_power -- in W / (litre / second), inclusive of any in use factors
            design_outdoor_air_flow_rate -- design outdoor air flow rate in m3/h
            simulation_time -- reference to Simulation time object
            energy_supply_conn -- Energy supply connection
            total_volume  -- Total zone volume (m3)
            altitude -- altitude of dwelling above sea level (m)
            ctrl_intermittent_MEV -- reference to Control objject with boolean schedule
                                    defining when the MechVent should be on.
            mvhr_eff -- MVHR efficiency
            theta_ctrl_sys -- Temperature variation based on control system (K)
        
        """
        # Hard coded variables
        self.__f_ctrl = 1 #From table B.4, for residential buildings, default f_ctrl = 1
        self.__f_sys = 1.1 #From table B.5, f_sys = 1.1
        self.__E_v = 1 #Section B.3.3.7 defaults E_v = 1 (this is the assumption for perfect mixing)
        self.theta_z_t = 0 #TODO get Thermal zone temperature - used for LOAD
        self.__sup_air_flw_ctrl = "ODA" #TODO currently hard coded until load comp implemented
        self.__sup_air_temp_ctrl = "NO_CTRL"#TODO currently hard coded until load comp implemented

        # Arguments
        self.__external_conditions = extcond
        self.__Q_H_des = Q_H_des
        self.__Q_C_des = Q_C_des
        self.__theta_ctrl_sys = theta_ctrl_sys
        self.__vent_type = vent_type
        self.vent_type = vent_type
        self.total_volume = total_volume
        self.__ctrl_intermittent_MEV = ctrl_intermittent_MEV
        self.__sfp = specific_fan_power
        self.__simtime = simulation_time
        self.__energy_supply_conn = energy_supply_conn
        self.__altitude = altitude
        self.design_outdoor_air_flow_rate_m3_h = design_outdoor_air_flow_rate # in m3/h
        self.__mvhr_eff = mvhr_eff

        # Calculated variables
        self.__qv_ODA_req_design = self.calculate_required_outdoor_air_flow_rate()
        self.__p_a_alt = adjust_air_density_for_altitude(altitude)

    def calculate_required_outdoor_air_flow_rate(self):
        """Calculate required outdoor ventilation air flow rates.
        Equation 9 from BS EN 16798-7."""
        # Required outdoor air flow rate in m3/h
        qv_ODA_req = ((self.__f_ctrl * self.__f_sys) / self.__E_v) * self.design_outdoor_air_flow_rate_m3_h
        return qv_ODA_req

    def calc_req_ODA_flow_rates_at_ATDs(self):
        """Calculate required outdoor air flow rates at the air terminal devices
        Equations 10-17 from BS EN 16798-7
        Adjusted to be based on ventilation type instead of vent_sys_op.
        """
        if self.__vent_type== "MVHR":
            qv_SUP_req = self.__qv_ODA_req_design
            qv_ETA_req = -self.__qv_ODA_req_design
            # NOTE: Calculation of effective flow rate of external air (in func
            # calc_mech_vent_air_flw_rates_req_to_supply_vent_zone) assumes that
            # supply and extract are perfectly balanced (as defined above), so
            # any future change to this assumption will need to be considered
            # with that in mind
        elif self.__vent_type == "Intermittent MEV" or self.__vent_type == 'Centralised continuous MEV' or self.__vent_type == 'Decentralised continuous MEV':
            qv_SUP_req = 0
            qv_ETA_req = -self.__qv_ODA_req_design
        elif self.__vent_type == "PIV":
            qv_SUP_req = self.__qv_ODA_req_design
            qv_ETA_req = 0
        else:
            sys.exit('Unrecognised ventilation system type')

        return qv_SUP_req, qv_ETA_req

    def __f_op_v(self):
        """ Returns the fraction of the timestep for which the ventilation is running """
        if self.__vent_type == "Intermittent MEV":
            f_op_V = self.__ctrl_intermittent_MEV.setpnt()
            if not 0 <= f_op_V <= 1:
                sys.exit('Error f_op_V is not between 0 and 1.')
        elif self.__vent_type in ("Decentralised continuous MEV", "Centralised continuous MEV", "MVHR"):
            #Assumed to operate continuously
            f_op_V = 1
        else:
            sys.exit('Unknown mechanical ventlation system type')
        return f_op_V

    def calc_mech_vent_air_flw_rates_req_to_supply_vent_zone(
            self,
            T_z,
            ):
        """Calculate the air flow rates to and from the ventilation zone required from mechanical ventilation.
        T_z -- thermal zone temperature (K)
        """
        T_e = Celcius2Kelvin(self.__external_conditions.air_temp())
        # Required air flow at air terminal devices
        qv_SUP_req, qv_ETA_req = self.calc_req_ODA_flow_rates_at_ATDs()

        # Amount of air flow depends on controls
        if self.__sup_air_flw_ctrl == "ODA":
            f_op_V = self.__f_op_v()

            # Based on Equation 18 and 19
            qv_SUP_dis_req = f_op_V * qv_SUP_req
            qv_ETA_dis_req = f_op_V * qv_ETA_req
        elif self.__sup_air_flw_ctrl == "LOAD":
            #TODO load not implementented
            Q_H_V_req = 10 #TODO get (Q_H_V_req) heat to be supplied to thermal zone by the ventilation system at current timestep.
            Q_C_V_req = 10 ##TODO get (Q_C_V_req) heat to be extracted from thermal zone by the ventilation system at current timestep.
            p_a = adjust_air_density_for_altitude(self.__altitude)
            theta_SUP_dis_out = 10 #TODO get supply air temperature from previous timestep.
            Q_H_V_emission_losses = Q_H_V_req * self.theta_ctrl_sys / (self.theta_z_t - self.theta_e_comb)
            Q_C_V_emission_losses = Q_C_V_req * self.theta_ctrl_sys / (self.theta_z_t - self.theta_e_comb)
            qv_SUP_dis_req = self.f_op_V * max(
                (Q_H_V_req + Q_H_V_emission_losses) / (p_a * c_a * (theta_SUP_dis_out - T_z + T_0_abs)),
                (Q_C_V_req + Q_C_V_emission_losses) / (p_a * c_a * (T_z - T_0_abs - theta_SUP_dis_out)),
                qv_SUP_req
                )
            qv_SUP_dis_max_des = max (
                self.Q_H_des / (p_a * c_a * (theta_SUP_dis_out - T_z + T_0_abs)),
                self.Q_C_des / (p_a * c_a * (T_z - T_0_abs - theta_SUP_dis_out)),
                self.qv_ODA_req_design
                )
            qv_SUP_dis_req = min(qv_SUP_dis_req, qv_SUP_dis_max_des)

            if self.vent_type == "MVHR":
                qv_ETA_dis_req = -qv_SUP_dis_req
            else:
                qv_ETA_dis_req = 0
        else:
            sys.exit('Unknown sup_air_flw_ctrl type')

        # Calculate effective flow rate of external air
        # NOTE: Technically, the MVHR system supplies air at a higher
        # temperature than the outside air. However, it is simpler to
        # account for the heat recovery effect using an "equivalent" or
        # "effective" flow rate of external air
        qv_effective_heat_recovery_saving = qv_SUP_dis_req * self.__mvhr_eff

        # Convert volume air flow rate to mass air flow rate
        qm_SUP_dis_req, qm_ETA_dis_req = convert_to_mass_air_flow_rate(
            qv_SUP_dis_req,
            qv_ETA_dis_req,
            T_e,
            T_z,
            self.__p_a_alt
        )
        qm_in_effective_heat_recovery_saving \
            = convert_volume_flow_rate_to_mass_flow_rate(
                qv_effective_heat_recovery_saving,
                T_e,
                self.__p_a_alt,
                )

        return qm_SUP_dis_req, qm_ETA_dis_req, qm_in_effective_heat_recovery_saving

    def fans(self, zone_volume, total_volume, throughput_factor=1.0):
        """ Calculate gains and energy use due to fans
        zone_volume -- volume of the zone (m3)
        total_volume -- volume of the dwelling (m3)
        vent_type -- one of "Intermittent MEV", "Centralised continuous MEV",
            "Decentralised continuous MEV", "MVHR" or "PIV".
        """
        # Calculate energy use by fans
        fan_power_W = (self.__sfp * (self.__qv_ODA_req_design / seconds_per_hour) * litres_per_cubic_metre)\
                      * (zone_volume / total_volume)
        fan_energy_use_kWh = (fan_power_W / W_per_kW) * self.__simtime.timestep() * self.__f_op_v()
        if self.__vent_type in ("Intermittent MEV", "Centralised continuous MEV", "Decentralised continuous MEV"):
            #Fan energy use = 0
            supply_fan_energy_use_kWh = 0.0
            extract_fan_energy_use_in_kWh = fan_energy_use_kWh
        elif self.__vent_type == 'MVHR':
            #Balanced, therefore split power between extract and supply fans
            supply_fan_energy_use_kWh = fan_energy_use_kWh / 2
            extract_fan_energy_use_in_kWh = fan_energy_use_kWh / 2
        elif self.__vent_type == 'PIV':
            #Positive input, supply fans only
            supply_fan_energy_use_kWh = fan_energy_use_kWh
            extract_fan_energy_use_in_kWh = 0
        else:
            sys.exit("Invalid ventilation type input.")
        self.__energy_supply_conn.demand_energy(supply_fan_energy_use_kWh)
        self.__energy_supply_conn.demand_energy(extract_fan_energy_use_in_kWh)

        return supply_fan_energy_use_kWh / (W_per_kW * self.__simtime.timestep())

class InfiltrationVentilation:
    """ A class to represent Infiltration and Ventilation object """
    def __init__(
            self,
            extcond,
            simulation_time,
            f_cross,
            shield_class,
            terrain_class,
            average_roof_pitch,
            windows,
            vents,
            leaks,
            combustion_appliances,
            ATDs,
            mech_vents,
            detailed_output_heating_cooling,
            altitude,
            total_volume,
            ):
        """
        Constructs a InfiltrationVentilation object

        Arguments:
            external_conditions -- reference to ExternalConditions object
            simulation_time -- reference to SimulationTime object
            f_cross -- cross-ventilation factor
            shield_class -- indicates the exposure to wind of an air flow path on a facade 
                (can can be open, normal and shielded)
            ventilation_zone_height -- height of ventilation zone (m)
            windows -- list of windows
            vents -- list of vents
            leaks -- required inputs for leaks
            ATDs -- list of air terminal devices
            mech_vents -- list of mech vents
            altitude -- altitude of dwelling above sea level (m)
            total_volume -- total zone volume
        """
        self.__external_conditions = extcond
        self.__simulation_time = simulation_time
        self.__f_cross = f_cross
        self.__shield_class = shield_class
        self.__terrain_class = terrain_class
        self.__C_rgh_site = ter_class_to_roughness_coeff(terrain_class)
        self.ventilation_zone_height = leaks["ventilation_zone_height"]
        self.__windows = []
        self.__vents = []
        self.__leaks = []
        self.__combustion_appliances = []
        self.__ATDs = []
        self.__mech_vents = []
        for window in windows.values():
            self.__windows.append(window)
        for vent in vents.values():
            self.__vents.append(vent)
        self.__leaks = self.make_leak_objects(leaks, average_roof_pitch)
        for combustion_appliance in combustion_appliances.values():
            self.__combustion_appliances.append(combustion_appliance)
        for ATD in ATDs.values():
            self.__ATDs.append(ATD)
        for mech_vent in mech_vents.values():
            self.__mech_vents.append(mech_vent)
        self.__detailed_output_heating_cooling = detailed_output_heating_cooling
        self.__ventilation_detailed_results = []
        self.__p_a_alt = adjust_air_density_for_altitude(altitude)
        self.total_volume = total_volume

    def temp_supply(self):
        """ Calculate supply temperature of the air flow element """
        # NOTE: Technically, the MVHR system supplies air at a higher temperature
        # than the outside air, i.e.:
        #     temp_supply = self.__efficiency * temp_int_air \
        #                 + (1 - self.__efficiency) * self.__external_conditions.air_temp()
        # However, calculating this requires the internal air temperature, which
        # has not been calculated yet. Calculating this properly would require
        # the equation above to be added to the heat balance solver. Therefore,
        # it is simpler to adjust the heat transfer coefficient h_ve to account
        # for the heat recovery effect using an "equivalent" flow rate of
        # external air, which is done elsewhere
        return self.__external_conditions.air_temp()

    def calculate_total_volume_air_flow_rate_in(self, qm_in, external_air_density):
        """Calculate total volume air flow rate entering ventilation zone 
        Equation 68 from BS EN 16798-7"""
        return qm_in / external_air_density #from weather file?

    def calculate_total_volume_air_flow_rate_out(self, qm_out, zone_air_density):
        """Calculate total volume air flow rate leaving ventilation zone
        Equation 69 from BS EN 16798-7"""
        return qm_out / zone_air_density

    def make_leak_objects(self, leak, average_roof_pitch,):
        """Distribute leaks around the dwelling according to Table B.12 from BS EN 16798-7.
        Create 5 leak objects:
            At 0.25*Height of the Ventilation Zone in the Windward facade
            At 0.25*Height of the Ventilation Zone in the Leeward facade
            At 0.75*Height of the Ventilation Zone in the Windward facade
            At 0.75*Height of the Ventilation Zone in the Leeward facade
            At the Height of the Ventilation Zone in the roof
        Arguments:
        leak - dict of leaks input data from JSON file
        average_roof_pitch - calculated in project.py, average pitch of all roof elements weighted by area (degrees)
        """
        h_path1_2 = 0.25 * leak["ventilation_zone_height"]
        h_path3_4 = 0.75 * leak["ventilation_zone_height"]
        h_path5 = leak["ventilation_zone_height"]
        h_path_list = [h_path1_2, h_path1_2, h_path3_4, h_path3_4, h_path5]

        if average_roof_pitch < 10:
            roof_pitch = 'Roof10'
        elif average_roof_pitch <= 30:
            roof_pitch = 'Roof10_30'
        elif average_roof_pitch < 60:
            roof_pitch = 'Roof30'
        facade_direction = ['Windward', 'Leeward', 'Windward', 'Leeward', roof_pitch]

        leaklist = []
        for i in range(0,5):
            leak_obj = Leaks(
                self.__external_conditions,
                h_path_list[i],
                leak['test_pressure'],
                leak['test_result'],
                facade_direction[i],
                leak['area_roof'],
                leak['area_facades'],
                leak['env_area'],
                leak['altitude']
            )
            leaklist.append(leak_obj)
        return leaklist

    def calculate_qv_pdu(self, qv_pdu, p_z_ref, T_z, h_z):
        """Implicit solver for qv_pdu"""
        return fsolve(self.implicit_formula_for_qv_pdu, qv_pdu, args = (p_z_ref, T_z, h_z)) #returns qv_pdu

    def implicit_formula_for_qv_pdu(self, qv_pdu, p_z_ref, T_z, h_z):
        """Implicit formula solving for qv_pdu as unknown.
        Equation 30 from BS EN 16798-7
        Arguments:
        qv_pdu -- volume flow rate from passive and hybrid ducts (m3/h)
        p_z_ref -- internal reference pressure (Pa)
        T_z -- thermal zone temperature (K)
        h_z -- height of ventilation zone (m)
        """
        T_e = Celcius2Kelvin(self.__external_conditions.air_temp())
        external_air_density = air_density_at_temp(T_e, self.__p_a_alt)
        zone_air_density = air_density_at_temp(T_z, self.__p_a_alt)

        #TODO Standard isn't clear if delta_p_ATD can be totalled or not. 
        delta_p_ATD_list = [atd.calculate_pressure_difference_atd(qv_pdu) for atd in self.__ATDs]
        delta_p_ATD = sum(delta_p_ATD_list)

        # Stack effect in passive and hybrid duct. As there is no air transfer 
        # between levels of the ventilation zone Equation B.1 is used.
        h_pdu_stack = h_z + 2

        #TODO include delta_p_dpu and delta_p_cowl in the return.
        return delta_p_ATD - p_z_ref - h_pdu_stack * g * (external_air_density - zone_air_density)
    
    def calculate_internal_reference_pressure(
                                            self,
                                            intial_p_z_ref_guess,
                                            temp_int_air,
                                            R_w_arg = None,
                                            ):
        """The root scalar function will iterate until it finds a value of p_z_ref 
        that satisfies the mass balance equation.
        The root scalar solver allows a range of intervals to be entered. 
        The loop begins with a small interval to start with and if no solution is 
        found or the boundary is too small for to cause a sign change then a wider 
        interval is used until a solution is found.
        """
        interval_expansion_list = [1, 5, 10, 15, 20, 40, 50, 100, 200]
        for interval_expansion in interval_expansion_list:
            try:
                with warnings.catch_warnings():
                    warnings.filterwarnings('error',category=RuntimeWarning)
                    sol = root_scalar(
                        self.implicit_mass_balance_for_internal_reference_pressure,
                        bracket=[intial_p_z_ref_guess - interval_expansion, intial_p_z_ref_guess + interval_expansion],
                        args=(temp_int_air, R_w_arg),
                        method='brentq'
                        )
                p_z_ref = sol.root
                return p_z_ref
            except ValueError as e:
                if str(e) == "f(a) and f(b) must have different signs":
                    continue
            except RuntimeWarning as e:
                continue
            except Exception as e:
                sys.exit("Solver failed: " +str(e))

    def implicit_mass_balance_for_internal_reference_pressure(
            self,
            p_z_ref,
            temp_int_air,
            R_w_arg_min_max,
            flag = None,
            ):
        """Used in calculate_internal_reference_pressure function for p_z_ref solve"""
        qm_in, qm_out, _ = self.__implicit_mass_balance_for_internal_reference_pressure_components(
            p_z_ref,
            temp_int_air,
            R_w_arg_min_max,
            flag,
            )
        return qm_in + qm_out

    def incoming_air_flow(
            self,
            p_z_ref,
            temp_int_air,
            R_w_arg_min_max,
            reporting_flag = None,
            report_effective_flow_rate = False,
            ):
        """ Calculate incoming air flow, in m3/hr, at specified conditions """
        qm_in, _, qm_effective_flow_rate \
            = self.__implicit_mass_balance_for_internal_reference_pressure_components(
                p_z_ref,
                temp_int_air,
                R_w_arg_min_max,
                reporting_flag,
                )
        if report_effective_flow_rate:
            qm_in -= qm_effective_flow_rate
        return convert_mass_flow_rate_to_volume_flow_rate(
            qm_in,
            Celcius2Kelvin(self.__external_conditions.air_temp()),
            self.__p_a_alt,
            )

    def __implicit_mass_balance_for_internal_reference_pressure_components(
            self,
            p_z_ref,
            temp_int_air,
            R_w_arg_min_max,
            reporting_flag = None,
            ):
        """Implicit mass balance for calculation of the internal reference pressure
        Equation 67 from BS EN 16798-7.

        Arguments:
        p_z_ref -- internal reference pressure (Pa)
        temp_int_air -- temperature of the intake air (K)
        reporting_flag -- flag used to give more detailed ventilation outputs (None = no additional reporting)

        Key Variables:
        qm_SUP_to_vent_zone - Supply air mass flow rate going to ventilation zone
        qm_ETA_from_vent_zone - Extract air mass flow rate from a ventilation zone
        qm_in_through_comb - Air mass flow rate entering through combustion appliances
        qm_out_through_comb - Air mass flow rate leaving through combustion appliances
        qm_in_through_passive_hybrid_ducts - Air mass flow rate entering through passive or hybrid duct
        qm_out_through_passive_hybrid_ducts - Air mass flow rate leaving through passive or hybrid duct
        qm_in_through_window_opening - Air mass flow rate entering through window opening
        qm_out_through_window_opening - Air mass flow rate leaving through window opening
        qm_in_through_vents - Air mass flow rate entering through vents (openings in the external envelope)
        qm_out_through_vents - Air mass flow rate leaving through vents (openings in the external envelope)
        qm_in_through_leaks - Air mass flow rate entering through envelope leakage
        qm_out_through_leaks - Air mass flow rate leaving through envelope leakage
        """
        wind_speed = self.__external_conditions.wind_speed()
        u_site = wind_speed_at_zone_level(self.__C_rgh_site, wind_speed)
        T_e = Celcius2Kelvin(self.__external_conditions.air_temp())
        T_z = Celcius2Kelvin(temp_int_air)
        qm_in_through_window_opening = 0
        qm_out_through_window_opening = 0
        qm_in_through_vents = 0
        qm_out_through_vents = 0
        qm_in_through_leaks = 0
        qm_out_through_leaks = 0
        qm_in_through_comb = 0
        qm_out_through_comb = 0
        qm_in_through_passive_hybrid_ducts = 0
        qm_out_through_passive_hybrid_ducts = 0
        qm_SUP_to_vent_zone = 0
        qm_ETA_from_vent_zone = 0
        qm_in_effective_heat_recovery_saving_total = 0.0

        for window in self.__windows:
            qm_in, qm_out = window.calculate_flow_from_internal_p(u_site, T_z, p_z_ref, self.__f_cross, self.__shield_class, R_w_arg_min_max)
            qm_in_through_window_opening += qm_in
            qm_out_through_window_opening += qm_out
        for vent in self.__vents:
            qm_in, qm_out = vent.calculate_flow_from_internal_p(u_site, T_z, p_z_ref, self.__f_cross, self.__shield_class)
            qm_in_through_vents += qm_in
            qm_out_through_vents += qm_out
        for leak in self.__leaks:
            qm_in, qm_out = leak.calculate_flow_from_internal_p(u_site, T_z, p_z_ref, self.__f_cross, self.__shield_class)
            qm_in_through_leaks += qm_in
            qm_out_through_leaks += qm_out
        for ATD in self.__ATDs:
            qv_pdu_initial = 0 #TODO get from prev timestep
            h_z = self.ventilation_zone_height
            qv_pdu = self.calculate_qv_pdu(qv_pdu_initial, p_z_ref, T_z, h_z)
            if qv_pdu >= 0:
                qv_pdu_in = qv_pdu
                qv_pdu_out = 0
            else:
                qv_pdu_in = 0
                qv_pdu_out = qv_pdu
            qm_in_through_phds, qm_out_through_phds = convert_to_mass_air_flow_rate(
                qv_pdu_in,
                qv_pdu_out,
                T_e,
                T_z,
                self.__p_a_alt
            )
            qm_in_through_passive_hybrid_ducts += qm_in_through_phds
            qm_out_through_passive_hybrid_ducts += qm_out_through_phds
        for combustion_appliance in self.__combustion_appliances:
            P_h_fi = 0 #TODO to work out from previous zone temperature? - Combustion appliance heating fuel input power
            f_op_comb = 1 #TODO work out what turns the appliance on or off. Schedule or Logic?
            qv_in, qv_out = combustion_appliance.calculate_air_flow_req_for_comb_appliance(f_op_comb, P_h_fi)
            qm_in_comb, qm_out_comb = convert_to_mass_air_flow_rate(
                qv_in,
                qv_out,
                T_e,
                T_z,
                self.__p_a_alt
            )
            qm_in_through_comb += qm_in_comb
            qm_out_through_comb += qm_out_comb
        for mech_vent in self.__mech_vents:
            qm_SUP, qm_ETA, qm_in_effective_heat_recovery_saving \
                = mech_vent.calc_mech_vent_air_flw_rates_req_to_supply_vent_zone(
                    T_z,
                    )
            qm_SUP_to_vent_zone += qm_SUP
            qm_ETA_from_vent_zone += qm_ETA
            qm_in_effective_heat_recovery_saving_total += qm_in_effective_heat_recovery_saving

        # Calculate the total mass airflow rate in and out of the zones
        qm_in \
            = qm_in_through_window_opening \
            + qm_in_through_vents \
            + qm_in_through_leaks \
            + qm_in_through_comb \
            + qm_in_through_passive_hybrid_ducts \
            + qm_SUP_to_vent_zone
        qm_out \
            = qm_out_through_window_opening \
            + qm_out_through_vents \
            + qm_out_through_leaks \
            + qm_out_through_comb \
            + qm_out_through_passive_hybrid_ducts \
            + qm_ETA_from_vent_zone
        
        # Output detailed ventilation file
        if self.__detailed_output_heating_cooling:
            incoming_air_flow = convert_mass_flow_rate_to_volume_flow_rate(
                qm_in,
                Celcius2Kelvin(self.__external_conditions.air_temp()),
                self.__p_a_alt,
                )
            air_changes_per_hour = incoming_air_flow / self.total_volume

            if reporting_flag is not None:
                self.__ventilation_detailed_results.append(
                    [self.__simulation_time.index(), reporting_flag,
                     incoming_air_flow, self.total_volume, air_changes_per_hour,
                     temp_int_air,p_z_ref,
                     qm_in_through_window_opening, qm_out_through_window_opening,
                     qm_in_through_vents, qm_out_through_vents,
                     qm_in_through_leaks, qm_out_through_leaks,
                     qm_in_through_comb, qm_out_through_comb,
                     qm_in_through_passive_hybrid_ducts, qm_out_through_passive_hybrid_ducts,
                     qm_SUP_to_vent_zone, qm_ETA_from_vent_zone, qm_in_effective_heat_recovery_saving_total,
                     qm_in , qm_out]
                    )

        return qm_in, qm_out, qm_in_effective_heat_recovery_saving_total

    def output_vent_results(self):
        """ Return the data dictionary containing detailed ventilation results"""
        return self.__ventilation_detailed_results

