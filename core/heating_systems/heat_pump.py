#!/usr/bin/env python3

"""
This module provides objects to represent heat pumps and heat pump test data.
The calculations are based on the DAHPSE method developed for generating PCDB
entries for SAP 2012 and SAP 10. DAHPSE was based on a draft of
BS EN 15316-4-2:2017 and is described in the SAP calculation method CALCM-01.
"""

# Standard library imports
import sys
from copy import deepcopy
from enum import Enum, auto

# Third-party imports
import numpy as np
from numpy.polynomial.polynomial import polyfit

# Local imports
from core.units import Celcius2Kelvin, Kelvin2Celcius, hours_per_day, W_per_kW, kJ_per_kWh, seconds_per_minute
from core.heating_systems.boiler import BoilerServiceSpace, BoilerServiceWaterRegular
from core.heating_systems.emitters import convert_flow_to_return_temp
from core.schedule import expand_schedule
from core.material_properties import WATER

# Constants
N_EXER = 3.0


# Data types

class SourceType(Enum):
    GROUND = auto()
    OUTSIDE_AIR = auto()
    EXHAUST_AIR_MEV = auto()
    EXHAUST_AIR_MVHR = auto()
    EXHAUST_AIR_MIXED = auto()
    WATER_GROUND = auto()
    WATER_SURFACE = auto()
    HEAT_NETWORK = auto()

    @classmethod
    def from_string(cls, strval):
        if strval == 'Ground':
            return cls.GROUND
        elif strval == 'OutsideAir':
            return cls.OUTSIDE_AIR
        elif strval == 'ExhaustAirMEV':
            return cls.EXHAUST_AIR_MEV
        elif strval == 'ExhaustAirMVHR':
            return cls.EXHAUST_AIR_MVHR
        elif strval == 'ExhaustAirMixed':
            return cls.EXHAUST_AIR_MIXED
        elif strval == 'WaterGround':
            return cls.WATER_GROUND
        elif strval == 'WaterSurface':
            return cls.WATER_SURFACE
        elif strval == 'HeatNetwork':
            return cls.HEAT_NETWORK
        else:
            sys.exit('SourceType (' + str(strval) + ') not valid.')
            # TODO Exit just the current case instead of whole program entirely?

    @classmethod
    def is_exhaust_air(cls, source_type):
        # If string has been provided, convert to SourceType before running check
        if isinstance(source_type, str):
            source_type = cls.from_string(source_type)

        if source_type in (cls.EXHAUST_AIR_MEV, cls.EXHAUST_AIR_MVHR, cls.EXHAUST_AIR_MIXED):
            return True
        elif source_type in (
            cls.GROUND,
            cls.OUTSIDE_AIR,
            cls.WATER_GROUND,
            cls.WATER_SURFACE,
            cls.HEAT_NETWORK,
            ):
            return False
        else:
            sys.exit('SourceType (' + str(source_type) + ') not defined as exhaust air or not.')

    @classmethod
    def source_fluid_is_air(cls, source_type):
        # If string has been provided, convert to SourceType before running check
        if isinstance(source_type, str):
            source_type = cls.from_string(source_type)

        if source_type \
        in (cls.OUTSIDE_AIR, cls.EXHAUST_AIR_MEV, cls.EXHAUST_AIR_MVHR, cls.EXHAUST_AIR_MIXED):
            return True
        elif source_type in (cls.GROUND, cls.WATER_GROUND, cls.WATER_SURFACE, cls.HEAT_NETWORK):
            return False
        else:
            sys.exit( 'SourceType (' + str(source_type) \
                    + ') not defined as having air as source fluid or not.')

    @classmethod
    def source_fluid_is_water(cls, source_type):
        # If string has been provided, convert to SourceType before running check
        if isinstance(source_type, str):
            source_type = cls.from_string(source_type)

        if source_type in (cls.GROUND, cls.WATER_GROUND, cls.WATER_SURFACE, cls.HEAT_NETWORK):
            return True
        elif source_type \
        in (cls.OUTSIDE_AIR, cls.EXHAUST_AIR_MEV, cls.EXHAUST_AIR_MVHR, cls.EXHAUST_AIR_MIXED):
            return False
        else:
            sys.exit( 'SourceType (' + str(source_type) \
                    + ') not defined as having water as source fluid or not.')


class SinkType(Enum):
    AIR = auto()
    WATER = auto()

    @classmethod
    def from_string(cls, strval):
        if strval == 'Air':
            return cls.AIR
        elif strval == 'Water':
            return cls.WATER
        else:
            sys.exit('SinkType (' + str(strval) + ') not valid.')
            # TODO Exit just the current case instead of whole program entirely?

class BackupCtrlType(Enum):
    NONE = auto()
    TOPUP = auto()
    SUBSTITUTE = auto()

    @classmethod
    def from_string(cls, strval):
        if strval == 'None':
            return cls.NONE
        elif strval == 'TopUp':
            return cls.TOPUP
        elif strval == 'Substitute':
            return cls.SUBSTITUTE
        else:
            sys.exit('BackupType (' + str(strval) + ') not valid.')
            # TODO Exit just the current case instead of whole program entirely?

class ServiceType(Enum):
    WATER = auto()
    SPACE = auto()


# Free functions

def carnot_cop(temp_source, temp_outlet, temp_diff_limit_low=None):
    """ Calculate Carnot CoP based on source and outlet temperatures (in Kelvin) """
    temp_diff = temp_outlet - temp_source
    if temp_diff_limit_low is not None:
        temp_diff = max (temp_diff, temp_diff_limit_low)
    return temp_outlet / temp_diff

def interpolate_exhaust_air_heat_pump_test_data(
        throughput_exhaust_air,
        hp_dict_test_data,
        source_type,
        ):
    """ Interpolate between test data records for different air flow rates
    
    Arguments:
    throughput_exhaust_air -- throughput (m3/h) of exhaust air
    hp_dict_test_data
        -- list of dictionaries of heat pump test data, each with the following elements:
                - air_flow_rate
                - test_letter
                - capacity
                - cop
                - degradation_coeff
                - design_flow_temp (in Celsius)
                - temp_outlet (in Celsius)
                - temp_source (in Celsius)
                - temp_test (in Celsius)
    """
    # Split test records into different lists by air flow rate
    test_data_by_air_flow_rate = {}
    for test_data_record in hp_dict_test_data:
        if test_data_record['air_flow_rate'] not in test_data_by_air_flow_rate.keys():
            # Initialise list for this air flow rate if it does not already exist
            test_data_by_air_flow_rate[test_data_record['air_flow_rate']] = []
        test_data_by_air_flow_rate[test_data_record['air_flow_rate']].append(test_data_record)

    # Check that all lists have same combinations of design flow temp and test letter
    fixed_temps_and_test_letters = None
    for air_flow_rate, test_data_record_list in test_data_by_air_flow_rate.items():
        # Find and save all the combinations of design flow temp and test letter
        # for this air flow rate
        fixed_temps_and_test_letters_this = []
        for test_data_record in test_data_record_list:
            fixed_temps_and_test_letters_this.append(
                ( test_data_record['design_flow_temp'],
                  test_data_record['test_letter'],
                  test_data_record['temp_outlet'],
                  test_data_record['temp_source'],
                  test_data_record['temp_test'],
                ))

        if fixed_temps_and_test_letters is None:
            # If we are on the first iteration of the loop, save the list of
            # design flow temps and test letters from this loop for comparison
            # in subsequent loops
            fixed_temps_and_test_letters = fixed_temps_and_test_letters_this
        else:
            # If we are not on the first iteration of the loop, check that same
            # design flow temps and test letters are present for this air flow
            # rate and the first one
            assert set(fixed_temps_and_test_letters) \
                == set(fixed_temps_and_test_letters_this)

    # Construct test data records interpolated by air flow rate
    air_flow_rates_ordered = sorted(test_data_by_air_flow_rate.keys())
    hp_dict_test_data_interp_by_air_flow_rate = []
    for design_flow_temp, test_letter, temp_outlet, temp_source, temp_test \
    in fixed_temps_and_test_letters:
        # Create lists of test data values ordered by air flow rate
        capacity_list = []
        cop_list = []
        degradation_coeff_list = []
        ext_air_ratio_list = []
        for air_flow_rate in air_flow_rates_ordered:
            for test_record in test_data_by_air_flow_rate[air_flow_rate]:
                if test_record['design_flow_temp'] == design_flow_temp \
                and test_record['test_letter'] == test_letter:
                    capacity_list.append(test_record['capacity'])
                    cop_list.append(test_record['cop'])
                    degradation_coeff_list.append(test_record['degradation_coeff'])
                    if source_type == SourceType.EXHAUST_AIR_MIXED:
                        ext_air_ratio_list.append(test_record['eahp_mixed_ext_air_ratio'])

        # Interpolate test data by air flow rate
        capacity = np.interp(throughput_exhaust_air, air_flow_rates_ordered, capacity_list)
        cop      = np.interp(throughput_exhaust_air, air_flow_rates_ordered, cop_list)
        degradation_coeff \
                 = np.interp(throughput_exhaust_air, air_flow_rates_ordered, degradation_coeff_list)
        if source_type == SourceType.EXHAUST_AIR_MIXED:
            ext_air_ratio \
                     = np.interp(throughput_exhaust_air, air_flow_rates_ordered, ext_air_ratio_list)
        else:
            ext_air_ratio = None

        # Construct interpolated test data record 
        hp_dict_test_data_interp_by_air_flow_rate.append({
            "test_letter": test_letter,
            "capacity": capacity,
            "cop": cop,
            "degradation_coeff": degradation_coeff,
            "design_flow_temp": design_flow_temp,
            "temp_outlet": temp_outlet,
            "temp_source": temp_source,
            "temp_test": temp_test,
            "ext_air_ratio": ext_air_ratio,
            })

    # Find lowest air flow rate in test data
    lowest_air_flow_rate_in_test_data = min(test_data_by_air_flow_rate.keys())

    return lowest_air_flow_rate_in_test_data, hp_dict_test_data_interp_by_air_flow_rate


# Classes

class BufferTank:
    """ A class for Buffer Tanks objects linked to heat pumps.
    ...
    """
    #part of the thermal losses transmitted to the room
    #Taken from Storage tank module: Method A from BS EN 15316-5:2017
    __f_sto_m = 0.75

    def __init__(
            self, 
            daily_losses,
            volume,
            pump_fixed_flow_rate, 
            pump_power_at_flow_rate,
            number_of_zones,
            simulation_time,
            contents=WATER,
            output_detailed_results=False,
            ):
        """ Construct a BufferTank object

        Arguments:
        daily_losses -- standing loss from the buffer tank under standard test condition in kWh/day
        volume       -- volume of the buffer tank in litres
        pump_fixed_flow_rate    -- flow rate of the buffer tank - emitters loop in l/min
        pump_power_at_flow_rate       -- pump power of the buffer tank - emitters loop in W
        number_of_zones    -- number of zones in the project
        simulation_time -- reference to SimulationTime object
        contents        -- reference to MaterialProperties object
        output_detailed_results -- if true, save detailed results from each timestep
                                   for later reporting
                                   TODO: Detailed reporting has not yet been implemented for Buffer Tanks
        """
        self.__daily_losses = daily_losses
        self.__volume = volume
        self.__pump_fixed_flow_rate = pump_fixed_flow_rate
        self.__pump_power_at_flow_rate = pump_power_at_flow_rate
        self.__number_of_zones = number_of_zones
        self.__simulation_time = simulation_time
        self.__contents=contents
        self.__service_results = []
        
        self.__buffer_emitter_circ_flow_rate = 0.0    
        self.__hp_buffer_circ_flow_rate = 0.0
        self.__print_warning = True # To avoid multiple warning messages for undersized pump
        #water specific heat in kWh/kg.K
        self.__Cp = contents.specific_heat_capacity_kWh()
        #volumic mass in kg/litre
        self.__rho = contents.density()
        
        self.__Q_heat_loss_buffer_rbl = 0.0
        self.__track_buffer_loss = 0.0
        
        # If detailed results are to be output, initialise list
        if output_detailed_results:
            self.__detailed_results = []
        else:
            self.__detailed_results = None

    def update_buffer_loss(self,buffer_loss):
        self.__track_buffer_loss = deepcopy(buffer_loss)

    def get_buffer_loss(self):
        return(self.__track_buffer_loss)
            
    def thermal_losses(self,temp_ave_buffer,temp_rm_prev):
        """Thermal losses are calculated with respect to the impact of the temperature set point"""
        temp_set_ref = 65
        temp_amb_ref = 20
        H_buffer_ls = (1000 * self.__daily_losses) / (24 * (temp_set_ref - temp_amb_ref)) #Specific loss in W/K
        heat_loss_buffer_W = (temp_ave_buffer - temp_rm_prev) * H_buffer_ls #Absolute loss in W
        heat_loss_buffer_kWh = heat_loss_buffer_W / 1000 * self.__simulation_time.timestep() / self.__number_of_zones 
        # TODO: update when zoning sorted out. Hopefully then there will be no need to divide by the number of zones. 

        #recoverable heat losses from buffer - kWh
        self.__Q_heat_loss_buffer_rbl = deepcopy(heat_loss_buffer_kWh * self.__f_sto_m)
        return heat_loss_buffer_kWh

    def internal_gains(self):
        """ Return the recoverable heat losses as internal gain for the current timestep in W"""
        return self.__Q_heat_loss_buffer_rbl * W_per_kW / self.__simulation_time.timestep() * self.__number_of_zones 
        # multiply by number of zones because this is only called once (not per zone like thermal_loss) for internal gains 

    def calc_buffer_tank(self,service_name,emitters_data_for_buffer_tank):
        temp_rm_prev = emitters_data_for_buffer_tank['temp_rm_prev']

        if emitters_data_for_buffer_tank['power_req_from_buffer_tank'] > 0.0:
            temp_emitter_req = emitters_data_for_buffer_tank['temp_emitter_req']
            self.__temp_ave_buffer = deepcopy(temp_emitter_req)
            
            # call to calculate thermal losses
            heat_loss_buffer_kWh = self.thermal_losses(temp_emitter_req,temp_rm_prev)
            
            # We are calculating the delta T needed to give the heat output required from the emitters
            # E = m C dT rearranged to make delta T the subject
            deltaT_buffer = (emitters_data_for_buffer_tank['power_req_from_buffer_tank'] \
                             / (self.__pump_fixed_flow_rate / seconds_per_minute)) / (self.__Cp * self.__rho * kJ_per_kWh)
                             
            deltaT_hp_to_buffer = emitters_data_for_buffer_tank['target_flow_temp'] - \
                                  convert_flow_to_return_temp(emitters_data_for_buffer_tank['target_flow_temp'])
                                  
            # Consistent with the 6/7ths rule of thumb used in 'emitters.py' to get the return temperature.
            # In future, deltaT could be made an input - since it is set by system commissioner. 
            # 1/7 gives values in the range 5-10K, which is consistent with MCS best practice guidance:
            # https://mcscertified.com/wp-content/uploads/2020/07/Heat-Pump-Guide.pdf#:~:text=URL%3A%20https%3A%2F%2Fmcscertified.com%2Fwp
            # 'target_flow_temp' is from the emitters calculation, but is a figure per zone. 
            # But this should be looked at again if the there is an update to zoning made later.
            buffer_flow_temp = temp_emitter_req + 0.5 * deltaT_buffer
            buffer_return_temp = temp_emitter_req - 0.5 * deltaT_buffer
            theoretical_hp_return_temp = buffer_return_temp
            theoretical_hp_flow_temp = theoretical_hp_return_temp + deltaT_hp_to_buffer
        
            flow_temp_increase_due_to_buffer = max(0,theoretical_hp_flow_temp - buffer_flow_temp)
            # TODO Consider circumstances where the flow rate on the hp-buffer loop exceed the flow rate in the 
            #      buffer-emitter loop. We think in this circumstances the flow temperature would match rather than
            #      the return temperature, changing the logic of the methodology.
            #      This would happen in situation of very cold weather which exceed the design conditions for the 
            #      buffer-emitters loop sizing/flow rate

            # If detailed results are to be output, save the results from the current timestep
            self.__service_results.append({
            'service_name': service_name + "_buffer_tank",
            'power_req_from_buffer_tank': emitters_data_for_buffer_tank['power_req_from_buffer_tank'],
            'temp_emitter_req': emitters_data_for_buffer_tank['temp_emitter_req'],
            'buffer_emitter_circ_flow_rate': self.__pump_fixed_flow_rate,
            'flow_temp_increase_due_to_buffer': flow_temp_increase_due_to_buffer,
            'pump_power_at_flow_rate': self.__pump_power_at_flow_rate,
            'heat_loss_buffer_kWh': heat_loss_buffer_kWh,
            })

            if self.__detailed_results is not None:
                self.__detailed_results.append(self.__service_results)
        else:
            # call to calculate cool down losses
            heat_loss_buffer_kWh = self.thermal_losses(self.__temp_ave_buffer,temp_rm_prev)
            heat_capacity_buffer = self.__volume * (self.__rho * self.__Cp * kJ_per_kWh)
            if self.__temp_ave_buffer is not None:
                temp_loss = heat_loss_buffer_kWh / (heat_capacity_buffer / kJ_per_kWh)
                new_temp_ave_buffer = self.__temp_ave_buffer - temp_loss
                
            else:
                new_temp_ave_buffer = convert_flow_to_return_temp(emitters_data_for_buffer_tank['design_flow_temp'])
                
            self.__temp_ave_buffer = deepcopy(new_temp_ave_buffer)
            
            self.__service_results.append({
                'service_name': service_name + "_buffer_tank",
                'power_req_from_buffer_tank': 0.0,
                'temp_emitter_req': emitters_data_for_buffer_tank['temp_emitter_req'],
                'buffer_emitter_circ_flow_rate': self.__pump_fixed_flow_rate,
                'flow_temp_increase_due_to_buffer': 0.0,
                'pump_power_at_flow_rate': 0.0,
                'heat_loss_buffer_kWh': heat_loss_buffer_kWh,
                })
            if self.__detailed_results is not None:
                self.__detailed_results.append(self.__service_results)
                
        return self.__service_results
        
    
class HeatPumpTestData:
    """ An object to represent EN 14825 test data for a heat pump.

    This object stores the data and provides functions to look up values from
    the correct data records for the conditions being modelled.
    """

    __test_letters_non_bivalent = ['A', 'B', 'C', 'D']
    __test_letters_all = ['A','B','C','D','F']

    def __init__(self, hp_testdata_dict_list):
        """ Construct a HeatPumpTestData object

        Arguments:
        hp_testdata_dict_list
            -- list of dictionaries of heat pump test data, each with the following elements:
                - test_letter
                - capacity
                - cop
                - degradation_coeff
                - design_flow_temp (in Celsius)
                - temp_outlet (in Celsius)
                - temp_source (in Celsius)
                - temp_test (in Celsius)
        """
        def duplicates(a, b):
            """ Determine whether records a and b are duplicates """
            return (a['temp_test'] == b['temp_test'] \
                and a['design_flow_temp'] == b['design_flow_temp'])

        # Keys will be design flow temps, values will be lists of dicts containing the test data
        self.__testdata = {}

        # A separate list of design flow temps is required because it can be
        # sorted, whereas the dict above can't be (at least before Python 3.7)
        self.__dsgn_flow_temps = []
        # Dict to count duplicate records for each design flow temp
        dupl = {}

        # Read the test data records
        # Work on a deep copy of the input data structure in case the original
        # is used to init other objects (or the same object multiple times
        # e.g. during testing)
        for hp_testdata_dict in deepcopy(hp_testdata_dict_list):
            dsgn_flow_temp = hp_testdata_dict['design_flow_temp']

            # When a new design flow temp is encountered, add it to the lists/dicts
            if dsgn_flow_temp not in self.__dsgn_flow_temps:
                self.__dsgn_flow_temps.append(dsgn_flow_temp)
                self.__testdata[dsgn_flow_temp] = []
                dupl[dsgn_flow_temp] = 0

            # Check for duplicate records
            duplicate = False
            for d in self.__testdata[dsgn_flow_temp]:
                if duplicates(hp_testdata_dict, d):
                    duplicate = True
                    # Increment count of number of duplicates for this design flow temp
                    # Handle records with same inlet temp
                    # Cannot process a row at the same inlet temperature (div
                    # by zero error during interpolation), so we add a tiny
                    # amount to the temperature (to 10DP) for each duplicate
                    # found.
                    # TODO Why do we need to alter the duplicate record? Can we
                    #      not just eliminate it?
                    hp_testdata_dict['temp_test']   += 0.0000000001
                    hp_testdata_dict['temp_source'] += 0.0000000001
                    # TODO The adjustment to temp_source is in the python
                    #      implementation of DAHPSE but not in the spreadsheet
                    #      implementation. Given that temp_source can be the
                    #      same for all test records anyway, is this adjustment
                    #      needed?
            # This increment has to be after loop to avoid multiple-counting
            # when there are 3 or more duplicates. E.g. if there are already 2
            # records that are the same, then when adding a third that is the
            # same, we only want to increment the counter by 1 (for the record
            # we are adding) and not 2 (the number of existing records the new
            # record duplicates).
            if duplicate:
                dupl[dsgn_flow_temp] += 1

            # Add the test record to the data structure, under the appropriate design flow temp
            self.__testdata[dsgn_flow_temp].append(hp_testdata_dict)

        # Check the number of test records is as expected
        # - Up to 4 design flow temps (EN 14825:2018 defines these as 35, 45, 55, 65)
        # - 4 or 5 distinct records for each flow temp
        # TODO Is there any reason the model couldn't handle more than 5 test
        #      records if data is available? Could/should we relax the
        #      restriction below?
        if len(self.__dsgn_flow_temps) < 1:
            sys.exit('No test data provided for heat pump performance')
        elif len(self.__dsgn_flow_temps) > 4:
            print( 'Warning: Test data for a maximum of 4 design flow temperatures is expected. '
                 + str(len(self.__dsgn_flow_temps)) + ' have been provided.'
                 )
        for dsgn_flow_temp, data in self.__testdata.items():
            if dupl[dsgn_flow_temp]:
                if (len(data) - dupl[dsgn_flow_temp]) != 4:
                    sys.exit('Expected 4 distinct records for each design flow temperature')
            elif len(data) != 5:
                sys.exit('Expected 5 records for each design flow temperature')

        # Check if test letters ABCDF are present as expected
        test_letter_array = []
        for temperature in self.__dsgn_flow_temps:
            for test_data in self.__testdata[temperature]:
                for test_letter in test_data['test_letter']:
                    test_letter_array.append(test_letter)
                if len(test_letter_array) == 5:
                    for test_letter_check in self.__test_letters_all:
                        if test_letter_check not in test_letter_array:
                            error_output = 'Expected test letter ' + test_letter_check + ' in ' + str(temperature) + ' degree temp data'
                            sys.exit(error_output)
                    test_letter_array = []

        # Sort the list of design flow temps
        self.__dsgn_flow_temps = sorted(self.__dsgn_flow_temps)

        # Sort the records in order of test temperature from low to high
        for dsgn_flow_temp, data in self.__testdata.items():
            data.sort(key=lambda sublist: sublist['temp_test'])

        # Calculate derived variables which are not time-dependent

        def ave_degradation_coeff():
            # The list average_deg_coeff will be in the same order as the
            # corresponding elements in self.__dsgn_flow_temps. This behaviour
            # is relied upon elsewhere.
            average_deg_coeff = []
            for dsgn_flow_temp in self.__dsgn_flow_temps:
                average_deg_coeff.append(
                    sum([
                        x['degradation_coeff']
                        for x in self.__testdata[dsgn_flow_temp]
                        if x['test_letter'] in self.__test_letters_non_bivalent
                        ]) \
                    / len(self.__test_letters_non_bivalent)
                    )
            return average_deg_coeff

        self.__average_deg_coeff = ave_degradation_coeff()

        def ave_capacity():
            # The list average_cap will be in the same order as the
            # corresponding elements in self.__dsgn_flow_temps. This behaviour
            # is relied upon elsewhere.
            average_cap = []
            for dsgn_flow_temp in self.__dsgn_flow_temps:
                average_cap.append(
                    sum([
                        x['capacity']
                        for x in self.__testdata[dsgn_flow_temp]
                        if x['test_letter'] in self.__test_letters_non_bivalent
                        ]) \
                    / len(self.__test_letters_non_bivalent)
                    )
            return average_cap

        self.__average_cap = ave_capacity()

        def init_temp_spread_test_conditions():
            """ List temp spread at test conditions for the design flow temps in the test data """
            dtheta_out_by_flow_temp = {20: 5.0, 35: 5.0, 45: 6.0, 55: 8.0, 65: 10.0}
            dtheta_out = []
            for dsgn_flow_temp in self.__dsgn_flow_temps:
                dtheta_out.append(dtheta_out_by_flow_temp[dsgn_flow_temp])
            return dtheta_out

        self.__temp_spread_test_conditions = init_temp_spread_test_conditions()

        def init_regression_coeffs():
            """ Calculate polynomial regression coefficients for test temperature vs. CoP """
            regression_coeffs = {}
            for dsgn_flow_temp in self.__dsgn_flow_temps:
                temp_test_list = [x['temp_test'] for x in self.__testdata[dsgn_flow_temp]]
                cop_list = [x['cop'] for x in self.__testdata[dsgn_flow_temp]]
                regression_coeffs[dsgn_flow_temp] = (list(polyfit(temp_test_list, cop_list, 2)))

            return regression_coeffs

        self.__regression_coeffs = init_regression_coeffs()

        # Calculate derived variables for each data record which are not time-dependent
        for dsgn_flow_temp in self.__dsgn_flow_temps:
            for data in self.__testdata[dsgn_flow_temp]:
                # Get the source and outlet temperatures from the test record
                temp_source = Celcius2Kelvin(data['temp_source'])
                temp_outlet = Celcius2Kelvin(data['temp_outlet'])

                # Calculate the Carnot CoP and add to the test record
                data['carnot_cop'] = carnot_cop(temp_source, temp_outlet)
                # Calculate the exergetic efficiency and add to the test record
                data['exergetic_eff'] = data['cop'] / data['carnot_cop']

            temp_source_cld = Celcius2Kelvin(self.__testdata[dsgn_flow_temp][0]['temp_source'])
            temp_outlet_cld = Celcius2Kelvin(self.__testdata[dsgn_flow_temp][0]['temp_outlet'])
            carnot_cop_cld = self.__testdata[dsgn_flow_temp][0]['carnot_cop']

            # Calculate derived variables that require values at coldest test temp as inputs
            for data in self.__testdata[dsgn_flow_temp]:
                # Get the source and outlet temperatures from the test record
                temp_source = Celcius2Kelvin(data['temp_source'])
                temp_outlet = Celcius2Kelvin(data['temp_outlet'])

                # Calculate the theoretical load ratio and add to the test record
                data['theoretical_load_ratio'] \
                    = ((data['carnot_cop'] / carnot_cop_cld) \
                    * (temp_outlet_cld * temp_source / (temp_source_cld * temp_outlet)) ** N_EXER)

    def average_degradation_coeff(self, flow_temp):
        """ Return average deg coeff for tests A-D, interpolated between design flow temps """
        if len(self.__dsgn_flow_temps) == 1:
            # If there is data for only one design flow temp, use that
            return self.__average_deg_coeff[0]

        flow_temp = Kelvin2Celcius(flow_temp)
        return np.interp(flow_temp, self.__dsgn_flow_temps, self.__average_deg_coeff)

    def average_capacity(self, flow_temp):
        """ Return average capacity for tests A-D, interpolated between design flow temps """
        if len(self.__dsgn_flow_temps) == 1:
            # If there is data for only one design flow temp, use that
            return self.__average_cap[0]

        flow_temp = Kelvin2Celcius(flow_temp)
        return np.interp(flow_temp, self.__dsgn_flow_temps, self.__average_cap)

    def temp_spread_test_conditions(self, flow_temp):
        """ Return temperature spread under test conditions, interpolated between design flow temps """
        if len(self.__dsgn_flow_temps) == 1:
            # If there is data for only one design flow temp, use that
            return self.__temp_spread_test_conditions[0]

        flow_temp = Kelvin2Celcius(flow_temp)
        return np.interp(flow_temp, self.__dsgn_flow_temps, self.__temp_spread_test_conditions)

    def __find_test_record_index(self, test_condition, dsgn_flow_temp):
        """ Find position of specified test condition in list """
        if test_condition == 'cld':
            # Coldest test condition is first in list
            return 0
        for index, test_record in enumerate(self.__testdata[dsgn_flow_temp]):
            if test_record['test_letter'] == test_condition:
                return index

    def __data_at_test_condition(self, data_item_name, test_condition, flow_temp):
        """ Return value at specified test condition, interpolated between design flow temps """
        # TODO What do we do if flow_temp is outside the range of design flow temps provided?

        if len(self.__dsgn_flow_temps) == 1:
            # If there is data for only one design flow temp, use that
            idx = self.__find_test_record_index(test_condition, self.__dsgn_flow_temps[0])
            return self.__testdata[self.__dsgn_flow_temps[0]][idx][data_item_name]

        # Interpolate between the values at each design flow temp
        data_list = []
        for dsgn_flow_temp in self.__dsgn_flow_temps:
            idx = self.__find_test_record_index(test_condition, dsgn_flow_temp)
            data_list.append(self.__testdata[dsgn_flow_temp][idx][data_item_name])

        flow_temp = Kelvin2Celcius(flow_temp)
        return np.interp(flow_temp, self.__dsgn_flow_temps, data_list)

    def carnot_cop_at_test_condition(self, test_condition, flow_temp):
        """
        Return Carnot CoP at specified test condition (A, B, C, D, F or cld),
        interpolated between design flow temps
        """
        return self.__data_at_test_condition('carnot_cop', test_condition, flow_temp)

    def outlet_temp_at_test_condition(self, test_condition, flow_temp):
        """
        Return outlet temp, in Kelvin, at specified test condition (A, B, C, D,
        F or cld), interpolated between design flow temps.
        """
        return Celcius2Kelvin(
            self.__data_at_test_condition('temp_outlet', test_condition, flow_temp)
            )

    def source_temp_at_test_condition(self, test_condition, flow_temp):
        """
        Return source temp, in Kelvin, at specified test condition (A, B, C, D,
        F or cld), interpolated between design flow temps.
        """
        return Celcius2Kelvin(
            self.__data_at_test_condition('temp_source', test_condition, flow_temp)
            )

    def capacity_at_test_condition(self, test_condition, flow_temp):
        """
        Return capacity, in kW, at specified test condition (A, B, C, D, F or
        cld), interpolated between design flow temps.
        """
        return self.__data_at_test_condition('capacity', test_condition, flow_temp)

    def lr_op_cond(self, flow_temp, temp_source, carnot_cop_op_cond):
        """ Return load ratio at operating conditions """
        lr_op_cond_list = []
        for dsgn_flow_temp in self.__dsgn_flow_temps:
            dsgn_flow_temp = Celcius2Kelvin(dsgn_flow_temp)
            temp_output_cld = self.outlet_temp_at_test_condition('cld', dsgn_flow_temp)
            temp_source_cld = self.source_temp_at_test_condition('cld', dsgn_flow_temp)
            carnot_cop_cld = self.carnot_cop_at_test_condition('cld', dsgn_flow_temp)

            lr_op_cond = (carnot_cop_op_cond / carnot_cop_cld) \
                       * ( temp_output_cld * temp_source
                         / (flow_temp * temp_source_cld)
                         ) \
                         ** N_EXER
            lr_op_cond_list.append(max(1.0, lr_op_cond))

        flow_temp = Kelvin2Celcius(flow_temp)
        return np.interp(flow_temp, self.__dsgn_flow_temps, lr_op_cond_list)

    def lr_eff_degcoeff_either_side_of_op_cond(self, flow_temp, exergy_lr_op_cond):
        """ Return test results either side of operating conditions.

        This function returns 6 results:
        - Exergy load ratio below operating conditions
        - Exergy load ratio above operating conditions
        - Exergy efficiency below operating conditions
        - Exergy efficiency above operating conditions
        - Degradation coeff below operating conditions
        - Degradation coeff above operating conditions

        Arguments:
        flow_temp         -- flow temperature, in Kelvin
        exergy_lr_op_cond -- exergy load ratio at operating conditions
        """
        load_ratios_below = []
        load_ratios_above = []
        efficiencies_below = []
        efficiencies_above = []
        degradation_coeffs_below = []
        degradation_coeffs_above = []

        # For each design flow temperature, find load ratios in test data
        # either side of load ratio calculated for operating conditions.
        # Note: Loop over sorted list of design flow temps and then index into
        #       self.__testdata, rather than looping over self.__testdata,
        #       which is unsorted and therefore may populate the lists in the
        #       wrong order.
        for dsgn_flow_temp in self.__dsgn_flow_temps:
            found = False
            dsgn_flow_temp_data = self.__testdata[dsgn_flow_temp]
            # Find the first load ratio in the test data that is greater than
            # or equal to than the load ratio at operating conditions - this
            # and the previous load ratio are the values either side of
            # operating conditions.
            for idx, test_record in enumerate(dsgn_flow_temp_data):
                # Note: Changed the condition below from ">=" to ">" because
                # otherwise when exergy_lr_op_cond == test_record['theoretical_load_ratio']
                # for the first record, idx == 0 which is not allowed
                if test_record['theoretical_load_ratio'] > exergy_lr_op_cond:
                    assert idx > 0
                    found = True
                    # Current value of idx will be used later, so break out of loop
                    break

            if not found:
                # Use the highest (list index -1) and second highest
                idx = -1

            # Look up correct load ratio and efficiency based on the idx found above
            load_ratios_below.append(dsgn_flow_temp_data[idx-1]['theoretical_load_ratio'])
            load_ratios_above.append(dsgn_flow_temp_data[idx]['theoretical_load_ratio'])
            efficiencies_below.append(dsgn_flow_temp_data[idx-1]['exergetic_eff'])
            efficiencies_above.append(dsgn_flow_temp_data[idx]['exergetic_eff'])
            degradation_coeffs_below.append(dsgn_flow_temp_data[idx-1]['degradation_coeff'])
            degradation_coeffs_above.append(dsgn_flow_temp_data[idx]['degradation_coeff'])

        if len(self.__dsgn_flow_temps) == 1:
            # If there is data for only one design flow temp, use that
            return load_ratios_below[0], load_ratios_above[0], \
                   efficiencies_below[0], efficiencies_above[0], \
                   degradation_coeffs_below[0], degradation_coeffs_above[0]

        # Interpolate between the values found for the different design flow temperatures
        flow_temp = Kelvin2Celcius(flow_temp)
        lr_below = np.interp(flow_temp, self.__dsgn_flow_temps, load_ratios_below)
        lr_above = np.interp(flow_temp, self.__dsgn_flow_temps, load_ratios_above)
        eff_below = np.interp(flow_temp, self.__dsgn_flow_temps, efficiencies_below)
        eff_above = np.interp(flow_temp, self.__dsgn_flow_temps, efficiencies_above)
        deg_below = np.interp(flow_temp, self.__dsgn_flow_temps, degradation_coeffs_below)
        deg_above = np.interp(flow_temp, self.__dsgn_flow_temps, degradation_coeffs_above)

        return lr_below, lr_above, eff_below, eff_above, deg_below, deg_above

    def cop_op_cond_if_not_air_source(
            self,
            temp_diff_limit_low,
            temp_ext_C,
            temp_source,
            temp_output,
            ):
        """ Calculate CoP at operating conditions when heat pump is not air-source

        Arguments:
        temp_diff_limit_low -- minimum temperature difference between source and sink
        temp_ext_C           -- external temperature, in Celsius. 
            Need to use Celsius for temp_ext_C because regression coeffs were calculated
            using temperature in Celsius.
        temp_source        -- source temperature, in Kelvin
        temp_output        -- output temperature, in Kelvin
        """

        # For each design flow temperature, calculate CoP at operating conditions
        # Note: Loop over sorted list of design flow temps and then index into
        #       self.__testdata, rather than looping over self.__testdata,
        #       which is unsorted and therefore may populate the lists in the
        #       wrong order.
        cop_op_cond = []
        for dsgn_flow_temp in self.__dsgn_flow_temps:
            dsgn_flow_temp_data = self.__testdata[dsgn_flow_temp]
            # Get the source and outlet temperatures from the coldest test record
            temp_outlet_cld = Celcius2Kelvin(dsgn_flow_temp_data[0]['temp_outlet'])
            temp_source_cld = Celcius2Kelvin(dsgn_flow_temp_data[0]['temp_source'])

            cop_operating_conditions \
                = ( self.__regression_coeffs[dsgn_flow_temp][0] \
                  + self.__regression_coeffs[dsgn_flow_temp][1] * temp_ext_C \
                  + self.__regression_coeffs[dsgn_flow_temp][2] * temp_ext_C ** 2 \
                  ) \
                * temp_output * (temp_outlet_cld - temp_source_cld) \
                / ( temp_outlet_cld * max( (temp_output - temp_source), temp_diff_limit_low))
            cop_op_cond.append(cop_operating_conditions)

        if len(self.__dsgn_flow_temps) == 1:
            # If there is data for only one design flow temp, use that
            return cop_op_cond[0]

        # Interpolate between the values found for the different design flow temperatures
        flow_temp = Kelvin2Celcius(temp_output)
        return np.interp(flow_temp, self.__dsgn_flow_temps, cop_op_cond)

    def capacity_op_cond_var_flow_or_source_temp(self, temp_output, temp_source, mod_ctrl):
        """ Calculate thermal capacity at operating conditions when flow temp
        during test was variable or source temp was variable
        
        Arguments:
        temp_source -- source temperature, in Kelvin
        temp_output -- output temperature, in Kelvin
        mod_ctrl -- boolean specifying whether or not the heat has controls
                    capable of varying the output (as opposed to just on/off
                    control)
        """
        # In eqns below, method uses condition A rather than coldest. From
        # CALCM-01 - DAHPSE - V2.0_DRAFT13, section 4.4:
        # The Temperature Operation Limit (TOL) is defined in EN14825 as
        # "the lowest outdoor temperature at which the unit can still
        # deliver heating capacity and is declared by the manufacturer.
        # Below this temperature the heat pump will not be able to
        # deliver any heating capacity."
        # The weather data used within this calculation method does not
        # feature a source temperature at or below the "TOL" test
        # temperature (which is -7C to -10C). Therefore, test data at
        # the TOL test condition is not used (Test condition "A" at -7C
        # is sufficient).
        # TODO The above implies that the TOL test temperature data may
        #      be needed if we change the weather data from that used in
        #      DAHPSE for SAP 2012/10.2
        therm_cap_op_cond = []

        if mod_ctrl:
            # For each design flow temperature, calculate capacity at operating conditions
            # Note: Loop over sorted list of design flow temps and then index into
            #       self.__testdata, rather than looping over self.__testdata,
            #       which is unsorted and therefore may populate the lists in the
            #       wrong order.
            for dsgn_flow_temp in self.__dsgn_flow_temps:
                dsgn_flow_temp_data = self.__testdata[dsgn_flow_temp]
                # Get the source and outlet temperatures from the coldest test record
                temp_outlet_cld = Celcius2Kelvin(dsgn_flow_temp_data[0]['temp_outlet'])
                temp_source_cld = Celcius2Kelvin(dsgn_flow_temp_data[0]['temp_source'])
                # Get the thermal capacity from the coldest test record
                thermal_capacity_cld = dsgn_flow_temp_data[0]['capacity']

                thermal_capacity_op_cond \
                    = thermal_capacity_cld \
                    * ( (temp_outlet_cld * temp_source) \
                      / (temp_output * temp_source_cld) \
                      ) \
                    ** N_EXER
                therm_cap_op_cond.append(thermal_capacity_op_cond)
        else:
            # For each design flow temperature, calculate capacity at operating conditions
            # Note: Loop over sorted list of design flow temps and then index into
            #       self.__testdata, rather than looping over self.__testdata,
            #       which is unsorted and therefore may populate the lists in the
            #       wrong order.
            for dsgn_flow_temp in self.__dsgn_flow_temps:
                dsgn_flow_temp_data = self.__testdata[dsgn_flow_temp]
                # Get the source and outlet temperatures from the coldest test record
                temp_outlet_cld = Celcius2Kelvin(dsgn_flow_temp_data[0]['temp_outlet'])
                temp_source_cld = Celcius2Kelvin(dsgn_flow_temp_data[0]['temp_source'])
                # Get the thermal capacity from the coldest test record
                thermal_capacity_cld = dsgn_flow_temp_data[0]['capacity']

                D_idx = self.__find_test_record_index('D', dsgn_flow_temp)
                # Get the source and outlet temperatures for test condition D
                temp_outlet_D = Celcius2Kelvin(dsgn_flow_temp_data[D_idx]['temp_outlet'])
                temp_source_D = Celcius2Kelvin(dsgn_flow_temp_data[D_idx]['temp_source'])
                # Get the thermal capacity for test condition D
                thermal_capacity_D = dsgn_flow_temp_data[D_idx]['capacity']

                temp_diff_cld = temp_outlet_cld - temp_source_cld
                temp_diff_D = temp_outlet_D - temp_source_D
                temp_diff_op_cond = temp_output - temp_source

                thermal_capacity_op_cond \
                    = thermal_capacity_cld \
                    + (thermal_capacity_D - thermal_capacity_cld) \
                    * ( (temp_diff_cld - temp_diff_op_cond) \
                      / (temp_diff_cld - temp_diff_D) \
                      )
                therm_cap_op_cond.append(thermal_capacity_op_cond)

        # Interpolate between the values found for the different design flow temperatures
        flow_temp = Kelvin2Celcius(temp_output)
        return np.interp(flow_temp, self.__dsgn_flow_temps, therm_cap_op_cond)

    def temp_spread_correction(
            self,
            temp_source,
            temp_output,
            temp_diff_evaporator,
            temp_diff_condenser,
            temp_spread_emitter,
            ):
        """ Calculate temperature spread correction factor

        Arguments:
        temp_source -- source temperature, in Kelvin
        temp_output -- output temperature, in Kelvin
        temp_diff_evaporator
            -- average temperature difference between heat transfer medium and
               refrigerant in evaporator, in deg C or Kelvin
        temp_diff_condenser
            -- average temperature difference between heat transfer medium and
               refrigerant in condenser, in deg C or Kelvin
        temp_spread_emitter
            -- temperature spread on condenser side in operation due to design
               of heat emission system
        """
        temp_spread_correction_list = []
        for i, dsgn_flow_temp in enumerate(self.__dsgn_flow_temps):
            temp_spread_test_cond = self.__temp_spread_test_conditions[i]
            temp_spread_correction \
                = 1.0 \
                - ((temp_spread_test_cond - temp_spread_emitter) / 2.0) \
                / ( temp_output - temp_spread_test_cond / 2.0 + temp_diff_condenser \
                  - temp_source + temp_diff_evaporator
                  )
            temp_spread_correction_list.append(temp_spread_correction)

        # Interpolate between the values found for the different design flow temperatures
        flow_temp = Kelvin2Celcius(temp_output)
        return np.interp(flow_temp, self.__dsgn_flow_temps, temp_spread_correction_list)


class HeatPumpService:
    """ A base class for objects representing services (e.g. water heating) provided by a heat pump.

    This object encapsulates the name of the service, meaning that the system
    consuming the energy does not have to specify this on every call, and
    helping to enforce that each service has a unique name.

    Derived objects provide a place to handle parts of the calculation (e.g.
    distribution flow temperature) that may differ for different services.

    Separate subclasses need to be implemented for different types of service
    (e.g. HW and space heating). These should implement the following functions:
    - demand_energy(self, energy_demand)
    """

    def __init__(self, heat_pump, service_name, control=None):
        """ Construct a HeatPumpService object

        Arguments:
        heat_pump    -- reference to the HeatPump object providing the service
        service_name -- name of the service demanding energy from the heat pump
        control -- reference to a control object which must implement is_on() func
        """
        self.__hp = heat_pump
        self.__service_name = service_name
        self.__control = control

    def is_on(self):
        if self.__control is not None:
            service_on = self.__control.is_on()
        else:
            service_on = True
        return service_on



class HeatPumpServiceWater(HeatPumpService):
    """ An object to represent a water heating service provided by a heat pump to e.g. a cylinder.

    This object contains the parts of the heat pump calculation that are
    specific to providing hot water.
    """

    __TIME_CONSTANT_WATER = 1560
    __SERVICE_TYPE = ServiceType.WATER

    def __init__(
            self,
            heat_pump,
            service_name,
            temp_hot_water,
            temp_limit_upper,
            cold_feed,
            control=None,
            boiler_service_water_regular=None,
            ):
        """ Construct a BoilerServiceWater object

        Arguments:
        heat_pump -- reference to the HeatPump object providing the service
        service_name -- name of the service demanding energy from the heat pump
        temp_hot_water -- temperature of the hot water to be provided, in deg C
        temp_limit_upper -- upper operating limit for temperature, in deg C
        cold_feed -- reference to ColdWaterSource object
        control -- reference to a control object which must implement is_on() func
        boiler_service_water_regular -- reference to the BoilerServiceWaterRegular
            object that is backing up or supplementing the heat pump service
        """
        super().__init__(heat_pump, service_name, control)

        self.__temp_hot_water = Celcius2Kelvin(temp_hot_water)
        # TODO Should temp_return_feed be calculated per timestep?
        self.__temp_limit_upper = Celcius2Kelvin(temp_limit_upper)
        self.__cold_feed = cold_feed
        self.__hybrid_boiler_service = boiler_service_water_regular

    def energy_output_max(self, temp_return):
        """ Calculate the maximum energy output of the HP, accounting for time
            spent on higher-priority services
        """
        temp_return_K = Celcius2Kelvin(temp_return)
        if not self.is_on():
            return 0.0

        return self._HeatPumpService__hp._HeatPump__energy_output_max(
            self.__temp_hot_water,
            temp_return_K,
            self.__hybrid_boiler_service,
            self.__SERVICE_TYPE,
            )

    def demand_energy(self, energy_demand, temp_return):
        """ Demand energy (in kWh) from the heat pump """
        temp_cold_water = Celcius2Kelvin(self.__cold_feed.temperature())
        temp_return_K = Celcius2Kelvin(temp_return)
        
        service_on = self.is_on()
        if not service_on:
            energy_demand = 0.0

        return self._HeatPumpService__hp._HeatPump__demand_energy(
            self._HeatPumpService__service_name,
            self.__SERVICE_TYPE,
            energy_demand,
            self.__temp_hot_water,
            temp_return_K,
            self.__temp_limit_upper,
            self.__TIME_CONSTANT_WATER,
            service_on,
            temp_used_for_scaling = temp_cold_water,
            hybrid_boiler_service = self.__hybrid_boiler_service,
            )


class HeatPumpServiceSpace(HeatPumpService):
    """ An object to represent a space heating service provided by a heat pump to e.g. radiators.

    This object contains the parts of the heat pump calculation that are
    specific to providing space heating.
    """

    __TIME_CONSTANT_SPACE = {
        SinkType.WATER: 1370,
        SinkType.AIR: 120,
        }
    __SERVICE_TYPE = ServiceType.SPACE

    def __init__(
            self,
            heat_pump,
            service_name,
            temp_limit_upper,
            temp_diff_emit_dsgn,
            control,
            volume_heated,
            boiler_service_space=None,
            ):
        """ Construct a BoilerServiceSpace object

        Arguments:
        heat_pump -- reference to the HeatPump object providing the service
        service_name -- name of the service demanding energy from the heat pump
        temp_limit_upper -- upper operating limit for temperature, in deg C
        temp_diff_emit_dsgn -- design temperature difference across the emitters, in deg C or K
        control -- reference to a control object which must implement is_on() and setpnt() funcs
        boiler_service_space -- reference to the BoilerServiceWaterSpace
            object that is backing up or supplementing the heat pump service
        volume_heated -- volume of zones heated (required for exhaust air HPs only), in m3
        """
        super().__init__(heat_pump, service_name, control)
        self.__temp_limit_upper = Celcius2Kelvin(temp_limit_upper)
        self.__temp_diff_emit_dsgn = temp_diff_emit_dsgn
        self.__hybrid_boiler_service = boiler_service_space
        self.__volume_heated = volume_heated
        
    def temp_setpnt(self):
        return self._HeatPumpService__control.setpnt()

    def in_required_period(self):
        return self._HeatPumpService__control.in_required_period()
    
    def energy_output_max(self, temp_output, temp_return_feed, emitters_data_for_buffer_tank=None):
        """ Calculate the maximum energy output of the HP, accounting for time
            spent on higher-priority services
        """
        if not self.is_on():
            return 0.0

        temp_output = Celcius2Kelvin(temp_output)
        return self._HeatPumpService__hp._HeatPump__energy_output_max(
            temp_output,
            temp_return_feed,
            self.__hybrid_boiler_service,
            self.__SERVICE_TYPE,
            temp_spread_correction = self.temp_spread_correction,
            emitters_data_for_buffer_tank = emitters_data_for_buffer_tank,
            service_name = self._HeatPumpService__service_name
            )
    
    def demand_energy(self, energy_demand, temp_flow, temp_return, emitters_data_for_buffer_tank=None):
        """ Demand energy (in kWh) from the heat pump

        Arguments:
        energy_demand -- space heating energy demand, in kWh
        temp_flow -- flow temperature for emitters, in deg C
        temp_return -- return temperature for emitters, in deg C
        """
        service_on = self.is_on()
        if not service_on:
            energy_demand = 0.0

        return self._HeatPumpService__hp._HeatPump__demand_energy(
            self._HeatPumpService__service_name,
            self.__SERVICE_TYPE,
            energy_demand,
            Celcius2Kelvin(temp_flow),
            Celcius2Kelvin(temp_return),
            self.__temp_limit_upper,
            self.__TIME_CONSTANT_SPACE[self._HeatPumpService__hp._HeatPump__sink_type],
            service_on,
            temp_spread_correction = self.temp_spread_correction,
            hybrid_boiler_service = self.__hybrid_boiler_service,
            emitters_data_for_buffer_tank = emitters_data_for_buffer_tank, 
            )

    def running_time_throughput_factor(
            self,
            space_heat_running_time_cumulative,
            energy_demand,
            temp_flow,
            temp_return,
            ):
        """ Return the cumulative running time and throughput factor (exhaust air HPs only) """
        service_on = self.is_on()
        if not service_on:
            energy_demand = 0.0

        return self._HeatPumpService__hp._HeatPump__running_time_throughput_factor(
            space_heat_running_time_cumulative,
            self._HeatPumpService__service_name,
            ServiceType.SPACE,
            energy_demand,
            Celcius2Kelvin(temp_flow),
            Celcius2Kelvin(temp_return),
            self.__temp_limit_upper,
            self.__TIME_CONSTANT_SPACE[self._HeatPumpService__hp._HeatPump__sink_type],
            service_on,
            self.__volume_heated,
            temp_spread_correction = self.temp_spread_correction,
            )

    def temp_spread_correction(self, temp_output, temp_source):
        """Calculate temperature spread correction """
        # Average temperature difference between heat transfer medium and
        # refrigerant in condenser
        temp_diff_condenser = 5.0

        # Average temperature difference between heat transfer medium and
        # refrigerant in evaporator
        # TODO Figures in BS EN ISO 15316-4-2:2017 are -15 and -10, but figures
        #      in BS EN ISO 15316-4-2:2008 were positive (although some were
        #      different numbers) and signs in temp_spread_correction equation
        #      have not changed, so need to check which is correct.
        #      Note: using negative numbers leads to divide-by-zero errors in
        #      the calculation which do not occur when using positive numbers.
        #      Given that the equation that uses these figures already has a
        #      minus sign in front of this variable (as written in the standard)
        #      this would seem to suggest that using positive numbers is correct
        if SourceType.source_fluid_is_air(self._HeatPumpService__hp._HeatPump__source_type):
            temp_diff_evaporator = 15.0
        elif SourceType.source_fluid_is_water(self._HeatPumpService__hp._HeatPump__source_type):
            temp_diff_evaporator = 10.0
        else:
            sys.exit('SourceType not recognised')

        # TODO The temp_spread_emitter input below (self.__temp_diff_emit_dsgn)
        #      is for no weather comp. Add weather comp case as well
        return self._HeatPumpService__hp._HeatPump__test_data.temp_spread_correction(
            temp_source,
            temp_output,
            temp_diff_evaporator,
            temp_diff_condenser,
            self.__temp_diff_emit_dsgn,
            )


class HeatPumpServiceSpaceWarmAir(HeatPumpServiceSpace):
    """ An object to represent a warm air space heating service provided by a heat pump.

    This object contains the parts of the heat pump calculation that are
    specific to providing space heating via warm air.
    """
    def __init__(
            self,
            heat_pump,
            service_name,
            temp_diff_emit_dsgn,
            control,
            temp_flow,
            frac_convective,
            volume_heated,
            ):
        """ Construct a HeatPumpServiceSpaceWarmAir object

        Arguments:
        heat_pump -- reference to the HeatPump object providing the service
        service_name -- name of the service demanding energy from the heat pump
        temp_diff_emit_dsgn -- design temperature difference across the emitters, in deg C or K
        control -- reference to a control object which must implement is_on() and setpnt() funcs
        temp_flow -- flow temperature, in deg C
        frac_convective -- convective fraction for heating
        volume_heated -- volume of zones heated (required for exhaust air HPs only), in m3
        """
        self.__frac_convective = frac_convective
        self.__temp_flow = temp_flow
        # Return temp won't be used in the relevant code paths anyway, so this is arbitrary
        self.__temp_return = temp_flow

        # Upper operating limit for temperature, in deg C
        temp_limit_upper = temp_flow

        super().__init__(
            heat_pump,
            service_name,
            temp_limit_upper,
            temp_diff_emit_dsgn,
            control,
            volume_heated,
            )

    def demand_energy(self, energy_demand):
        """ Demand energy (in kWh) from the heat pump

        Arguments:
        energy_demand -- space heating energy demand, in kWh
        """
        return HeatPumpServiceSpace.demand_energy(
            self,
            energy_demand,
            self.__temp_flow,
            self.__temp_return,
            )

    def running_time_throughput_factor(self, energy_demand, space_heat_running_time_cumulative):
        """ Return the cumulative running time and throughput factor (exhaust air HPs only)

        Arguments:
        energy_demand -- in kWh
        space_heat_running_time_cumulative
            -- running time spent on higher-priority space heating services
        """
        return HeatPumpServiceSpace.running_time_throughput_factor(
            self,
            space_heat_running_time_cumulative,
            energy_demand,
            self.__temp_flow,
            self.__temp_return,
            )

    def frac_convective(self):
        return self.__frac_convective


class HeatPump:
    """ An object to represent an electric heat pump """

    # From CALCM-01 - DAHPSE - V2.0_DRAFT13, section 4.5.3:
    # A minimum temperature difference of 6K between the source and sink
    # temperature is applied to prevent very high Carnot COPs entering the
    # calculation. This only arises when the temperature difference and heating
    # load is small and is unlikely to affect the calculated SPF.
    __temp_diff_limit_low = 6.0 # Kelvin

    # Fraction of the energy input dedicated to auxiliaries when on
    # TODO This is always zero for electric heat pumps, but if we want to deal
    #      with non-electric heat pumps then this will need to be altered.
    __f_aux = 0.0

    def __init__(
            self,
            hp_dict,
            energy_supply,
            energy_supply_conn_name_auxiliary,
            simulation_time,
            external_conditions,
            number_of_zones,
            throughput_exhaust_air=None,
            heat_network=None,
            output_detailed_results=False,
            boiler=None,
            cost_schedule_hybrid_hp=None,
            project=None,
            ):
        """ Construct a HeatPump object

        Arguments:
        hp_dict -- dictionary of heat pump characteristics, with the following elements:
            - test_data -- EN 14825 test data (list of dictionaries)
            - SourceType -- string specifying heat source type, one of:
                - "Ground"
                - "OutsideAir"
                - "ExhaustAirMEV"
                - "ExhaustAirMVHR"
                - "ExhaustAirMixed"
                - "WaterGround"
                - "WaterSurface"
                - "HeatNetwork"
            - SinkType -- string specifying heat distribution type, one of:
                - "Air"
                - "Water"
            - BackupCtrlType -- string specifying control arrangement for backup
                                heater, one of:
                - "None" -- backup heater disabled or not present
                - "TopUp" -- when heat pump has insufficient capacity, backup
                             heater will supplement the heat pump
                - "Substitute" -- when heat pump has insufficient capacity, backup
                                  heater will provide all the heat required, and
                                  heat pump will switch off
            - time_delay_backup -- time after which the backup heater will activate
                                   if demand has not been satisfied
            - modulating_control -- boolean specifying whether or not the heat
                                    has controls capable of varying the output
                                    (as opposed to just on/off control)
            - time_constant_onoff_operation
                -- a characteristic parameter of the heat pump, due to the
                   inertia of the on/off transient
            - temp_return_feed_max -- maximum allowable temperature of the
                                      return feed, in Celsius
            - temp_lower_operating_limit
                -- minimum source temperature at which the heat pump can operate,
                   in Celsius
            - min_temp_diff_flow_return_for_hp_to_operate
                -- minimum difference between flow and return temperatures
                   required for the HP to operate, in Celsius or Kelvin
            - var_flow_temp_ctrl_during_test
                -- boolean specifying whether or not variable flow temperature
                   control was enabled during the EN 14825 tests
            - power_heating_circ_pump -- power (kW) of central heating circulation pump
            - power_source_circ_pump
                -- power (kW) of source ciculation pump or fan circulation when not
                   implicit in CoP measurements
            - power_standby -- power (kW) consumption in standby mode
            - power_crankcase_heater -- power (kW) consumption in crankcase heater mode
            - power_off -- power (kW) consumption in off mode
            - power_max_backup -- max. power (kW) of backup heater
            - temp_distribution_heat_network
                  -- distribution temperature of the heat network (for HPs that use heat
                     network as heat source)
        energy_supply -- reference to EnergySupply object
        energy_supply_conn_name_auxiliary
            -- name to be used for EnergySupplyConnection object for auxiliary energy
        simulation_time -- reference to SimulationTime object
        external_conditions -- reference to ExternalConditions object
        number_of_zones  -- number of zones in the project
        throughput_exhaust_air -- throughput (m3/h) of exhaust air
        heat_network -- reference to EnergySupply object representing heat network
                        (for HPs that use heat network as heat source)
        output_detailed_results -- if true, save detailed results from each timestep
                                   for later reporting
        boiler -- reference to Boiler object representing boiler
                  (for hybrid heat pumps)
        cost_schedule_hybrid_hp -- dict of cost schedule for heat pump and boiler

        Other variables:
        energy_supply_connections
            -- dictionary with service name strings as keys and corresponding
               EnergySupplyConnection objects as values
        energy_supply_connection_aux -- EnergySupplyConnection object for auxiliary energy
        test_data -- HeatPumpTestData object
        """
        self.__energy_supply = energy_supply
        self.__simulation_time = simulation_time
        self.__external_conditions = external_conditions
        self.__throughput_exhaust_air = throughput_exhaust_air

        self.__energy_supply_connections = {}
        self.__energy_supply_connection_aux \
            = self.__energy_supply.connection(energy_supply_conn_name_auxiliary)

        self.__service_results = []
        self.__total_time_running_current_timestep = 0.0
        self.__time_running_continuous = 0.0

        # Assign hp_dict elements to member variables of this class
        self.__source_type = SourceType.from_string(hp_dict['source_type'])
        self.__sink_type = SinkType.from_string(hp_dict['sink_type'])
        self.__backup_ctrl = BackupCtrlType.from_string(hp_dict['backup_ctrl_type'])
        if self.__backup_ctrl != BackupCtrlType.NONE:
            self.__time_delay_backup = float(hp_dict['time_delay_backup'])
            if boiler is None:
                self.__power_max_backup = hp_dict['power_max_backup']
            else:
                self.__power_max_backup = 0
        else:
            self.__power_max_backup = 0
        self.__modulating_ctrl = bool(hp_dict['modulating_control'])
        self.__time_constant_onoff_operation = float(hp_dict['time_constant_onoff_operation'])
        if self.__sink_type != SinkType.AIR:
            self.__temp_return_feed_max = Celcius2Kelvin(float(hp_dict['temp_return_feed_max']))
        self.__temp_lower_op_limit = Celcius2Kelvin(float(hp_dict['temp_lower_operating_limit']))
        self.__temp_diff_flow_return_min \
            = float(hp_dict['min_temp_diff_flow_return_for_hp_to_operate'])
        self.__var_flow_temp_ctrl_during_test = bool(hp_dict['var_flow_temp_ctrl_during_test'])
        self.__power_heating_circ_pump = 0.0 if 'power_heating_circ_pump' not in hp_dict else hp_dict['power_heating_circ_pump']
        self.__power_heating_warm_air_fan = 0.0 if 'power_heating_warm_air_fan' not in hp_dict else hp_dict['power_heating_warm_air_fan']

        self.__power_source_circ_pump = hp_dict['power_source_circ_pump']
        self.__power_standby = hp_dict['power_standby']
        self.__power_crankcase_heater_mode = hp_dict['power_crankcase_heater']
        self.__power_off_mode = hp_dict['power_off']
        self.__boiler = boiler
        self.__cost_schedule_hybrid_hp = cost_schedule_hybrid_hp
        if cost_schedule_hybrid_hp is not None:
            self.__cost_schedule_hp = expand_schedule(float, cost_schedule_hybrid_hp["cost_schedule_hp"], "main", False) 
            self.__cost_schedule_boiler = expand_schedule(float, cost_schedule_hybrid_hp["cost_schedule_boiler"], "main", False)
            self.__cost_schedule_start_day = cost_schedule_hybrid_hp["cost_schedule_start_day"]
            self.__cost_schedule_time_series_step = cost_schedule_hybrid_hp["cost_schedule_time_series_step"]
        self.__project = project

        # HPs that use heat network as heat source require different/additional
        # initialisation, which is implemented here
        if self.__source_type == SourceType.HEAT_NETWORK:
            if heat_network is None:
                sys.exit('If HP uses heat network as source, then heat network must be specified')
            self.__temp_distribution_heat_network = hp_dict['temp_distribution_heat_network']
            self.__energy_supply_HN = heat_network
            self.__energy_supply_HN_connections = {}

        # Exhaust air HP requires different/additional initialisation, which is implemented here
        if SourceType.is_exhaust_air(self.__source_type):
            lowest_air_flow_rate_in_test_data, hp_dict['test_data'] \
                = interpolate_exhaust_air_heat_pump_test_data(
                    throughput_exhaust_air,
                    hp_dict['test_data'],
                    self.__source_type,
                    )
            self.__overvent_ratio = max(
                1.0,
                lowest_air_flow_rate_in_test_data / throughput_exhaust_air,
                )
            self.__volume_heated_all_services = 0.0
        else:
            self.__overvent_ratio = 1.0
        # TODO For now, disable exhaust air heat pump when conditions are out of
        #      range of test data. Decision to be made on whether to allow this
        #      and if so how to handle it in the calculation
        if self.__overvent_ratio > 1.0:
            sys.exit("Exhaust air heat pump: Flow rate for associated ventilation "
                     "system is lower than flow rate in any of the test data provided.")

        # Check there is no remaining test data specific to an air flow rate
        # For exhaust air HPs, this should have been eliminated in the
        # interpolation above and for other HPs, it should not be present in the
        # first place.
        for test_data_record in hp_dict['test_data']:
            if 'air_flow_rate' in test_data_record:
                sys.exit('Unexpected test data specific to an air flow rate')

        # Parse and initialise heat pump test data
        self.__test_data = HeatPumpTestData(hp_dict['test_data'])

        # Mixed exhaust air heat pumps
        if self.__source_type == SourceType.EXHAUST_AIR_MIXED:
            self.__eahp_mixed_max_temp = hp_dict['eahp_mixed_max_temp']
            self.__eahp_mixed_min_temp = hp_dict['eahp_mixed_min_temp']
            ext_air_ratio_list = [test_data['ext_air_ratio'] for test_data in hp_dict['test_data']]
            same_ext_air_ratio = len(set(ext_air_ratio_list))
            if same_ext_air_ratio > 1:
                sys.exit('Error: More than one unique external air ratio entered.')
            self.__ext_air_ratio = sum(ext_air_ratio_list) / len(ext_air_ratio_list)

        if self.__modulating_ctrl:
            if self.__sink_type == SinkType.AIR:
                self.__temp_min_modulation_rate_low = 20.0
                self.__min_modulation_rate_low = float(hp_dict['min_modulation_rate_20'])
            elif self.__sink_type == SinkType.WATER:
                self.__temp_min_modulation_rate_low = 35.0
                self.__min_modulation_rate_low = float(hp_dict['min_modulation_rate_35'])
            else:
                sys.exit('Sink type not recognised')

            if 55.0 in self.__test_data._HeatPumpTestData__dsgn_flow_temps:
                self.__temp_min_modulation_rate_high = 55.0
                self.__min_modulation_rate_55 = float(hp_dict['min_modulation_rate_55'])

        # If detailed results are to be output, initialise list
        if output_detailed_results:
            self.__detailed_results = []
        else:
            self.__detailed_results = None
            
        # Add Buffer Tank object to HP if defined in input file    
        if 'BufferTank' in hp_dict :
            self.__buffer_tank = BufferTank (
                hp_dict['BufferTank']['daily_losses'],
                hp_dict['BufferTank']['volume'],
                hp_dict['BufferTank']['pump_fixed_flow_rate'],
                hp_dict['BufferTank']['pump_power_at_flow_rate'],
                number_of_zones,
                self.__simulation_time,
                WATER,
                output_detailed_results,
                )
        else:
            self.__buffer_tank = None
        
    def source_is_exhaust_air(self):
        return SourceType.is_exhaust_air(self.__source_type)

    def buffer_int_gains(self):
        if self.__buffer_tank:
            return self.__buffer_tank.internal_gains()
        else:
            return 0.0

    def __create_service_connection(self, service_name):
        """ Return a HeatPumpService object """
        # Check that service_name is not already registered
        if service_name in self.__energy_supply_connections.keys():
            sys.exit("Error: Service name already used: "+service_name)
            # TODO Exit just the current case instead of whole program entirely?

        # Set up EnergySupplyConnection for this service
        self.__energy_supply_connections[service_name] = \
            self.__energy_supply.connection(service_name)

        # If HP uses heat network as source, then set up connection
        if self.__source_type == SourceType.HEAT_NETWORK:
            self.__energy_supply_HN_connections[service_name] \
                = self.__energy_supply_HN.connection(service_name)

    def create_service_hot_water_combi(
            self,
            boiler_data,
            service_name,
            temp_hot_water,
            cold_feed,
            ):

        if self.__boiler == None:
            sys.exit('Error: Missing Boiler object for hybrid heat pump. '
                     'Non-hybrid heat pumps cannot be used as an '
                     'instantaneous hot water system')

        return self.__boiler.create_service_hot_water_combi(
                boiler_data,
                service_name,
                temp_hot_water,
                cold_feed,
                )

    def create_service_hot_water(
            self,
            service_name,
            temp_hot_water,
            temp_limit_upper,
            cold_feed,
            control=None,
            ):
        """ Return a HeatPumpServiceWater object and create an EnergySupplyConnection for it

        Arguments:
        service_name -- name of the service demanding energy from the boiler
        temp_hot_water -- temperature of the hot water to be provided, in deg C
        temp_limit_upper -- upper operating limit for temperature, in deg C
        cold_feed -- reference to ColdWaterSource object
        control -- reference to a control object which must implement is_on() func
        """
        self.__create_service_connection(service_name)

        if self.__boiler is not None:
            boiler_service = self.__boiler.create_service_hot_water_regular(
                service_name,
                temp_hot_water,
                cold_feed,
                control
                )
        else:
            boiler_service = None

        return HeatPumpServiceWater(
            self,
            service_name,
            temp_hot_water,
            temp_limit_upper,
            cold_feed,
            control,
            boiler_service,
            )

    def create_service_space_heating(
            self,
            service_name,
            temp_limit_upper,
            temp_diff_emit_dsgn,
            control,
            volume_heated,
            ):
        """ Return a HeatPumpServiceSpace object and create an EnergySupplyConnection for it

        Arguments:
        service_name -- name of the service demanding energy from the heat pump
        temp_limit_upper -- upper operating limit for temperature, in deg C
        temp_diff_emit_dsgn -- design temperature difference across the emitters, in deg C or K
        control -- reference to a control object which must implement is_on() func
        volume_heated -- volume of zones heated (required for exhaust air HPs only), in m3
        """
        if self.__boiler is not None:
            boiler_service = self.__boiler.create_service_space_heating(service_name, control)
        else:
            boiler_service = None
        if self.source_is_exhaust_air():
            self.__volume_heated_all_services += volume_heated
        self.__create_service_connection(service_name)
        return HeatPumpServiceSpace(
            self,
            service_name,
            temp_limit_upper,
            temp_diff_emit_dsgn,
            control,
            volume_heated,
            boiler_service,
            )

    def create_service_space_heating_warm_air(
            self,
            service_name,
            control,
            frac_convective,
            volume_heated,
            ):
        """ Return a HeatPumpServiceSpaceWarmAir object and create an EnergySupplyConnection for it

        Arguments:
        service_name -- name of the service demanding energy from the heat pump
        control -- reference to a control object which must implement is_on() func
        frac_convective -- convective fraction for heating
        volume_heated -- volume of zones heated (required for exhaust air HPs only), in m3
        """
        if self.__sink_type != SinkType.AIR:
            sys.exit('Warm air space heating service requires heat pump with sink type Air')

        if self.__boiler is not None:
            # TODO More evidence is required before warm air space heating systems can work with hybrid heat pumps
            sys.exit('Missing Boiler object')
        if self.source_is_exhaust_air():
            self.__volume_heated_all_services += volume_heated

        # Use low temperature test data for space heating - set flow temp such
        # that it matches the one used in the test
        temp_flow = self.__test_data._HeatPumpTestData__dsgn_flow_temps[0]

        # Design temperature difference across the emitters, in deg C or K
        temp_diff_emit_dsgn = max(
            temp_flow / 7.0,
            self.__test_data.temp_spread_test_conditions(temp_flow),
            )

        self.__create_service_connection(service_name)
        return HeatPumpServiceSpaceWarmAir(
            self,
            service_name,
            temp_diff_emit_dsgn,
            control,
            temp_flow,
            frac_convective,
            volume_heated,
            )

    def __get_temp_source(self):
        """ Get source temp according to rules in CALCM-01 - DAHPSE - V2.0_DRAFT13, 3.1.1 """
        if self.__source_type == SourceType.GROUND:
            # Subject to max source temp of 8 degC and min of 0 degC
            temp_ext = self.__external_conditions.air_temp()
            temp_source = max(0, min(8, temp_ext * 0.25806 + 2.8387))
        elif self.__source_type == SourceType.OUTSIDE_AIR:
            temp_source = self.__external_conditions.air_temp()
        elif self.__source_type == SourceType.EXHAUST_AIR_MEV:
            temp_source = self.__project.temp_internal_air()
        elif self.__source_type == SourceType.EXHAUST_AIR_MVHR:
            temp_source = self.__project.temp_internal_air()
        elif self.__source_type == SourceType.EXHAUST_AIR_MIXED:
            # Mixed exhaust air heat pumps use proportion of airflow from outside
            # combined with the internal air. However, when the external temperature
            # is above a maximum temperature threshold only the internal air is used
            # and when the mixed air temperature is below a minimum temperature threshold. 
            temp_ext = self.__external_conditions.air_temp()
            temp_int = self.__project.temp_internal_air()
            temp_mixed = self.__ext_air_ratio * temp_ext + (1 - self.__ext_air_ratio) * temp_int
            if  temp_ext > self.__eahp_mixed_max_temp \
                or temp_mixed < self.__eahp_mixed_min_temp:
                temp_source =  temp_int
            else:
                temp_source = temp_mixed

        elif self.__source_type == SourceType.WATER_GROUND:
            # Where water is extracted from the ground and re-injected into the
            # ground or discharged at the surface, the surface temperature is assumed
            # to be constant and equal to the annual average air temperature.
            temp_source = self.__external_conditions.air_temp_annual()
        elif self.__source_type == SourceType.WATER_SURFACE:
            # Where water is extracted from surface water, such as rivers and lakes,
            # it is assumed that this extraction does not substantively effect the
            # average temperature of the water volume, thus it must have a sufficient
            # thermal capacity. The source temperature is taken as the monthly air
            # temperature.
            temp_source = self.__external_conditions.air_temp_monthly()
        elif self.__source_type == SourceType.HEAT_NETWORK:
            temp_source = self.__temp_distribution_heat_network
        else:
            # If we reach here, then earlier input validation has failed, or a
            # SourceType option is missing above.
            sys.exit('SourceType not valid.')

        return Celcius2Kelvin(temp_source)

    def __thermal_capacity_op_cond(self, temp_output, temp_source):
        """ Calculate the thermal capacity of the heat pump at operating conditions

        Based on CALCM-01 - DAHPSE - V2.0_DRAFT13, section 4.4
        """
        if not self.__source_type == SourceType.OUTSIDE_AIR \
        and not self.__var_flow_temp_ctrl_during_test:
            thermal_capacity_op_cond = self.__test_data.average_capacity(temp_output)
        else:
            thermal_capacity_op_cond \
                = self.__test_data.capacity_op_cond_var_flow_or_source_temp(
                    temp_output,
                    temp_source,
                    self.__modulating_ctrl,
                    )

        return thermal_capacity_op_cond

    def __backup_energy_output_max(
            self,
            temp_output,
            temp_return_feed,
            time_available,
            hybrid_boiler_service=None,
            ):

        if self.__backup_ctrl == BackupCtrlType.TOPUP:
            # For top up mode, the time passed to the functions
            # should be set to None, as the time available is independant.
            time_elapsed_hp = None
        elif self.__backup_ctrl == BackupCtrlType.SUBSTITUTE \
            or self.__backup_ctrl == BackupCtrlType.NONE: \
            time_elapsed_hp = self.__total_time_running_current_timestep
        else:
            sys.exit('Invalid BackupCtrlType')

        if hybrid_boiler_service is not None:
            if isinstance(hybrid_boiler_service, BoilerServiceWaterRegular):
                energy_output_max = hybrid_boiler_service.energy_output_max(
                    Kelvin2Celcius(temp_return_feed),
                    time_elapsed_hp,
                    )
            elif isinstance(hybrid_boiler_service, BoilerServiceSpace):
                energy_output_max = hybrid_boiler_service.energy_output_max(
                    Kelvin2Celcius(temp_output),
                    Kelvin2Celcius(temp_return_feed),
                    time_elapsed_hp,
                    )

        else:
            energy_output_max = self.__power_max_backup * time_available

        return energy_output_max

    def __energy_output_max(
            self,
            temp_output,
            temp_return_feed,
            hybrid_boiler_service,
            service_type=None,
            temp_spread_correction=1.0,
            emitters_data_for_buffer_tank=None,
            service_name=None,
            ): 
        """ Calculate the maximum energy output of the HP, accounting for time
            spent on higher-priority services

        Note: Call via a HeatPumpService object, not directly.
        """
        timestep = self.__simulation_time.timestep()
        time_available = timestep - self.__total_time_running_current_timestep
        temp_source = self.__get_temp_source()

        # If there is buffer tank, energy output max will be affected by the temperature lift required
        # by the tank 
        heat_loss_buffer_kWh = 0.0
        if self.__buffer_tank is not None and emitters_data_for_buffer_tank is not None:
            buffer_tank_results = self.__buffer_tank.calc_buffer_tank(
                        service_name,
                        emitters_data_for_buffer_tank,
                        )

            if 'flow_temp_increase_due_to_buffer' in buffer_tank_results[-1]:
                flow_temp_increase_due_to_buffer = buffer_tank_results[-1]['flow_temp_increase_due_to_buffer']
                temp_output += flow_temp_increase_due_to_buffer
            if 'heat_loss_buffer_kWh' in buffer_tank_results:
                heat_loss_buffer_kWh = self.__buffer_tank.get_buffer_loss() + buffer_tank_results['heat_loss_buffer_kWh']
                
            emitters_data_for_buffer_tank['results'] = buffer_tank_results[-1]

        if self.__cost_schedule_hybrid_hp is not None:
            cop_op_cond, _ = self.__cop_deg_coeff_op_cond(
                service_type,
                temp_output,
                temp_source,
                temp_spread_correction,
                )
            energy_output_max_boiler = self.__backup_energy_output_max(
                temp_output,
                temp_return_feed,
                time_available,
                hybrid_boiler_service,
                )
            boiler_eff = self.__boiler.calc_boiler_eff(
                service_type,
                Kelvin2Celcius(temp_return_feed),
                energy_output_max_boiler,
                timestep
                )
            hp_cost_effective = self.__is_heat_pump_cost_effective(cop_op_cond, boiler_eff)
        else:
            hp_cost_effective = True

        if not hp_cost_effective:
            energy_max = self.__backup_energy_output_max(
                temp_output,
                temp_return_feed,
                time_available,
                hybrid_boiler_service,
                )
        else:
            if self.__outside_operating_limits(temp_return_feed):
                power_max_HP = 0.0
            else:
                power_max_HP = self.__thermal_capacity_op_cond(temp_output, temp_source)

            if self.__backup_ctrl == BackupCtrlType.NONE \
            or not self.__backup_heater_delay_time_elapsed():
                energy_max = power_max_HP * time_available
            elif self.__backup_ctrl == BackupCtrlType.TOPUP:
                energy_max_backup = self.__backup_energy_output_max(
                    temp_output,
                    temp_return_feed,
                    time_available,
                    hybrid_boiler_service,
                    )
                energy_max = (power_max_HP * time_available) + energy_max_backup
            elif self.__backup_ctrl == BackupCtrlType.SUBSTITUTE:
                energy_max_backup = self.__backup_energy_output_max(
                    temp_output,
                    temp_return_feed,
                    time_available,
                    hybrid_boiler_service,
                    )
                energy_max = max(power_max_HP * time_available, energy_max_backup)

        if self.__buffer_tank is not None and emitters_data_for_buffer_tank is not None:
            return energy_max - heat_loss_buffer_kWh, emitters_data_for_buffer_tank
        else:
            return energy_max - heat_loss_buffer_kWh

    def __cop_deg_coeff_op_cond(
            self,
            service_type,
            temp_output, # Kelvin
            temp_source, # Kelvin
            temp_spread_correction,
            ):
        """ Calculate CoP and degradation coefficient at operating conditions """
        if callable(temp_spread_correction):
            temp_spread_correction_factor = temp_spread_correction(temp_output, temp_source)
        else:
            temp_spread_correction_factor = temp_spread_correction

        # TODO Make if/elif/else chain exhaustive?
        if not self.__source_type == SourceType.OUTSIDE_AIR \
        and not self.__var_flow_temp_ctrl_during_test:
            cop_op_cond \
                = temp_spread_correction_factor \
                * self.__test_data.cop_op_cond_if_not_air_source(
                    self.__temp_diff_limit_low,
                    self.__external_conditions.air_temp(),
                    temp_source,
                    temp_output,
                    )
            deg_coeff_op_cond = self.__test_data.average_degradation_coeff(temp_output)
        else:
            carnot_cop_op_cond = carnot_cop(temp_source, temp_output, self.__temp_diff_limit_low)
            # Get exergy load ratio at operating conditions and exergy load ratio,
            # exergy efficiency and degradation coeff at test conditions above and
            # below operating conditions
            lr_op_cond = self.__test_data.lr_op_cond(temp_output, temp_source, carnot_cop_op_cond)
            lr_below, lr_above, eff_below, eff_above, deg_coeff_below, deg_coeff_above \
                = self.__test_data.lr_eff_degcoeff_either_side_of_op_cond(temp_output, lr_op_cond)

            # CALCM-01 - DAHPSE - V2.0_DRAFT13, section 4.5.4
            # Get exergy efficiency by interpolating between figures above and
            # below operating conditions
            exer_eff_op_cond \
                = eff_below \
                + (eff_below - eff_above) \
                * (lr_op_cond - lr_below) \
                / (lr_below - lr_above)

            # CALCM-01 - DAHPSE - V2.0_DRAFT13, section 4.5.5
            # Note: DAHPSE method document section 4.5.5 doesn't have
            # temp_spread_correction_factor in formula below. However, section 4.5.7
            # states that the correction factor is to be applied to the CoP.
            cop_op_cond = max(
                1.0,
                exer_eff_op_cond * carnot_cop_op_cond * temp_spread_correction_factor,
                )

            if self.__sink_type == SinkType.AIR and service_type != ServiceType.WATER:
                limit_upper = 0.25
            else:
                limit_upper = 1.0

            if self.__sink_type == SinkType.AIR and service_type != ServiceType.WATER:
                limit_lower = 0.0
            else:
                limit_lower = 0.9

            if lr_below == lr_above:
                deg_coeff_op_cond = deg_coeff_below
            else:
                deg_coeff_op_cond \
                    = deg_coeff_below \
                    + (deg_coeff_below - deg_coeff_above) \
                    * (lr_op_cond - lr_below) \
                    / (lr_below - lr_above)

            deg_coeff_op_cond = max(min(deg_coeff_op_cond, limit_upper), limit_lower)

        return cop_op_cond, deg_coeff_op_cond

    def __energy_output_limited(
            self,
            energy_output_required,
            temp_output,
            temp_used_for_scaling,
            temp_limit_upper
            ):
        """ Calculate energy output limited by upper temperature """
        if temp_output > temp_limit_upper:
        # If required output temp is above upper limit
            if temp_output == temp_used_for_scaling:
            # If flow and return temps are equal
                return energy_output_required
            else:
            # If flow and return temps are not equal
                if (temp_limit_upper - temp_used_for_scaling) >= self.__temp_diff_flow_return_min:
                # If max. achievable temp diff is at least the min required
                # for the HP to operate.
                    return \
                          energy_output_required \
                        * (temp_limit_upper - temp_used_for_scaling) \
                        / (temp_output - temp_used_for_scaling)
                else:
                # If max. achievable temp diff is less than the min required
                # for the HP to operate.
                    return 0.0
        else:
        # If required output temp is below upper limit
            return energy_output_required

    def __backup_heater_delay_time_elapsed(self):
        """ Check if backup heater is available or still in delay period """
        return self.__time_running_continuous >= self.__time_delay_backup

    def __outside_operating_limits(self, temp_return_feed):
        """ Check if heat pump is outside operating limits """
        temp_source = self.__get_temp_source()
        below_min_ext_temp = temp_source <= self.__temp_lower_op_limit

        if self.__sink_type == SinkType.WATER:
            above_temp_return_feed_max = temp_return_feed > self.__temp_return_feed_max
        elif self.__sink_type == SinkType.AIR:
            above_temp_return_feed_max = False
        else:
            sys.exit('Return feed temp check not defined for sink type')

        return below_min_ext_temp or above_temp_return_feed_max

    def __inadequate_capacity(
            self,
            energy_output_required,
            thermal_capacity_op_cond,
            temp_output,
            time_available,
            temp_return_feed,
            hybrid_boiler_service,
            ):
        """ Check if heat pump has adequate capacity to meet demand """
        timestep = self.__simulation_time.timestep()

        # For top-up backup heater, use backup if delay time has elapsed.
        # For substitute backup heater, use backup if delay time has elapsed and
        # backup heater can provide more energy than heat pump. This assumption
        # is required to make the maximum energy output of the system
        # predictable before the demand is known.
        energy_max_backup = self.__backup_energy_output_max(
                temp_output,
                temp_return_feed,
                time_available,
                hybrid_boiler_service,
                )

        if (   self.__backup_ctrl == BackupCtrlType.TOPUP 
           and self.__backup_heater_delay_time_elapsed()
           ) \
        or (   self.__backup_ctrl == BackupCtrlType.SUBSTITUTE
           and self.__backup_heater_delay_time_elapsed()
           and energy_max_backup > thermal_capacity_op_cond * time_available
           ):
            inadequate_capacity = energy_output_required > thermal_capacity_op_cond * timestep
        else:
            inadequate_capacity = False

        return inadequate_capacity

    def __is_heat_pump_cost_effective(
            self,
            cop_op_cond,
            boiler_eff,
            ):

        cost_hp = self.__cost_schedule_hp[self.__simulation_time.time_series_idx(self.__cost_schedule_start_day, self.__cost_schedule_time_series_step)]
        cost_boiler = self.__cost_schedule_boiler[self.__simulation_time.time_series_idx(self.__cost_schedule_start_day, self.__cost_schedule_time_series_step)]

        if cost_hp / cop_op_cond  <= cost_boiler / boiler_eff:
            hp_cost_effective  = True
        else:
            hp_cost_effective = False
        return hp_cost_effective

    def __use_backup_heater_only(
            self,
            cop_op_cond,
            energy_output_required,
            thermal_capacity_op_cond,
            temp_output,
            time_available,
            temp_return_feed,
            hybrid_boiler_service=None,
            boiler_eff=None,
            ):
        """ Evaluate boolean conditions that may trigger backup heater """
        outside_operating_limits = self.__outside_operating_limits(temp_return_feed)
        inadequate_capacity \
            = self.__inadequate_capacity(
                energy_output_required,
                thermal_capacity_op_cond,
                temp_output,
                time_available,
                temp_return_feed,
                hybrid_boiler_service,
                )

        hp_not_cost_effective = False
        if self.__cost_schedule_hybrid_hp is not None:
            hp_not_cost_effective = not self.__is_heat_pump_cost_effective(
                cop_op_cond,
                boiler_eff,
                )

        return  self.__backup_ctrl != BackupCtrlType.NONE \
            and (  outside_operating_limits \
                or (inadequate_capacity and self.__backup_ctrl == BackupCtrlType.SUBSTITUTE) \
                or hp_not_cost_effective
                )

    def __run_demand_energy_calc(
            self,
            service_name,
            service_type,
            energy_output_required,
            temp_output, # Kelvin
            temp_return_feed, # Kelvin
            temp_limit_upper, # Kelvin
            time_constant_for_service,
            service_on, # bool - is service allowed to run?
            temp_spread_correction=1.0,
            temp_used_for_scaling=None,
            hybrid_boiler_service=None,
            boiler_eff = None,
            additional_time_unavailable=0.0,
            emitters_data_for_buffer_tank=None
            ):
        """ Calculate energy required by heat pump to satisfy demand for the service indicated.

        Note: Call via the __demand_energy func, not directly.
              This function should not save any results to member variables of
              this class, because it may need to be run more than once (e.g. for
              exhaust air heat pumps). Results should be returned to the
              __demand_energy function which calls this one and will save results
              when appropriate.
        Note: The optional variable additional_time_unavailable is used for
              calculating running time without needing to update any state - the
              variable contains the time already committed to other services
              where the running time has not been added to
              self.__total_time_running_current_timestep
        """
        flow_temp_increase_due_to_buffer = 0.0    
        power_buffer_tank_pump = 0.0
        heat_loss_buffer_kWh = 0.0
        if self.__buffer_tank is not None and emitters_data_for_buffer_tank is not None:
            buffer_tank_results = emitters_data_for_buffer_tank['results']
            if 'flow_temp_increase_due_to_buffer' in buffer_tank_results:
                flow_temp_increase_due_to_buffer = buffer_tank_results['flow_temp_increase_due_to_buffer']
            if 'pump_power_at_flow_rate' in buffer_tank_results:
                power_buffer_tank_pump = buffer_tank_results['pump_power_at_flow_rate']
            if 'heat_loss_buffer_kWh' in buffer_tank_results:
                heat_loss_buffer_kWh = self.__buffer_tank.get_buffer_loss() + buffer_tank_results['heat_loss_buffer_kWh']
            # If the service is not working, all the buffer loss is accumulated for the next time step
            if not service_on:
                self.__buffer_tank.update_buffer_loss(heat_loss_buffer_kWh)
                heat_loss_buffer_kWh = 0.0

        # Adding buffer tank losses to the energy required to be delivered by the HP
        energy_output_required += heat_loss_buffer_kWh

        temp_output += flow_temp_increase_due_to_buffer
        
        if temp_used_for_scaling is None:
            temp_used_for_scaling = temp_return_feed

        timestep = self.__simulation_time.timestep()

        energy_output_limited = self.__energy_output_limited(
            energy_output_required,
            temp_output,
            temp_used_for_scaling,
            temp_limit_upper,
            )

        temp_source = self.__get_temp_source() # Kelvin
        # From here onwards, output temp to be used is subject to the upper limit
        temp_output = min(temp_output, temp_limit_upper) # Kelvin

        # Get thermal capacity, CoP and degradation coeff at operating conditions
        thermal_capacity_op_cond = self.__thermal_capacity_op_cond(temp_output, temp_source)
        cop_op_cond, deg_coeff_op_cond = self.__cop_deg_coeff_op_cond(
            service_type,
            temp_output,
            temp_source,
            temp_spread_correction,
            )

        # Calculate running time of HP
        time_required = energy_output_limited / thermal_capacity_op_cond
        time_available = timestep - self.__total_time_running_current_timestep - additional_time_unavailable
        time_running_current_service = min(time_required, time_available)

        # Calculate load ratio
        load_ratio = time_running_current_service / timestep
        if self.__modulating_ctrl:
            if 55.0 in self.__test_data._HeatPumpTestData__dsgn_flow_temps:
                load_ratio_continuous_min \
                    = ( min(
                            max(temp_output, self.__temp_min_modulation_rate_low),
                            self.__temp_min_modulation_rate_high
                            )
                      - self.__temp_min_modulation_rate_low
                      ) \
                    / (self.__temp_min_modulation_rate_high - self.__temp_min_modulation_rate_low) \
                    * self.__min_modulation_rate_low \
                    + ( 1.0 
                      - ( min(
                              max(temp_output, self.__temp_min_modulation_rate_low),
                              self.__temp_min_modulation_rate_high
                          )
                        - self.__temp_min_modulation_rate_low
                        )
                      ) \
                    / (self.__temp_min_modulation_rate_high - self.__temp_min_modulation_rate_low) \
                    * self.__min_modulation_rate_55
            else:
                load_ratio_continuous_min = self.__min_modulation_rate_low
        else:
            # On/off heat pump cannot modulate below maximum power
            load_ratio_continuous_min = 1.0

        # Determine whether or not HP is operating in on/off mode
        hp_operating_in_onoff_mode = (load_ratio > 0.0 and load_ratio < load_ratio_continuous_min)

        compressor_power_full_load = thermal_capacity_op_cond / cop_op_cond

        # CALCM-01 - DAHPSE - V2.0_DRAFT13, section 4.5.10, step 1:
        compressor_power_min_load \
            = compressor_power_full_load * load_ratio_continuous_min
        # CALCM-01 - DAHPSE - V2.0_DRAFT13, section 4.5.10, step 2:
        # TODO The value calculated here is only used in some code branches
        #      later. In future, rearrange this function so that this is only
        #      calculated when needed.
        if load_ratio >= load_ratio_continuous_min:
            power_used_due_to_inertia_effects = 0.0
        else:
            power_used_due_to_inertia_effects \
                = compressor_power_min_load \
                * self.__time_constant_onoff_operation \
                * load_ratio \
                * (1.0 - load_ratio) \
                / time_constant_for_service

        # TODO Consider moving some of these checks earlier or to HeatPumpService
        #      classes. May be able to skip a lot of the calculation.
        if hybrid_boiler_service is not None:
            boiler_eff = self.__boiler.calc_boiler_eff(
                service_type,
                Kelvin2Celcius(temp_return_feed),
                energy_output_required,
                timestep
                )

        use_backup_heater_only = self.__use_backup_heater_only(
            cop_op_cond,
            energy_output_required,
            thermal_capacity_op_cond,
            temp_output,
            time_available,
            temp_return_feed,
            hybrid_boiler_service,
            boiler_eff,
            )

        # Calculate energy delivered by HP and energy input
        if use_backup_heater_only or not service_on:
            energy_delivered_HP = 0.0
            energy_input_HP = 0.0
            energy_input_HP_divisor = None
            # If using back heater only then heat pump does not run
            # therefore set time running to zero.
            time_running_current_service = 0.0
        else:
            # Backup heater not providing entire energy requirement
            energy_delivered_HP = thermal_capacity_op_cond * time_running_current_service

            if hp_operating_in_onoff_mode:
                # TODO Why does the divisor below differ for DHW from warm air HPs?
                if service_type == ServiceType.WATER and self.__sink_type == SinkType.AIR:
                    energy_input_HP_divisor \
                        = 1.0 \
                        - deg_coeff_op_cond \
                        * (1.0 - load_ratio / load_ratio_continuous_min)
                else:
                    energy_input_HP_divisor = 1.0

                # TODO In the eqn below, should compressor_power_full_load actually
                #      be compressor_power_min_load? In this code branch the HP is
                #      operating in on/off mode, so presumably cycling between off
                #      and min load rather than full load.
                # Note: Energy_ancillary_when_off should also be included in the
                # energy input for on/off operation, but at this stage we have
                # not calculated whether a lower-priority service will run
                # instead, so this will need to be calculated later and
                # (energy_ancillary_when_off / eqn_denom) added to the energy
                # input
                energy_input_HP \
                    = ( ( compressor_power_full_load * (1.0 + self.__f_aux) \
                        + power_used_due_to_inertia_effects \
                        ) \
                      * time_running_current_service \
                      ) \
                    / energy_input_HP_divisor
            else:
                # If not operating in on/off mode
                energy_input_HP = energy_delivered_HP / cop_op_cond
                energy_input_HP_divisor = None

        # Calculate energy delivered by backup heater
        if self.__backup_ctrl == BackupCtrlType.NONE \
        or not self.__backup_heater_delay_time_elapsed() and not use_backup_heater_only \
        or not service_on:
            energy_delivered_backup = 0.0
        elif self.__backup_ctrl == BackupCtrlType.TOPUP \
        or self.__backup_ctrl == BackupCtrlType.SUBSTITUTE:
            energy_max_backup = self.__backup_energy_output_max(
                        temp_output,
                        temp_return_feed,
                        time_available,
                        hybrid_boiler_service,
                        )
            energy_delivered_backup = max(
                min(
                    energy_max_backup,
                    energy_output_required - energy_delivered_HP,
                ),
                0.0,
                )
        else:
            sys.exit('Invalid BackupCtrlType')

        # Calculate energy input to backup heater
        # TODO Account for backup heater efficiency, or call another heating
        #      system object. For now, assume 100% efficiency
        energy_input_backup = energy_delivered_backup

        if hybrid_boiler_service is not None:
            energy_output_required_boiler = energy_delivered_backup
            energy_delivered_backup = 0.0
            energy_input_backup = 0.0
        else:
            energy_output_required_boiler = 0.0

        # Energy used by pumps
        energy_source_circ_pump \
            = time_running_current_service * self.__power_source_circ_pump
        if service_type == ServiceType.SPACE and self.__sink_type == SinkType.AIR:
            # If warm air distribution add electricity for warm air fan
            energy_heating_warm_air_fan \
                = time_running_current_service * self.__power_heating_warm_air_fan
            energy_heating_circ_pump = 0
        else:
            # if wet distribution add electricity for wet distribution. 
            energy_heating_warm_air_fan = 0
            energy_heating_circ_pump \
                = time_running_current_service * (self.__power_heating_circ_pump + power_buffer_tank_pump)
           
        # Calculate total energy delivered and input
        energy_delivered_total = energy_delivered_HP + energy_delivered_backup
        energy_input_total \
            = energy_input_HP + energy_input_backup \
            + energy_heating_circ_pump + energy_source_circ_pump + energy_heating_warm_air_fan
            
        if self.__buffer_tank is not None:
            if energy_delivered_total >= heat_loss_buffer_kWh:
                energy_delivered_total -= heat_loss_buffer_kWh
                self.__buffer_tank.update_buffer_loss(0.0)
            else:
                hp_covered_buffer_losses = energy_delivered_total
                energy_delivered_total = 0.0
                self.__buffer_tank.update_buffer_loss(heat_loss_buffer_kWh - hp_covered_buffer_losses)

        return {
            'service_name': service_name,
            'service_type': service_type,
            'service_on': service_on,
            'energy_output_required': energy_output_required,
            'temp_output': temp_output,
            'temp_source': temp_source,
            'cop_op_cond': cop_op_cond,
            'thermal_capacity_op_cond': thermal_capacity_op_cond,
            'time_running': time_running_current_service,
            'deg_coeff_op_cond': deg_coeff_op_cond,
            'compressor_power_min_load': compressor_power_min_load,
            'load_ratio_continuous_min': load_ratio_continuous_min,
            'load_ratio': load_ratio,
            'use_backup_heater_only': use_backup_heater_only,
            'hp_operating_in_onoff_mode': hp_operating_in_onoff_mode,
            'energy_input_HP_divisor': energy_input_HP_divisor,
            'energy_input_HP': energy_input_HP,
            'energy_delivered_HP': energy_delivered_HP,
            'energy_input_backup': energy_input_backup,
            'energy_delivered_backup': energy_delivered_backup,
            'energy_input_total': energy_input_total,
            'energy_delivered_total': energy_delivered_total,
            'energy_heating_circ_pump': energy_heating_circ_pump,
            'energy_source_circ_pump': energy_source_circ_pump,
            'energy_output_required_boiler':energy_output_required_boiler,
            'energy_heating_warm_air_fan' : energy_heating_warm_air_fan,
            }

    def __demand_energy(
            self,
            service_name,
            service_type,
            energy_output_required,
            temp_output, # Kelvin
            temp_return_feed, # Kelvin
            temp_limit_upper, # Kelvin
            time_constant_for_service,
            service_on, # bool - is service allowed to run?
            temp_spread_correction=1.0,
            temp_used_for_scaling=None,
            hybrid_boiler_service = None,
            emitters_data_for_buffer_tank=None
            ):
        """ Calculate energy required by heat pump to satisfy demand for the service indicated.

        Note: Call via a HeatPumpService object, not directly.
        """
        service_results = self.__run_demand_energy_calc(
                service_name,
                service_type,
                energy_output_required,
                temp_output,
                temp_return_feed,
                temp_limit_upper,
                time_constant_for_service,
                service_on,
                temp_spread_correction = temp_spread_correction,
                temp_used_for_scaling = temp_used_for_scaling,
                hybrid_boiler_service = hybrid_boiler_service,
                emitters_data_for_buffer_tank = emitters_data_for_buffer_tank
                )


        if hybrid_boiler_service is not None:
            # Call demand function for boiler and 
            # return the boiler running time so that it can be added to the total running time
            hybrid_service_bool = True
            
            if self.__backup_ctrl == BackupCtrlType.TOPUP:
                # For top up mode, the time passed to the functions
                # should be set to None, as the time available is independant.
                time_elapsed_hp = None
            elif self.__backup_ctrl == BackupCtrlType.SUBSTITUTE \
                or self.__backup_ctrl == BackupCtrlType.NONE: \
                time_elapsed_hp = self.__total_time_running_current_timestep
            else:
                sys.exit('Invalid BackupCtrlType')

            if isinstance(hybrid_boiler_service, BoilerServiceWaterRegular):
                service_results['energy_output_delivered_boiler'], time_running_boiler = \
                    hybrid_boiler_service.demand_energy(
                        service_results['energy_output_required_boiler'],
                        Kelvin2Celcius(temp_return_feed),
                        hybrid_service_bool,
                        time_elapsed_hp,
                    )
            elif isinstance(hybrid_boiler_service, BoilerServiceSpace):
                service_results['energy_output_delivered_boiler'], time_running_boiler = \
                    hybrid_boiler_service.demand_energy(
                        service_results['energy_output_required_boiler'],
                        Kelvin2Celcius(temp_output),
                        Kelvin2Celcius(temp_return_feed),
                        hybrid_service_bool,
                        time_elapsed_hp,
                        )
            if self.__backup_ctrl == BackupCtrlType.SUBSTITUTE:
                self.__total_time_running_current_timestep += time_running_boiler
        else:
            service_results['energy_output_delivered_boiler'] = 0.0

        # Save results that are needed later (in the timestep_end function)
        self.__service_results.append(service_results)
        self.__total_time_running_current_timestep \
            += service_results['time_running']

        # Feed/return results to other modules
        self.__energy_supply_connections[service_name].demand_energy(
            service_results['energy_input_total']
            )
        return service_results['energy_delivered_total'] + service_results['energy_output_delivered_boiler']

    def __running_time_throughput_factor(
            self,
            space_heat_running_time_cumulative,
            service_name,
            service_type,
            energy_output_required,
            temp_output,
            temp_return_feed,
            temp_limit_upper,
            time_constant_for_service,
            service_on,
            volume_heated_by_service,
            temp_spread_correction,
            ):
        """ Return the cumulative running time and throughput factor (exhaust air HPs only) """

        # TODO Run HP calculation to get total running time incl space heating,
        #      but do not save space heating running time
        service_results = self.__run_demand_energy_calc(
                service_name,
                service_type,
                energy_output_required,
                temp_output,
                temp_return_feed,
                temp_limit_upper,
                time_constant_for_service,
                service_on,
                temp_spread_correction,
                additional_time_unavailable=space_heat_running_time_cumulative,
                )

        throughput_factor = self.__calc_throughput_factor(service_results['time_running'])

        # Adjust throughput factor to match simplifying assumption that all
        # extra ventilation required due to space heating demand in the current
        # zone is assigned to current zone
        throughput_factor_zone \
            = (throughput_factor - 1.0) \
            * self.__volume_heated_all_services / volume_heated_by_service \
            + 1.0

        return service_results['time_running'], throughput_factor_zone

    def __calc_throughput_factor(self, time_running):
        timestep = self.__simulation_time.timestep()

        # Apply overventilation ratio to part of timestep where HP is running
        # to calculate throughput_factor.
        throughput_factor \
            = ( (timestep - time_running) \
              + self.__overvent_ratio * time_running \
              ) \
            / timestep
        return throughput_factor

    def throughput_factor(self):
        return self.__calc_throughput_factor(self.__total_time_running_current_timestep)

    def __calc_ancillary_energy(self, timestep, time_remaining_current_timestep):
        """ Calculate ancillary energy for each service """
        for service_no, service_data in enumerate(self.__service_results):
            # Unpack results of previous calculations for this service
            service_name = service_data['service_name']
            service_type = service_data['service_type']
            service_on = service_data['service_on']
            time_running_current_service = service_data['time_running']
            deg_coeff_op_cond = service_data['deg_coeff_op_cond']
            compressor_power_min_load = service_data['compressor_power_min_load']
            load_ratio_continuous_min = service_data['load_ratio_continuous_min']
            load_ratio = service_data['load_ratio']
            use_backup_heater_only = service_data['use_backup_heater_only']
            hp_operating_in_onoff_mode = service_data['hp_operating_in_onoff_mode']
            energy_input_HP_divisor = service_data['energy_input_HP_divisor']

            time_running_subsequent_services \
                = sum([ \
                    data['time_running'] \
                    for data in self.__service_results[service_no + 1 :] \
                    ])

            if service_on \
            and time_running_current_service > 0.0 and not time_running_subsequent_services > 0.0 \
            and not (self.__sink_type == SinkType.AIR and service_type == ServiceType.WATER):
                energy_ancillary_when_off \
                    = (1.0 - deg_coeff_op_cond) \
                    * (compressor_power_min_load / load_ratio_continuous_min) \
                    * max(
                        ( time_remaining_current_timestep \
                        - load_ratio / load_ratio_continuous_min * timestep
                        ),
                        0.0
                        )
            else:
                energy_ancillary_when_off = 0.0

            if not use_backup_heater_only and hp_operating_in_onoff_mode:
                energy_input_HP = energy_ancillary_when_off / energy_input_HP_divisor
            else:
                energy_input_HP = 0.0

            self.__energy_supply_connections[service_name].demand_energy(energy_input_HP)
            self.__service_results[service_no]['energy_input_HP'] += energy_input_HP
            self.__service_results[service_no]['energy_input_total'] += energy_input_HP

    def __calc_auxiliary_energy(self, timestep, time_remaining_current_timestep):
        """ Calculate auxiliary energy according to CALCM-01 - DAHPSE - V2.0_DRAFT13, section 4.7 """

        # Retrieve control settings for this timestep
        heating_profile_on = False
        water_profile_on = False
        for service_data in self.__service_results:
            if service_data['service_type'] == ServiceType.SPACE:
                heating_profile_on = service_data['service_on']
            elif service_data['service_type'] == ServiceType.WATER:
                water_profile_on = service_data['service_on']
            else:
                sys.exit('ServiceType not recognised')

        # Energy used in standby and crankcase heater mode
        # TODO Crankcase heater mode appears to be relevant only when HP is
        #      available to provide space heating. Therefore, it could be added
        #      to space heating energy consumption instead of auxiliary
        # TODO Standby power is only relevant when at least one service is
        #      available. Therefore, it could be split between the available
        #      services rather than treated as auxiliary
        energy_off_mode = 0.0
        energy_standby = 0.0
        energy_crankcase_heater_mode = 0.0
        if heating_profile_on:
            energy_standby = time_remaining_current_timestep * self.__power_standby
            energy_crankcase_heater_mode \
                = time_remaining_current_timestep * self.__power_crankcase_heater_mode
        elif not heating_profile_on and water_profile_on:
            energy_standby = time_remaining_current_timestep * self.__power_standby
        # Energy used in off mode
        elif not heating_profile_on and not water_profile_on:
            energy_off_mode = timestep * self.__power_off_mode
        else:
            sys.exit() # Should never get here.

        energy_aux = energy_standby + energy_crankcase_heater_mode + energy_off_mode
        self.__energy_supply_connection_aux.demand_energy(energy_aux)
        return energy_standby, energy_crankcase_heater_mode, energy_off_mode

    def __extract_energy_from_source(self):
        """ If HP uses heat network as source, calculate energy extracted from heat network """
        for service_data in self.__service_results:
            service_name = service_data['service_name']
            energy_delivered_HP = service_data['energy_delivered_HP']
            energy_input_HP = service_data['energy_input_HP']
            energy_extracted_HP = energy_delivered_HP - energy_input_HP
            self.__energy_supply_HN_connections[service_name].demand_energy(energy_extracted_HP)

    def timestep_end(self):
        """ Calculations to be done at the end of each timestep """
        timestep = self.__simulation_time.timestep()
        time_remaining_current_timestep = timestep - self.__total_time_running_current_timestep

        if time_remaining_current_timestep == 0.0:
            self.__time_running_continuous += self.__total_time_running_current_timestep
        else:
            self.__time_running_continuous = 0.0

        self.__calc_ancillary_energy(timestep, time_remaining_current_timestep)
        energy_standby, energy_crankcase_heater_mode, energy_off_mode \
               = self.__calc_auxiliary_energy(timestep, time_remaining_current_timestep)

        if self.__source_type == SourceType.HEAT_NETWORK:
            self.__extract_energy_from_source()

        # If detailed results are to be output, save the results from the current timestep
        if self.__detailed_results is not None:
            self.__service_results.append({
                'energy_standby': energy_standby,
                'energy_crankcase_heater_mode': energy_crankcase_heater_mode,
                'energy_off_mode': energy_off_mode,
                })
            self.__detailed_results.append(self.__service_results)

        # Variables below need to be reset at the end of each timestep.
        self.__total_time_running_current_timestep = 0.0
        self.__service_results = []

    def output_detailed_results(self, hot_water_energy_output):
        """ Output detailed results of heat pump calculation """

        # Define parameters to output
        # Second element of each tuple controls whether item is summed for annual total
        output_parameters = [
            ('service_name', None, False),
            ('service_type', None, False),
            ('service_on', None, False),
            ('energy_output_required', 'kWh', True),
            ('temp_output', 'K', False),
            ('temp_source', 'K', False),
            ('thermal_capacity_op_cond', 'kW', False),
            ('cop_op_cond', None, False),
            ('time_running', 'hours', True),
            ('load_ratio', None, False),
            ('hp_operating_in_onoff_mode', None, False),
            ('energy_delivered_HP', 'kWh', True),
            ('energy_delivered_backup', 'kWh', True),
            ('energy_delivered_total', 'kWh', True),
            ('energy_input_HP', 'kWh', True),
            ('energy_input_backup', 'kWh', True),
            ('energy_heating_circ_pump', 'kWh', True),
            ('energy_source_circ_pump', 'kWh', True),
            ('energy_heating_warm_air_fan', 'kWh', True),
            ('energy_input_total', 'kWh', True),
            ('energy_output_delivered_boiler','kWh', True)
            ]
        aux_parameters = [
            ('energy_standby', 'kWh', True),
            ('energy_crankcase_heater_mode', 'kWh', True),
            ('energy_off_mode', 'kWh', True),
            ]

        results_per_timestep = {'auxiliary': {}}
        # Report auxiliary parameters (not specific to a service)
        for parameter, param_unit, _ in aux_parameters:
            results_per_timestep['auxiliary'][(parameter, param_unit)] = []
            for t_idx, service_results in enumerate(self.__detailed_results):
                result = service_results[-1][parameter]
                results_per_timestep['auxiliary'][(parameter, param_unit)].append(result)
        # For each service, report required output parameters
        for service_idx, service_name in enumerate(self.__energy_supply_connections.keys()):
            results_per_timestep[service_name] = {}
            # Look up each required parameter
            for parameter, param_unit, _ in output_parameters:
                results_per_timestep[service_name][(parameter, param_unit)] = []
                # Look up value of required parameter in each timestep
                for t_idx, service_results in enumerate(self.__detailed_results):
                    result = service_results[service_idx][parameter]
                    results_per_timestep[service_name][(parameter, param_unit)].append(result)
            # For water heating service, record hot water energy delivered from tank
            if self.__detailed_results[0][service_idx]['service_type'] == ServiceType.WATER :
                # For DHW, need to include storage and primary circuit losses.
                # Can do this by replacing H4 numerator with total energy
                # draw-off from hot water cylinder.
                # TODO Note that the below assumes that there is only one water
                #      heating service and therefore that all hot water energy
                #      output is assigned to that service. If the model changes in
                #      future to allow more than one hot water system, this code may
                #      need to be revised to handle that scenario.
                results_per_timestep[service_name][('energy_delivered_H4', 'kWh')] \
                    = hot_water_energy_output
            else:
                # TODO Note that the below assumes there is no buffer tank for
                #      space heating, which is not currently included in the
                #      model. If this is included in future, this code will need
                #      to be revised.
                results_per_timestep[service_name][('energy_delivered_H4', 'kWh')] \
                    = results_per_timestep[service_name][('energy_delivered_total', 'kWh')]

        results_annual = {
            'Overall': {
                (parameter, param_units): 0.0
                for parameter, param_units, incl_in_annual in output_parameters
                if incl_in_annual
                },
            'auxiliary': {},
            }
        results_annual['Overall'][('energy_delivered_H4', 'kWh')] = 0.0
        # Report auxiliary parameters (not specific to a service)
        for parameter, param_unit, incl_in_annual in aux_parameters:
            if incl_in_annual:
                results_annual['auxiliary'][(parameter, param_unit)] \
                    = sum(results_per_timestep['auxiliary'][(parameter, param_unit)])
        # For each service, report required output parameters
        for service_idx, service_name in enumerate(self.__energy_supply_connections.keys()):
            results_annual[service_name] = {}
            for parameter, param_unit, incl_in_annual in output_parameters:
                if incl_in_annual:
                    parameter_annual_total \
                        = sum(results_per_timestep[service_name][(parameter, param_unit)])
                    results_annual[service_name][(parameter, param_unit)] = parameter_annual_total
                    results_annual['Overall'][(parameter, param_unit)] += parameter_annual_total
            results_annual[service_name][('energy_delivered_H4', 'kWh')] \
                = sum(results_per_timestep[service_name][('energy_delivered_H4', 'kWh')])
            results_annual['Overall'][('energy_delivered_H4', 'kWh')] \
                += results_annual[service_name][('energy_delivered_H4', 'kWh')]
            # For each service, calculate CoP at different system boundaries
            self.__calc_service_cop(results_annual[service_name])

        # Calculate overall CoP for all services combined
        self.__calc_service_cop(results_annual['Overall'], results_annual['auxiliary'])

        return results_per_timestep, results_annual

    def __calc_service_cop(self, results_totals, results_auxiliary=None):
        """ Calculate CoP for whole simulation period for the given service (or overall) """
        # Add auxiliary energy to overall CoP
        if results_auxiliary is not None:
            energy_auxiliary = sum(result for result in results_auxiliary.values())
        else:
            energy_auxiliary = 0.0

        # Calculate CoP at different system boundaries
        cop_h1_numerator = results_totals[('energy_delivered_HP', 'kWh')]
        cop_h1_denominator = results_totals[('energy_input_HP', 'kWh')] + energy_auxiliary
        cop_h2_numerator = cop_h1_numerator
        cop_h2_denominator \
            = cop_h1_denominator + results_totals[('energy_source_circ_pump', 'kWh')]
        cop_h3_numerator \
            = cop_h2_numerator + results_totals[('energy_delivered_backup', 'kWh')]
        cop_h3_denominator \
            = cop_h2_denominator + results_totals[('energy_input_backup', 'kWh')]
        cop_h4_numerator = results_totals[('energy_delivered_H4', 'kWh')]
        cop_h4_denominator \
            = cop_h3_denominator \
            + results_totals[('energy_heating_circ_pump', 'kWh')] \
            + results_totals[('energy_heating_warm_air_fan', 'kWh')]

        cop_h4_note = 'Note: For water heating services, only valid when HP is only heat source'

        if cop_h1_denominator == 0.0:
            results_totals[('CoP (H1)', None)] = 0.0
        else:
            results_totals[('CoP (H1)', None)] = cop_h1_numerator / cop_h1_denominator

        if cop_h2_denominator == 0.0:
            results_totals[('CoP (H2)', None)] = 0.0
        else:
            results_totals[('CoP (H2)', None)] = cop_h2_numerator / cop_h2_denominator

        if cop_h3_denominator == 0.0:
            results_totals[('CoP (H3)', None)] = 0.0
        else:
            results_totals[('CoP (H3)', None)] = cop_h3_numerator / cop_h3_denominator

        if cop_h4_denominator == 0.0:
            results_totals[('CoP (H4)', None)] = 0.0
        else:
            results_totals[('CoP (H4)', cop_h4_note)] = cop_h4_numerator / cop_h4_denominator

        return results_totals


class HeatPump_HWOnly:
    """ An object to represent an electric hot-water-only heat pump, tested to EN 16147 """

    def __init__(
            self,
            power_max,
            test_data,
            vol_daily_average,
            tank_volume,
            daily_losses,
            heat_exchanger_surface_area,
            in_use_factor_mismatch,
            tank_volume_declared,
            heat_exchanger_surface_area_declared,
            daily_losses_declared,
            energy_supply_conn,
            simulation_time,
            control=None,
            ):
        """ Construct a HeatPump_HWOnly object

        Arguments:
        power_max -- in kW
        test_data -- dictionary with keys denoting tapping profile letter (M or L)
                     and values being another dictionary containing the following:
                     - cop_dhw -- CoP measured during EN 16147 test
                     - hw_tapping_prof_daily_total -- daily energy requirement
                         (kWh/day) for tapping profile used for test
                     - energy_input_measured -- electrical input energy (kWh)
                         measured in EN 16147 test over 24 hrs
                     - power_standby -- standby power (kW) measured in EN 16147 test
                     - hw_vessel_loss_daily -- daily hot water vessel heat loss
                         (kWh/day) for a 45 K temperature difference between vessel
                         and surroundings, tested in accordance with BS 1566 or
                         EN 12897 or any equivalent standard. Vessel must be same
                         as that used during EN 16147 test
        vol_daily_average -- annual average hot water use for the dwelling, in litres / day
        tank_volume -- volume of tank in litres
        daily_losses --daily losses in kWh/day
        heat_exchanger_surface_area -- surface area of heat exchanger in m2 
        in_use_factor_mismatch -- in use factor to be applied to heat pump efficiency
        tank_volume_declared -- tank volume stored in the database in litres
        heat_exchanger_surface_area_declared 
            -- surface area of heat exchanger stored in the database in m2
        daily_losses_declared -- standing heat loss in kWh/day
        energy_supply_conn -- reference to EnergySupplyConnection object
        simulation_time -- reference to SimulationTime object
        control -- reference to a control object which must implement is_on() func
        """

        self.__pwr = power_max
        self.__energy_supply_conn = energy_supply_conn
        self.__simulation_time = simulation_time
        self.__control = control
        self.__tank_volume = tank_volume
        self.__daily_losses = daily_losses
        self.__heat_exchanger_surface_area = heat_exchanger_surface_area
        self.__in_use_factor_mismatch = in_use_factor_mismatch
        self.__tank_volume_declared = tank_volume_declared
        self.__heat_exchanger_surface_area_declared = heat_exchanger_surface_area_declared
        self.__daily_losses_declared = daily_losses_declared
        self.__vol_daily_average = vol_daily_average

        def init_efficiency_tapping_profile(
                cop_dhw,
                hw_tapping_prof_daily_total,
                energy_input_measured,
                power_standby,
                hw_vessel_loss_daily,
                ):
            """ Calculate efficiency for given test condition (tapping profile) """
            # CALCM-01 - DAHPSE - V2.0_DRAFT13, section 4.2
            temp_factor = 0.6 * 0.9
            energy_input_hw_vessel_loss = hw_vessel_loss_daily / cop_dhw * temp_factor
            energy_input_standby = power_standby * hours_per_day * temp_factor
            energy_input_test \
                = energy_input_measured - energy_input_standby + energy_input_hw_vessel_loss
            energy_demand_test = hw_tapping_prof_daily_total + hw_vessel_loss_daily * temp_factor
            return energy_demand_test / energy_input_test

        # Calculate efficiency for each tapping profile
        # TODO Check that expected tapping profiles have been provided
        efficiencies = {}
        for profile_name, profile_data in test_data.items():
            efficiencies[profile_name] = init_efficiency_tapping_profile(
                profile_data['cop_dhw'],
                profile_data['hw_tapping_prof_daily_total'],
                profile_data['energy_input_measured'],
                profile_data['power_standby'],
                profile_data['hw_vessel_loss_daily'],
                )

        def init_efficiency():
            """ Calculate initial efficiency based on SAP 10.2 section N3.7 b) and c) """
            if len(efficiencies) == 1 and 'M' in efficiencies.keys():
                # If efficiency for tapping profile M only has been provided, use it
                eff = efficiencies['M']
            elif len(efficiencies) == 2 \
            and 'M' in efficiencies.keys() \
            and 'L' in efficiencies.keys():
                # If efficiencies for tapping profiles M and L have been provided, interpolate
                vol_daily_limit_lower = 100.2
                vol_daily_limit_upper = 199.8
                if self.__vol_daily_average <= vol_daily_limit_lower :
                    eff = efficiencies['M']
                elif self.__vol_daily_average >= vol_daily_limit_upper :
                    eff = efficiencies['L']
                else:
                    eff_M = efficiencies['M']
                    eff_L = efficiencies['L']
                    eff = eff_M + (eff_L - eff_M) \
                        / (vol_daily_limit_upper - vol_daily_limit_lower) \
                        * (vol_daily_average - vol_daily_limit_lower)
            else:
                sys.exit('Unrecognised combination of tapping profiles in test data')
            return eff

        self.__initial_efficiency = init_efficiency()

    def calc_efficiency(self):
        """
        Calculate efficiency after applying in use factor if entered tank characteristics
        do not meet criteria of data in database.
        """
        if  self.__tank_volume < self.__tank_volume_declared or \
            self.__heat_exchanger_surface_area < self.__heat_exchanger_surface_area_declared or \
            self.__daily_losses > self.__daily_losses_declared:
            # heat pump does not meet criteria then in use factor applied
            in_use_factor_mismatch = self.__in_use_factor_mismatch
        else:
            # no in use factor is applied
            in_use_factor_mismatch = 1.0

        return self.__initial_efficiency * in_use_factor_mismatch

    def demand_energy(self, energy_demand, temp_return):
        """ Demand energy (in kWh) from the heat pump """
        # Account for time control where present. If no control present, assume
        # system is always active (except for basic thermostatic control, which
        # is implicit in demand calculation).
        if self.__control is None or self.__control.is_on():
            # Energy that heater is able to supply is limited by power rating
            energy_supplied = min(energy_demand, self.__pwr * self.__simulation_time.timestep())
        else:
            energy_supplied = 0.0

        energy_required = energy_supplied / self.calc_efficiency()
        self.__energy_supply_conn.demand_energy(energy_required)
        return energy_supplied

    def energy_output_max(self, temp_return):
        """ Calculate the maximum energy output (in kWh) from the heater """

        # Account for time control where present. If no control present, assume
        # system is always active (except for basic thermostatic control, which
        # is implicit in demand calculation).
        if self.__control is None or self.__control.is_on():
            # Energy that heater is able to supply is limited by power rating
            energy_max = self.__pwr * self.__simulation_time.timestep()
        else:
            energy_max = 0.0

        return energy_max
