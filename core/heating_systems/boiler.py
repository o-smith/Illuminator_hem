#!/usr/bin/env python3

"""
This module provides objects to represent boilers. The calculations are based
on the hourly method developed for calculating efficiency improvements for
different classes of boiler controls in SAP 10 Appendix D2.2 (described in the
SAP 10 technical paper S10TP-12). This in turn was based on a combination of
the DAHPSE method for heat pumps (itself based on a draft of
BS EN 15316-4-2:2017 and described in the SAP calculation method CALCM-01) and
the Energy Balance Validation (EBV) method (described in the SAP 2009 technical
paper STP09/B02).
"""

# Standard library imports
import sys
from enum import Enum, auto

#Local imports
from core.energy_supply.energy_supply import Fuel_code
from core.material_properties import WATER
import core.units as units
from numpy import interp

class ServiceType(Enum):
    WATER_COMBI = auto()
    WATER_REGULAR = auto()
    SPACE = auto()


class Boiler_HW_test(Enum):
    M_L = auto()
    M_S = auto()
    M_only = auto()
    No_additional_tests = auto()

    @classmethod
    def from_string(cls, strval):
        if strval == 'M&L':
            return cls.M_L
        elif strval == 'M&S':
            return cls.M_S
        elif strval == 'M_only':
            return cls.M_only
        elif strval == 'No_additional_tests':
            return cls.No_additional_tests
        else:
            sys.exit('Hot water test ('+ str(strval) + ') not valid')

class BoilerService:
    """ A base class for objects representing services (e.g. water heating) provided by a boiler.

    This object encapsulates the name of the service, meaning that the system
    consuming the energy does not have to specify this on every call, and
    helping to enforce that each service has a unique name.

    Derived objects provide a place to handle parts of the calculation (e.g.
    distribution flow temperature) that may differ for different services.

    Separate subclasses need to be implemented for different types of service
    (e.g. HW and space heating). These should implement the following functions:
    - demand_energy(self, energy_demand)
    """

    def __init__(self, boiler, service_name, control=None):
        """ Construct a BoilerService object

        Arguments:
        boiler       -- reference to the Boiler object providing the service
        service_name -- name of the service demanding energy from the boiler
        control -- reference to a control object which must implement is_on() func
        """
        self._boiler = boiler
        self._service_name = service_name
        self.__control = control

    def is_on(self):
        if self.__control is not None:
            service_on = self.__control.is_on()
        else:
            service_on = True
        return service_on


class BoilerServiceWaterCombi(BoilerService):
    """ An object to represent a water heating service provided by a combi boiler.

    This object contains the parts of the boiler calculation that are
    specific to providing hot water.
    """

    def __init__(self,
                 boiler,
                 boiler_data,
                 service_name,
                 temp_hot_water,
                 cold_feed,
                 simulation_time
                ):
        """ Construct a BoilerServiceWaterCombi object

        Arguments:
        boiler       -- reference to the Boiler object providing the service
        boiler_data       -- combi boiler heating properties
        service_name -- name of the service demanding energy from the boiler_data
        temp_hot_water -- temperature of the hot water to be provided, in deg C
        cold_feed -- reference to ColdWaterSource object
        simulation_time -- reference to SimulationTime object
        """
        super().__init__(boiler, service_name)
        
        self.__temp_hot_water = temp_hot_water
        self.__cold_feed = cold_feed
        self.__service_name = service_name
        self.__simulation_time = simulation_time

        hw_tests = boiler_data["separate_DHW_tests"]
        self.__separate_DHW_tests = Boiler_HW_test.from_string(hw_tests)

        if (self.__separate_DHW_tests == Boiler_HW_test.M_L) \
            or (self.__separate_DHW_tests == Boiler_HW_test.M_S):
            #tapping cycle M and S, or M and L
            self.__rejected_energy_1 = boiler_data["rejected_energy_1"]
            self.__storage_loss_factor_2 = boiler_data["storage_loss_factor_2"]
            self.__rejected_factor_3 = boiler_data["rejected_factor_3"]

        elif self.__separate_DHW_tests == Boiler_HW_test.M_only:
            #tapping cycle M only test results
            self.__rejected_energy_1 = boiler_data["rejected_energy_1"]
            self.__storage_loss_factor_2 = boiler_data["storage_loss_factor_2"]

        #TODO feed in actual daily HW usage
        self.__daily_HW_usage = boiler_data["daily_HW_usage"]

    def get_cold_water_source(self):
        return self.__cold_feed

    def get_temp_hot_water(self):
        return self.__temp_hot_water

    def demand_hot_water(self, volume_demanded_target):
        """ Demand volume from boiler. Currently combi only """
        timestep = self.__simulation_time.timestep()
        return_temperature = 60 
        
        energy_content_kWh_per_litre = WATER.volumetric_energy_content_kWh_per_litre(
            self.__temp_hot_water,
            self.__cold_feed.temperature()
            )
        
        if 'temp_hot_water' in volume_demanded_target:
            volume_demanded = volume_demanded_target['temp_hot_water']['warm_vol']
        else:
            volume_demanded = 0.0

        energy_demand = volume_demanded * energy_content_kWh_per_litre 

        combi_loss = self.boiler_combi_loss(energy_demand, timestep)
        energy_demand = energy_demand + combi_loss

        return self._boiler._Boiler__demand_energy(
            self.__service_name,
            ServiceType.WATER_COMBI,
            energy_demand,
            return_temperature
            )

    def boiler_combi_loss(self, energy_demand, timestep):
        # daily hot water usage factor
        fu = 1.0
        threshold_volume = 100 # litres/day
        if self.__daily_HW_usage < threshold_volume:
            fu = self.__daily_HW_usage / threshold_volume

        # Equivalent hot water litres at 60C for HW load profiles
        hw_litres_S_profile = 36.0
        hw_litres_M_profile = 100.2
        hw_litres_L_profile = 199.8

        daily_vol_factor = hw_litres_M_profile - self.__daily_HW_usage
        if self.__separate_DHW_tests == Boiler_HW_test.M_S \
            and self.__daily_HW_usage < hw_litres_S_profile:
            daily_vol_factor = 64.2
        elif (self.__separate_DHW_tests == Boiler_HW_test.M_L \
            and self.__daily_HW_usage < hw_litres_M_profile) \
            or (self.__separate_DHW_tests == Boiler_HW_test.M_S \
            and self.__daily_HW_usage > hw_litres_M_profile):
            daily_vol_factor = 0
        elif self.__separate_DHW_tests == Boiler_HW_test.M_L \
            and self.__daily_HW_usage > hw_litres_L_profile:
            daily_vol_factor = -99.6

        combi_loss = 0.0
        if (self.__separate_DHW_tests == Boiler_HW_test.M_L) \
            or (self.__separate_DHW_tests == Boiler_HW_test.M_S):
            #combi loss calculation with tapping cycle M and S, or M and L
            combi_loss = (energy_demand * \
                          (self.__rejected_energy_1 + daily_vol_factor * self.__rejected_factor_3)) * fu \
                          + self.__storage_loss_factor_2 * (timestep / units.hours_per_day)

        elif self.__separate_DHW_tests == Boiler_HW_test.M_only:
            #combi loss calculation with tapping cycle M only test results
            combi_loss = (energy_demand * (self.__rejected_energy_1)) * fu \
                + self.__storage_loss_factor_2 * (timestep / units.hours_per_day)

        elif self.__separate_DHW_tests == Boiler_HW_test.No_additional_tests:
            # when no additional hot water test has been done
            default_combi_loss = 600 # annual default (kWh/day)
            combi_loss = default_combi_loss / units.days_per_year \
                        * (timestep / units.hours_per_day)

        else:
            exit('Invalid hot water test option')

        self.__combi_loss = combi_loss
        return combi_loss

    def internal_gains(self):
        # TODO Fraction of hot water energy resulting in internal gains should
        #      ideally be defined in one place, but it is duplicated here and in
        #      main hot water demand calculation for now.
        frac_dhw_energy_internal_gains = 0.25
        gain_internal \
            = frac_dhw_energy_internal_gains * self.__combi_loss \
            * units.W_per_kW / self.__simulation_time.timestep()
        self.__combi_loss = 0.0
        return gain_internal

    def energy_output_max(self):
        """ Calculate the maximum energy output of the boiler"""
        return self._boiler._Boiler__energy_output_max(self.__temp_hot_water)


class BoilerServiceWaterRegular(BoilerService):
    """ An object to represent a water heating service provided by a regular boiler.

    This object contains the parts of the boiler calculation that are
    specific to providing hot water.
    """

    def __init__(self,
                 boiler,
                 service_name,
                 temp_hot_water,
                 cold_feed,
                 simulation_time,
                 control=None,
                ):
        """ Construct a BoilerServiceWaterRegular object

        Arguments:
        boiler       -- reference to the Boiler object providing the service
        service_name -- name of the service demanding energy from the boiler_data
        temp_hot_water -- temperature of the hot water to be provided, in deg C
        cold_feed -- reference to ColdWaterSource object
        simulation_time -- reference to SimulationTime object
        control -- reference to a control object which must implement is_on() func
        """
        super().__init__(boiler, service_name, control)
        
        self.__temp_hot_water = temp_hot_water
        self.__cold_feed = cold_feed
        self.__service_name = service_name
        self.__simulation_time = simulation_time


    def demand_energy(
            self,
            energy_demand,
            temp_return,
            hybrid_service_bool=False,
            time_elapsed_hp = None,
            ):
        """ Demand energy (in kWh) from the boiler """
        service_on = self.is_on()
        if not service_on:
            energy_demand = 0.0
        
        return self._boiler._Boiler__demand_energy(
            self.__service_name,
            ServiceType.WATER_REGULAR,
            energy_demand,
            temp_return,
            hybrid_service_bool,
            time_elapsed_hp,
            )


    def energy_output_max(self, temp_return, time_elapsed_hp = None):
        """ Calculate the maximum energy output of the boiler"""
        service_on = self.is_on()
        if not service_on:
            return 0.0

        return self._boiler._Boiler__energy_output_max(self.__temp_hot_water, time_elapsed_hp)


class BoilerServiceSpace(BoilerService):
    """ An object to represent a space heating service provided by a boiler to e.g. a cylinder.

    This object contains the parts of the boiler calculation that are
    specific to providing space heating-.
    """
    def __init__(self, boiler, service_name, control):
        """ Construct a BoilerServiceSpace object

        Arguments:
        boiler       -- reference to the Boiler object providing the service
        service_name -- name of the service demanding energy from the boiler
        control -- reference to a control object which must implement is_on() and setpnt() funcs
        """
        super().__init__(boiler, service_name, control)
        self.__service_name = service_name
        self.__control = control

    def temp_setpnt(self):
        return self.__control.setpnt()

    def in_required_period(self):
        return self.__control.in_required_period()

    def demand_energy(
            self,
            energy_demand,
            temp_flow,
            temp_return,
            hybrid_service_bool=False,
            time_elapsed_hp = None,
            ):
        """ Demand energy (in kWh) from the boiler """
        if not self.is_on():
            if not hybrid_service_bool:
                # If not a hybrid hp system
                return 0.0
            else:
                # If hybrid hp system is off
                return 0.0, 0.0

        return self._boiler._Boiler__demand_energy(
            self.__service_name,
            ServiceType.SPACE,
            energy_demand,
            temp_return,
            hybrid_service_bool,
            time_elapsed_hp,
            )

    def energy_output_max(self, temp_output, temp_return_feed, time_elapsed_hp = None):
        """ Calculate the maximum energy output of the boiler"""
        if not self.is_on():
            return 0.0

        return self._boiler._Boiler__energy_output_max(temp_output, time_elapsed_hp)


class Boiler:
    """ An object to represent a boiler """

    def __init__(self, 
                boiler_dict,
                energy_supply,
                energy_supply_conn_aux,
                simulation_time,
                ext_cond, 
                ):
        """ Construct a Boiler object

        Arguments:
        boiler_dict -- dictionary of boiler characteristics, with the following elements:
            TODO List elements and their definitions
        energy_supply -- reference to EnergySupply object
        simulation_time -- reference to SimulationTime object
        external_conditions -- reference to ExternalConditions object
        
        Other variables:
        energy_supply_connections
            -- dictionary with service name strings as keys and corresponding
               EnergySupplyConnection objects as values
        """
        self.__energy_supply = energy_supply
        self.__simulation_time = simulation_time
        self.__external_conditions = ext_cond
        self.__energy_supply_connections = {}
        self.__energy_supply_connection_aux = energy_supply_conn_aux
        self.__service_results = []

        # boiler properties
        self.__boiler_location = boiler_dict["boiler_location"]
        self.__min_modulation_load = boiler_dict["modulation_load"]
        if self.__min_modulation_load > 1:
            sys.exit('Minimum modulation ratio cannot be greater than 1')
        self.__boiler_power = boiler_dict["rated_power"]
        full_load_gross = boiler_dict["efficiency_full_load"]
        part_load_gross = boiler_dict["efficiency_part_load"]
        self.__fuel_code = self.__energy_supply.fuel_type()

        # electricity properties
        self.__power_circ_pump = boiler_dict["electricity_circ_pump"]
        self.__power_part_load = boiler_dict["electricity_part_load"]
        self.__power_full_load = boiler_dict["electricity_full_load"]
        self.__power_standby = boiler_dict["electricity_standby"]
        self.__total_time_running_current_timestep = 0.0
        
        # high value correction 
        net_to_gross = self.net_to_gross()
        full_load_net = full_load_gross / net_to_gross
        part_load_net = part_load_gross / net_to_gross
        corrected_full_load_net = self.high_value_correction_full_load(full_load_net)
        corrected_part_load_net = self.high_value_correction_part_load(part_load_net)
        self.__corrected_full_load_gross = corrected_full_load_net * net_to_gross
        corrected_part_load_gross = corrected_part_load_net * net_to_gross

        #SAP model properties
        self.__room_temp = 19.5 #TODO use actual room temp instead of hard coding

        #30 is the nominal temperature difference between boiler and test room 
        #during standby loss test (EN15502-1 or EN15034)
        self.__temp_rise_standby_loss = 30.0
        #boiler standby heat loss power law index
        self.__sby_loss_idx = 1.25 

        #Calculate offset for EBV curves
        average_measured_eff = (corrected_part_load_gross + self.__corrected_full_load_gross) / 2.0
        # test conducted at return temperature 30C
        temp_part_load_test = 30.0 
        # test conducted at return temperature 60C
        temp_full_load_test = 60.0 
        offset_for_theoretical_eff = 0.0
        theoretical_eff_part_load = self.effvsreturntemp(temp_part_load_test, \
                                                         offset_for_theoretical_eff)
        theoretical_eff_full_load = self.effvsreturntemp(temp_full_load_test, \
                                                         offset_for_theoretical_eff)
        average_theoretical_eff = (theoretical_eff_part_load + theoretical_eff_full_load)/ 2.0
        self.__offset = average_theoretical_eff - average_measured_eff 

    def __create_service_connection(self, service_name):
        """ Create an EnergySupplyConnection for the service name given """
        # Check that service_name is not already registered
        if service_name in self.__energy_supply_connections.keys():
            sys.exit("Error: Service name already used: "+service_name)
            # TODO Exit just the current case instead of whole program entirely?

        # Set up EnergySupplyConnection for this service
        self.__energy_supply_connections[service_name] = \
            self.__energy_supply.connection(service_name)

    def create_service_hot_water_combi(
            self,
            boiler_data,
            service_name,
            temp_hot_water,
            cold_feed
            ):
        """ Return a BoilerServiceWater object and create an EnergySupplyConnection for it
        
        Arguments:
        boiler_data -- boiler hot water heating properties
        service_name -- name of the service demanding energy from the boiler
        temp_hot_water -- temperature of the hot water to be provided, in deg C
        cold_feed -- reference to ColdWaterSource object
        """
        self.__create_service_connection(service_name)
        return BoilerServiceWaterCombi(
            self,
            boiler_data,
            service_name,
            temp_hot_water,
            cold_feed,
            self.__simulation_time
            )
    
    def create_service_hot_water_regular(
            self,
            service_name,
            temp_hot_water,
            cold_feed,
            control=None,
            ):
            """ Return a BoilerServiceWaterRegular object and create an EnergySupplyConnection for it

            Arguments:
            service_name -- name of the service demanding energy from the boiler
            temp_hot_water -- temperature of the hot water to be provided, in deg C
            temp_limit_upper -- upper operating limit for temperature, in deg C
            cold_feed -- reference to ColdWaterSource object
            control -- reference to a control object which must implement is_on() func
            """
            
            self.__create_service_connection(service_name)
            return BoilerServiceWaterRegular(
                self,
                service_name,
                temp_hot_water,
                cold_feed,
                self.__simulation_time,
                control,
                )

    def create_service_space_heating(
            self,
            service_name,
            control,
            ):
        """ Return a BoilerServiceSpace object and create an EnergySupplyConnection for it

        Arguments:
        service_name -- name of the service demanding energy from the boiler
        control -- reference to a control object which must implement is_on() and setpnt() funcs
        """
        self.__create_service_connection(service_name)
        return BoilerServiceSpace(
            self,
            service_name,
            control,
            )

    def __cycling_adjustment(self, temp_return_feed, standing_loss, prop_of_timestep_at_min_rate, temp_boiler_loc):
        ton_toff = (1.0 - prop_of_timestep_at_min_rate) / prop_of_timestep_at_min_rate
        cycling_adjustment = standing_loss \
                             * ton_toff \
                             * ((temp_return_feed - temp_boiler_loc) \
                             / (self.__temp_rise_standby_loss) \
                             ) \
                             ** self.__sby_loss_idx
        
        return cycling_adjustment


    def location_adjustment(self, temp_return_feed, standing_loss, temp_boiler_loc):
        location_adjustment \
            = max((standing_loss * \
                    ((temp_return_feed - self.__room_temp))**self.__sby_loss_idx \
                    - (temp_return_feed - temp_boiler_loc)**self.__sby_loss_idx)\
                    , 0.0
                 )
        return location_adjustment

    def __calc_current_boiler_power(self, energy_output_provided, time_available):

        if time_available <= 0:
            return 0.0

        min_power = self.__boiler_power * self.__min_modulation_load

        current_boiler_power = max(energy_output_provided / time_available , min_power)

        return current_boiler_power

    def calc_boiler_eff(
            self,
            service_type,
            temp_return_feed,
            energy_output_required,
            time_elapsed_hp,
            ):

        time_available = self.__time_available(time_elapsed_hp)
        energy_output_provided = self.__calc_energy_output_provided(energy_output_required, time_available)

        current_boiler_power = self.__calc_current_boiler_power(energy_output_provided, time_available)

        # The efficiency of the boiler depends on whether it cycles on/off.
        # If this occurs, an adjustment is calculated for the calculation 
        # timestep as follows (when the boiler is firing continuously no 
        # adjustment is necessary so cycling_adjustment=0).
        if time_available <= 0.0:
            prop_of_timestep_at_min_rate = 0.0
        else:
            prop_of_timestep_at_min_rate = min(energy_output_required \
                               / (self.__boiler_power * self.__min_modulation_load * time_available)
                               ,1.0)

        # Default value for the stand-by heat losses as a function of the current boiler power
        # Equation 5 in EN15316-4-1
        # fgen = (c5*(Pn)^c6)/100
        # where c5 = 4.0, c6 = -0.4 and Pn is the current boiler power
        if current_boiler_power == 0.0:
            standing_loss = 0.0
        else:
            standing_loss = (4.0 * (current_boiler_power)** - 0.4) / 100.0 

        #use weather temperature at timestep
        outside_temp = self.__external_conditions.air_temp()

        # A boiler’s efficiency reduces when installed outside due to an increase in case heat loss.
        # The following adjustment is made when the boiler is located outside 
        # (when installed inside no adjustment is necessary so location_adjustment=0)
        if self.__boiler_location == "external":
            temp_boiler_loc = outside_temp
        elif self.__boiler_location == "internal":
            temp_boiler_loc = self.__room_temp
        else:
            sys.exit('boiler location ('+ str(self.__boiler_location) + ') not valid')

        # Calculate location adjustment
        location_adjustment = 0.0
        if self.__boiler_location == "external":
            location_adjustment = self.location_adjustment(temp_return_feed, standing_loss, 
                temp_boiler_loc)

        # Calculate cycling adjustment
        cycling_adjustment = 0.0
        if (0.0 < prop_of_timestep_at_min_rate < 1.0) and service_type != ServiceType.WATER_COMBI:
            cycling_adjustment = self.__cycling_adjustment(temp_return_feed, 
                standing_loss,
                prop_of_timestep_at_min_rate,
                temp_boiler_loc)

        # Calculate combined cyclic and location adjustment
        cyclic_location_adjustment = cycling_adjustment + location_adjustment

        # Calculate boiler efficiency based on the return temperature and offset
        boiler_eff = self.effvsreturntemp(temp_return_feed, self.__offset)

        # If boiler starts cycling use the corrected full load efficiency 
        # as the boiler eff before cycling adjustment is applied.
        if cycling_adjustment > 0.0:
            boiler_eff = self.__corrected_full_load_gross

        # Calculate the final boiler efficiency
        blr_eff_final = 1.0 / ((1.0 / boiler_eff) + cyclic_location_adjustment)

        return blr_eff_final


    def __calc_energy_output_provided(self, energy_output_required, time_available):
        energy_output_max_power = self.__boiler_power * time_available
        energy_output_provided = min(energy_output_required, energy_output_max_power)

        return energy_output_provided

    def __time_available(self, time_elapsed_hp = None):
        timestep = self.__simulation_time.timestep()
        if time_elapsed_hp is not None:
            time_available = timestep - time_elapsed_hp
        else:
            time_available = timestep - self.__total_time_running_current_timestep
        return time_available

    def __demand_energy(
            self,
            service_name,
            service_type,
            energy_output_required,
            temp_return_feed,
            hybrid_service_bool = False,
            time_elapsed_hp = None,
            ):
        """ Calculate energy required by boiler to satisfy demand for the service indicated."""
        # Account for time control where present. If no control present, assume
        # system is always active (except for basic thermostatic control, which
        # is implicit in demand calculation).
        # if self.__control is None or self.__control.is_on():
        #     # Energy that heater is able to supply is limited by power rating
        #     energy_output_provided = min(energy_output_required, self.__boiler_power * self.__simulation_time.timestep())
        # else:
        #     energy_output_provided = 0.0

        time_available = self.__time_available(time_elapsed_hp)
        
        # If there is no demand on the boiler or no remaining time then no energy should be provided
        if energy_output_required <= 0.0 or time_available <= 0.0:
            energy_output_provided = 0.0
            fuel_demand = 0.0
            self.__energy_supply_connections[service_name].demand_energy(fuel_demand)

        blr_eff_final = self.calc_boiler_eff(
            service_type,
            temp_return_feed,
            energy_output_required,
            time_elapsed_hp,
            )

        energy_output_provided = self.__calc_energy_output_provided(energy_output_required, time_available)
        current_boiler_power = self.__calc_current_boiler_power(energy_output_provided, time_available)

        fuel_demand = energy_output_provided / blr_eff_final

        self.__energy_supply_connections[service_name].demand_energy(fuel_demand)

        # Calculate running time of Boiler
        if current_boiler_power <= 0.0:
            time_running_current_service = 0.0
        else:
            time_running_current_service = min(
                energy_output_provided / current_boiler_power,
                time_available
                )
        self.__total_time_running_current_timestep += time_running_current_service
        
        # Save results that are needed later (in the timestep_end function)
        self.__service_results.append({
            'service_name': service_name,
            'time_running': time_running_current_service,
            'current_boiler_power': current_boiler_power
            })
        
        if hybrid_service_bool:
            return energy_output_provided, time_running_current_service
        else:
            return energy_output_provided
    
    def __calc_auxiliary_energy(self, timestep, time_remaining_current_timestep):
        """Calculation of boiler electrical consumption"""

        #Energy used by circulation pump
        energy_aux = self.__total_time_running_current_timestep \
            * self.__power_circ_pump
        
        #Energy used in standby mode
        energy_aux += self.__power_standby * time_remaining_current_timestep
        
        #Energy used by flue fan
        #Flue fan electricity for on-off boilers
        elec_energy_flue_fan = self.__total_time_running_current_timestep \
            * self.__power_full_load

        #Overwrite flue fan electricity if boiler modulates
        #TODO does cycling below part load decrease elec consumption
        for service_no, service_data in enumerate(self.__service_results):
            modulation_ratio = min(service_data['current_boiler_power'] / self.__boiler_power, 1.0)
            if self.__min_modulation_load < 1:
                x_axis = [0.3, 1.0]
                y_axis = [self.__power_part_load, self.__power_full_load]

                # Uses numpy interp
                flue_fan_el = interp(modulation_ratio, x_axis, y_axis)
                elec_energy_flue_fan = service_data['time_running'] \
                                     * flue_fan_el
                energy_aux += elec_energy_flue_fan

        self.__energy_supply_connection_aux.demand_energy(energy_aux)

    def timestep_end(self):
        """" Calculations to be done at the end of each timestep"""
        timestep = self.__simulation_time.timestep()
        time_remaining_current_timestep = timestep - self.__total_time_running_current_timestep

        self.__calc_auxiliary_energy(timestep, time_remaining_current_timestep)

        #Variabales below need to be reset at the end of each timestep
        self.__total_time_running_current_timestep = 0.0
        self.__service_results = []
        
    def __energy_output_max(self, temp_output, time_elapsed_hp = None):
        time_available = self.__time_available(time_elapsed_hp)
        return self.__boiler_power * time_available
    
    def effvsreturntemp(self, return_temp, offset):
        """ Return boiler efficiency at different return temperatures """
        mains_gas_dewpoint = 52.2
        lpg_dewpoint = 48.3
        #TODO: add remaining fuels 
        if self.__fuel_code == Fuel_code.MAINS_GAS:
            if return_temp < mains_gas_dewpoint:
                theoretical_eff = -0.00007 * (return_temp)**2 + 0.0017 * return_temp + 0.979 
            else:
                theoretical_eff = -0.0006 * return_temp + 0.9129
        elif (self.__fuel_code == Fuel_code.LPG_BULK) or \
             (self.__fuel_code == Fuel_code.LPG_BOTTLED) or \
             (self.__fuel_code == Fuel_code.LPG_CONDITION_11F):
            if return_temp < lpg_dewpoint:
                theoretical_eff = -0.00006 * (return_temp)**2 + 0.0013 * return_temp + 0.9859
            else:
                theoretical_eff = -0.0006 * return_temp + 0.933
        else:
            exit('Fuel code does not exist')
        blr_theoretical_eff = theoretical_eff - offset

        return blr_theoretical_eff

    def high_value_correction_part_load(self, net_efficiency_part_load):
        """ Return a Boiler efficiency corrected for high values """
        if self.__fuel_code == Fuel_code.MAINS_GAS:
            maximum_part_load_eff = 1.08
        elif (self.__fuel_code == Fuel_code.LPG_BULK) or \
             (self.__fuel_code == Fuel_code.LPG_BOTTLED) or \
             (self.__fuel_code == Fuel_code.LPG_CONDITION_11F):
            maximum_part_load_eff = 1.06
        else:
            exit('Unknown fuel code '+str(self.__fuel_code))
        corrected_net_efficiency_part_load = min(net_efficiency_part_load \
                                                 - 0.213 \
                                                 * (net_efficiency_part_load - 0.966), \
                                                 maximum_part_load_eff)
        return corrected_net_efficiency_part_load

    def high_value_correction_full_load(self, net_efficiency_full_load):
        corrected_net_efficiency_full_load = min(net_efficiency_full_load \
                                                 - 0.673 * (net_efficiency_full_load - 0.955), \
                                                 0.98)
        return corrected_net_efficiency_full_load

    def net_to_gross(self):
        """ Returns net to gross factor """
        if self.__fuel_code == Fuel_code.MAINS_GAS:
            net_to_gross = 0.901
        elif (self.__fuel_code == Fuel_code.LPG_BULK) or \
             (self.__fuel_code == Fuel_code.LPG_BOTTLED) or \
             (self.__fuel_code == Fuel_code.LPG_CONDITION_11F):
            net_to_gross = 0.921
        else:
            exit('Unknown fuel code '+str(self.__fuel_code))
        return net_to_gross
                

