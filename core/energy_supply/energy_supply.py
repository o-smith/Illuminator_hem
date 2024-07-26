#!/usr/bin/env python3

"""
This module contains objects that represent energy supplies such as mains gas,
mains electricity or other fuels (e.g. LPG, wood pellets).
"""

# Standard library inputs
import sys
import numpy as np
from enum import Enum,auto

class Fuel_code(Enum):
    MAINS_GAS = auto()
    ELECTRICITY = auto()
    UNMET_DEMAND = auto()
    CUSTOM = auto()
    LPG_BULK = auto()
    LPG_BOTTLED = auto()
    LPG_CONDITION_11F = auto()

    @classmethod
    def from_string(cls, strval):
        if strval == 'mains_gas':
            return cls.MAINS_GAS
        elif strval == 'electricity':
            return cls.ELECTRICITY
        elif strval == 'unmet_demand':
            return cls.UNMET_DEMAND
        elif strval == 'custom':
            return cls.CUSTOM
        elif strval == 'LPG_bulk':
            return cls.LPG_BULK
        elif strval == 'LPG_bottled':
            return cls.LPG_BOTTLED
        elif strval == 'LPG_condition_11F':
            return cls.LPG_CONDITION_11F
        else:
            sys.exit('fuel code ('+ str(strval) + ') not valid')

class EnergySupplyConnection:
    """ An object to represent the connection of a system that consumes energy to the energy supply

    This object encapsulates the name of the connection, meaning that the
    system consuming the energy does not have to specify these on every call,
    and helping to enforce that each connection to a single supply has a unique
    name.
    """

    def __init__(self, energy_supply, end_user_name):
        """ Construct an EnergySupplyConnection object

        Arguments:
        energy_supply -- reference to the EnergySupply object that the connection is to
        end_user_name -- name of the system (and end use, where applicable)
                         consuming energy from this connection
        """
        self.__energy_supply = energy_supply
        self.__end_user_name = end_user_name

    def energy_out(self, amount_demanded):
        """ Forwards the amount of energy out (in kWh) to the relevant EnergySupply object """
        self.__energy_supply._EnergySupply__energy_out(self.__end_user_name, amount_demanded)

    def demand_energy(self, amount_demanded):
        """ Forwards the amount of energy demanded (in kWh) to the relevant EnergySupply object """
        self.__energy_supply._EnergySupply__demand_energy(self.__end_user_name, amount_demanded)

    def supply_energy(self, amount_produced):
        """ Forwards the amount of energy produced (in kWh) to the relevant EnergySupply object """
        self.__energy_supply._EnergySupply__supply_energy(self.__end_user_name, amount_produced)

    def fuel_type(self):
        return self.__energy_supply.fuel_type()

class EnergySupply:
    """ An object to represent an energy supply, and to report energy consumption """
    # TODO Do we need a subclass for electricity supply specifically, to
    #      account for generators? Or do we just handle it in this object and
    #      have an empty list of generators when not electricity?

    def __init__(self, fuel_type, simulation_time, elec_battery = None, priority = None, is_export_capable = True):
        """ Construct an EnergySupply object

        Arguments:
        fuel_type          -- string denoting type of fuel
                              TODO Consider replacing with fuel_type object
        simulation_time    -- reference to SimulationTime object
        elec_battery       -- reference to an ElectricBattery object
        is_export_capable  -- denotes that this Energy Supply can export its surplus supply

        Other variables:
        demand_total       -- list to hold total demand on this energy supply at each timestep
        demand_by_end_user -- dictionary of lists to hold demand from each end user on this
                              energy supply at each timestep
        """
        self.__fuel_type          = Fuel_code.from_string(fuel_type)
        self.__simulation_time    = simulation_time
        self.__elec_battery       = elec_battery
        self.__diverter           = None
        self.__priority           = priority
        self.__is_export_capable  = is_export_capable

        self.__demand_total       = self.__init_demand_list()
        self.__demand_by_end_user = {}
        self.__energy_out_by_end_user = {}
        self.__beta_factor = self.__init_demand_list() #this would be multiple columns if multiple beta factors
        self.__supply_surplus = self.__init_demand_list()
        self.__demand_not_met = self.__init_demand_list()
        self.__energy_into_battery = self.__init_demand_list()
        self.__energy_out_of_battery = self.__init_demand_list()
        self.__energy_diverted = self.__init_demand_list()
        self.__energy_generated_consumed = self.__init_demand_list()

    def __init_demand_list(self):
        """ Initialise zeroed list of demand figures (one list entry for each timestep) """
        # TODO Consider moving this function to SimulationTime object if it
        #      turns out to be more generally useful.
        return [0] * self.__simulation_time.total_steps()

    def connection(self, end_user_name):
        """ Return an EnergySupplyConnection object and initialise list for the end user demand """
        # Check that end_user_name is not already registered/connected
        if end_user_name in self.__demand_by_end_user.keys():
            sys.exit("Error: End user name already used: "+end_user_name)
            # TODO Exit just the current case instead of whole program entirely?

        self.__demand_by_end_user[end_user_name] = self.__init_demand_list()
        self.__energy_out_by_end_user[end_user_name] = self.__init_demand_list()
        return EnergySupplyConnection(self, end_user_name)

    def __energy_out(self, end_user_name, amount_demanded):
        # Check that end_user_name is already connected/registered
        if end_user_name not in self.__demand_by_end_user.keys():
            sys.exit("Error: End user name ("+end_user_name+
                     ") not already registered by calling connection function.")
            # TODO Exit just the current case instead of whole program entirely?

        t_idx = self.__simulation_time.index()
        self.__energy_out_by_end_user[end_user_name][t_idx] \
            = self.__energy_out_by_end_user[end_user_name][t_idx] \
            + amount_demanded

    def connect_diverter(self, diverter):
        if self.__diverter is not None:
            sys.exit('Diverter already connected.')
        self.__diverter = diverter

    def __demand_energy(self, end_user_name, amount_demanded):
        """ Record energy demand (in kWh) for the end user specified.

        Note: Call via an EnergySupplyConnection object, not directly.
        """
        # Check that end_user_name is already connected/registered
        if end_user_name not in self.__demand_by_end_user.keys():
            sys.exit("Error: End user name ("+end_user_name+
                     ") not already registered by calling connection function.")
            # TODO Exit just the current case instead of whole program entirely?

        t_idx = self.__simulation_time.index()
        self.__demand_total[t_idx] = self.__demand_total[t_idx] + amount_demanded
        self.__demand_by_end_user[end_user_name][t_idx] \
            = self.__demand_by_end_user[end_user_name][t_idx] \
            + amount_demanded

    def __supply_energy(self, end_user_name, amount_produced):
        """ Record energy produced (in kWh) for the end user specified.

        Note: this is energy generated so it is subtracted from demand.
        Treat as negative
        """
        #energy produced in kWh as 'negative demand'
        amount_produced = amount_produced * -1
        self.__demand_energy(end_user_name, amount_produced)

    def results_total(self):
        """ Return list of the total demand on this energy source for each timestep """
        return self.__demand_total

    def results_by_end_user(self):
        """ Return the demand from each end user on this energy source for each timestep.

        Returns dictionary of lists, where dictionary keys are names of end users.
        """
        # If the keys do match then we will just return the demand by end users
        if self.__demand_by_end_user.keys() != self.__energy_out_by_end_user.keys():
            return self.__demand_by_end_user

        all_results_by_end_user = {}
        for demand, energy_out in zip(self.__demand_by_end_user.items(), self.__energy_out_by_end_user.items()):
            if demand[0] == energy_out[0]:
                user_name = demand[0]  # Can use either demand[0] or energy_out[0] to retrieve end user name
                all_results_by_end_user[user_name] = np.array(demand[1]) + np.array(energy_out[1])

        return all_results_by_end_user

    def get_energy_import(self):
        return self.__demand_not_met

    def get_energy_export(self):
        return self.__supply_surplus

    def get_energy_generated_consumed(self):
        """ Return the amount of generated energy consumed in the building for all timesteps """
        return self.__energy_generated_consumed

    def get_energy_to_from_battery(self):
        """ Return the amount of generated energy sent to battery and drawn from battery """
        return self.__energy_into_battery, self.__energy_out_of_battery

    def get_energy_diverted(self):
        """ Return the amount of generated energy diverted to minimise export """
        return self.__energy_diverted

    def get_beta_factor(self):
        return self.__beta_factor

    def calc_energy_import_export_betafactor(self):
        """
        calculate how much of that supply can be offset against demand.
        And then calculate what demand and supply is left after offsetting, which are the amount exported imported
        """

        supplies=[]
        demands=[]
        t_idx = self.__simulation_time.index()
        for user in self.__demand_by_end_user.keys():
            demand = self.__demand_by_end_user[user][t_idx]
            if demand < 0.0:
                # if energy is negative that means its actually a supply, we
                # need to separate the two for beta factor calc. If we had
                # multiple different supplies they would have to be separated
                # here
                supplies.append(demand)
            else:
                demands.append(demand)

        self.__beta_factor[t_idx] = self.beta_factor_function(- sum(supplies), sum(demands), 'PV')

        # PV elec consumed within dwelling in absence of battery storage or diverter (kWh)
        # if there were multiple sources they would each have their own beta factors
        supply_consumed = sum(supplies) * self.__beta_factor[t_idx]
        # Surplus PV elec generation (kWh) - ie amount to be exported to the grid or batteries
        supply_surplus = sum(supplies) * (1 - self.__beta_factor[t_idx])
        # Elec demand not met by PV (kWh) - ie amount to be imported from the grid or batteries
        demand_not_met = sum(demands) + supply_consumed
        
        if self.__priority is None and self.__elec_battery is not None:
            #See if the battery can deal with excess supply/demand for this timestep
            #supply_surplus is -ve by convention and demand_not_met is +ve
            #TODO: assumption made here that supply is done before demand, could
            # revise in future if more evidence becomes available.
            energy_out_of_battery = self.__elec_battery.charge_discharge_battery(supply_surplus)
            supply_surplus -= energy_out_of_battery
            self.__energy_into_battery[t_idx] = - energy_out_of_battery
            energy_out_of_battery = self.__elec_battery.charge_discharge_battery(demand_not_met)
            demand_not_met -= energy_out_of_battery
            self.__energy_out_of_battery[t_idx] = - energy_out_of_battery
            
        if self.__priority is None and self.__diverter is not None:
            # Divert as much surplus energy as possible, and calculate remaining surplus
            self.__energy_diverted[t_idx] = self.__diverter.divert_surplus(supply_surplus)
            supply_surplus += self.__energy_diverted[t_idx]
        
        # If the priority order of energy surplus is specified,calculate the same according to the order
        if self.__priority is not None:
            for item in self.__priority:
                if item == 'ElectricBattery' and self.__elec_battery:
                    energy_out_of_battery = self.__elec_battery.charge_discharge_battery(supply_surplus)
                    supply_surplus -= energy_out_of_battery
                    self.__energy_into_battery[t_idx] = - energy_out_of_battery
                    energy_out_of_battery = self.__elec_battery.charge_discharge_battery(demand_not_met)
                    demand_not_met -= energy_out_of_battery
                    self.__energy_out_of_battery[t_idx] = - energy_out_of_battery   
                if item == 'diverter'and self.__diverter is not None:
                    self.__energy_diverted[t_idx] = self.__diverter.divert_surplus(supply_surplus)
                    supply_surplus += self.__energy_diverted[t_idx]
                    
        if self.__is_export_capable: self.__supply_surplus[t_idx] += supply_surplus
        self.__demand_not_met[t_idx] += demand_not_met
        # Report energy generated and consumed as positive number, so subtract negative number
        self.__energy_generated_consumed[t_idx] -= supply_consumed                        

    def beta_factor_function(self,supply,demand,beta_factor_function):
        """
        wrapper that applies relevant function to obtain
        beta factor from energy supply+demand at a given timestep
        """

        if supply == 0.0:
            beta_factor = 1.0
            return beta_factor

        if demand == 0.0:
            beta_factor = 0.0
            return beta_factor


        demand_ratio = float(supply) / float(demand)
        if beta_factor_function=='PV':
            # Equation for beta factor below is based on hourly data from four
            # dwellings, which gives a similar monthly beta factor to that
            # calculated from the beta factor equation in SAP 10.2, which was
            # based on monthly data from 15 dwellings.
            # TODO: come up with better fit curve for PV
            beta_factor = min(0.6748 *pow(demand_ratio,-0.703),1.0)
        # TODO
        # elif function=='wind':
        #     beta_factor=1.0
        else:
            sys.exit('Invalid value for beta_factor_function')

        """
        predicted beta should not be greater than 1/demand_ratio, otherwise
        we might predict demand fulfilled by PV/generation to be greater than
        total demand.
        """

        beta_factor = min(beta_factor,1/demand_ratio)

        return beta_factor

    def fuel_type(self):
        return self.__fuel_type
    
    