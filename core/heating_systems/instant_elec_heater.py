#!/usr/bin/env python3

"""
This module provides object(s) to model the behaviour of instantaneous electric
room heaters.
"""

class InstantElecHeater:
    """ Class to represent instantaneous electric heaters """

    def __init__(
            self,
            rated_power,
            frac_convective,
            energy_supply_conn,
            simulation_time,
            control,
            ):
        """ Construct an InstantElecHeater object

        Arguments:
        rated_power        -- in kW
        frac_convective    -- convective fraction for heating
        energy_supply_conn -- reference to EnergySupplyConnection object
        simulation_time    -- reference to SimulationTime object
        control -- reference to a control object which must implement is_on() and setpnt() funcs
        """
        self.__pwr                = rated_power
        self.__frac_convective    = frac_convective
        self.__energy_supply_conn = energy_supply_conn
        self.__simulation_time    = simulation_time
        self.__control            = control

    def temp_setpnt(self):
        return self.__control.setpnt()

    def in_required_period(self):
        return self.__control.in_required_period()

    def frac_convective(self):
        return self.__frac_convective

    def demand_energy(self, energy_demand):
        """ Demand energy (in kWh) from the heater """

        # Account for time control where present. If no control present, assume
        # system is always active (except for basic thermostatic control, which
        # is implicit in demand calculation).
        # TODO Account for manual (or smart) control where heater may be left
        #      on for longer than it would be under a simple thermostatic
        #      control?
        if self.__control is None or self.__control.is_on():
            # Energy that heater is able to supply is limited by power rating
            energy_supplied = min(energy_demand, self.__pwr * self.__simulation_time.timestep())
        else:
            energy_supplied = 0.0

        self.__energy_supply_conn.demand_energy(energy_supplied)
        
        return energy_supplied
