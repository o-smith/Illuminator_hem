#!/usr/bin/env python3

"""
This module provides objects to air conditioning.
"""


class AirConditioning:
    """ An object to model an air conditioning system """

    def __init__(
            self,
            cooling_capacity,
            efficiency,
            frac_convective,
            energy_supply_conn,
            simulation_time,
            control,
            ):
        """ Construct an air conditioning object

        Arguments:
        cooling_capacity -- maximum cooling capacity of the system, in kW
        efficiency -- SEER
        frac_convective    -- convective fraction for cooling
        energy_supply_conn -- reference to EnergySupplyConnection object
        simulation_time    -- reference to SimulationTime object
        control -- reference to a control object which must implement is_on() and setpnt() funcs
        """
        self.__cooling_capacity = cooling_capacity
        self.__efficiency = efficiency
        self.__frac_convective = frac_convective
        self.__energy_supply_conn = energy_supply_conn
        self.__simulation_time = simulation_time
        self.__control = control

    def temp_setpnt(self):
        return self.__control.setpnt()

    def in_required_period(self):
        return self.__control.in_required_period()

    def frac_convective(self):
        return self.__frac_convective

    def demand_energy(self, cooling_demand):
        """ Demand energy (in kWh) from the cooling system """

        # Account for time control where present. If no control present, assume
        # system is always active (except for basic thermostatic control, which
        # is implicit in demand calculation).
        # TODO Account for manual (or smart) control where cooling system may be left
        #      on for longer than it would be under a simple thermostatic
        #      control?
        if self.__control is None or self.__control.is_on():
            # Energy that cooling system is able to supply is limited by power rating
            # Note: negate cooling capacity because cooling demand is negative by convention
            cooling_supplied = max(
                cooling_demand,
                - self.__cooling_capacity * self.__simulation_time.timestep()
                )
        else:
            cooling_supplied = 0.0

        # Negate cooling_supplied because cooling demand is negative by convention,
        # but demand on energy_supply_conn must be positive
        self.__energy_supply_conn.demand_energy(- cooling_supplied / self.__efficiency)

        return cooling_supplied