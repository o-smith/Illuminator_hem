#!/usr/bin/env python3

"""
This module provides object(s) to model the behaviour of electric storage heaters.
"""

# Third-party imports
import sys
from enum import Enum, auto
from scipy.integrate import solve_ivp
import numpy as np
from scipy.integrate._ivp.ivp import OdeResult
import types
from typing import Union

# Local imports
import core.units as units
from core.space_heat_demand.zone import Zone
from core.energy_supply.energy_supply import EnergySupplyConnection
from core.simulation_time import SimulationTime
from core.controls.time_control import ToUChargeControl, SetpointTimeControl


class AirFlowType(Enum):
    FAN_ASSISTED = auto()
    DAMPER_ONLY = auto()

    @classmethod
    def from_string(cls, strval):
        if strval == 'fan-assisted':
            return cls.FAN_ASSISTED
        elif strval == 'damper-only':
            return cls.DAMPER_ONLY
        else:
            sys.exit('AirFlowType (' + str(strval) + ') not valid.')
            # TODO Exit just the current case instead of whole program entirely?


class ElecStorageHeater:
    """ Class to represent electric storage heaters """

    def __init__(
        self,
        rated_power: float,
        rated_power_instant: float,
        air_flow_type: str,
        temp_dis_safe: float,
        thermal_mass: float,
        frac_convective: float,
        U_ins: float,
        temp_charge_cut: float,
        mass_core: float,
        c_pcore: float,
        temp_core_target: float,
        A_core: float,
        c_wall: float,
        n_wall: float,
        thermal_mass_wall: float,
        fan_pwr: float,
        n_units: int,
        zone: Zone,
        energy_supply_conn: EnergySupplyConnection,
        simulation_time: SimulationTime,
        control: SetpointTimeControl,
        charge_control: ToUChargeControl
    ):
        """Construct an ElecStorageHeater object

        Arguments:
        rated_power          -- in kW (Charging)
        rated_power_instant  -- in kW (Instant backup)
        air_flow_type        -- string specifying type of Electric Storage Heater:
                             -- fan-assisted
                             -- damper-only
        temp_dis_safe        -- safe temperature to discharge hot air from device (60 degC)
        thermal_mass         -- thermal mass of emitters, in kWh / K
        frac_convective      -- convective fraction for heating (TODO: Check if necessary)
        U_ins                -- U-value insulation between core and case [W/m^2/K]
        temp_charge_cut      -- Room temperature at which, if sensed during a charging hour, the heater won't charge
        mass_core            -- kg mass core material [kg]
        c_pcore              -- thermal capacity of core material [J/kg/K]
        temp_core_target     -- target temperature for the core material on charging mode
                             -- this might include weather compensation with future more
                             -- advances controls
        A_core               -- Transfer area between the core and air [m2]
        c_wall               -- constant from characteristic equation of emitters (e.g. derived from BS EN 442 tests)
        n_wall               -- exponent from characteristic equation of emitters (e.g. derived from BS EN 442 tests)
        thermal_mass_wall    -- thermal mass of the case
        fan_pwr              -- Fan power [W]
        n_units              -- number of units install in zone
        zone                 -- zone where the unit(s) is/are installed
        energy_supply_conn   -- reference to EnergySupplyConnection object
        simulation_time      -- reference to SimulationTime object
        control              -- reference to a control object which must implement is_on() and setpnt() funcs
        """

        self.__pwr_in: float = (rated_power * units.W_per_kW)
        self.__pwr_instant: float = (rated_power_instant * units.W_per_kW)
        self.__air_flow_type: str = AirFlowType.from_string(air_flow_type)
        self.__t_dis_safe: float = temp_dis_safe
        self.__thermal_mass: float = thermal_mass
        self.__frac_convective: float = frac_convective
        self.__Uins: float = U_ins
        self.__temp_charge_cut: float = temp_charge_cut
        self.__n_units: int = n_units
        self.__zone: Zone = zone
        self.__energy_supply_conn: EnergySupplyConnection = energy_supply_conn
        self.__simtime: SimulationTime = simulation_time
        self.__control: SetpointTimeControl = control
        self.__charge_control: ToUChargeControl = charge_control
        self.__mass: float = mass_core  # 180.0  # kg of core
        self.__c_pcore: float = c_pcore  # 920.0  # J/kg/K core material specific heat
        self.__t_core_target: float = temp_core_target
        self.__A: float = A_core  # 4.0  # m^2 transfer area between core and case or wall
        self.__c: float = c_wall
        self.__n: float = n_wall
        self.__thermal_mass_wall = thermal_mass_wall
        self.__fan_pwr = fan_pwr

        # Power for driving fan
        # TODO: Modify fan energy calculation to SFP
        # self.__sfp: float = 0.0

        # Initialising other variables
        # Parameters

        self.__c_p: float = 1.0054  # J/kg/K air specific heat

        """
        The value of R (resistance of air gap) depends on several factors such as the thickness of the air gap, the
        temperature difference across the gap, and the air flow characteristics.
        Here are some typical values for R per inch of air gap thickness:
        Still air: 0.17 m²·K/W
        Air movement (average): 0.07 m²·K/W
        Air movement (high): 0.04 m²·K/W
        Note: These values are for a temperature difference of 24°F (14°C) and a pressure difference of 1 inch of water
        column. The values will change with changes in temperature and pressure differences.
        """

        # 0.17  # Thermal resistance of the air layer between the insulation and the wall of the device.
        # Assuming no change of resistance with temp diff. (Rough aprox)
        self.__Rair_off: float = 0.17
        self.__Rair_on: float = 0.07  # Same as above when the damper is on

        self.temp_air: float = self.__zone.temp_internal_air()  # °C Room temperature
        # case/wall c and n parameters as emitter.

        # This parameter specicify the opening ratio for the damper of the storage heater. It's set to 1.0 by
        # default but there might be control strategies using this to configure diferent levels of release
        self.damper_fraction: float = 1.0
        # labs_test for electric storage heater reaching 300 degC
        # This represents the temperature difference between the core and the room on the first column
        # and the fraction of air flow relating to the nominal as defined above on the second column
        self.labs_tests_400: list = [
            [85.07, 1.6],
            [91.65, 1.6],
            [98.73, 1.6],
            [106.36, 1.6],
            [114.57, 2.6],
            [123.42, 2.62],
            [133.25, 2.68],
            [144.29, 2.75],
            [156.79, 2.83],
            [171.03, 2.92],
            [187.41, 3.02],
            [206.41, 3.13],
            [228.63, 3.25],
            [254.93, 3.4],
            [286.52, 3.58],
            [324.91, 3.77]
        ]

        self.labs_tests_400_fan: list = [
            [0.0, 0.0],
            [8.31, 3.1],
            [21.1, 3.5],
            [34.84, 3.8],
            [49.84, 3.5],
            [106.53, 4.53],
            [235.47, 4.53],
            [347.15, 3.53],
            [463.39, 2.53],
            [584.73, 2.53],
            [713.24, 0.03]
        ]
        self.labs_tests_400_mod: list = [
            [2.62, 4.02],
            [3.49, 4.03],
            [4.66, 4.02],
            [6.22, 4.03],
            [8.31, 4.03],
            [11.1, 4.03],
            [14.84, 4.03],
            [19.84, 4.03],
            [26.53, 4.03],
            [35.47, 4.03],
            [47.42, 4.03],
            [63.39, 4.03],
            [84.73, 4.03],
            [113.24, 4.03],
            [151.33, 4.03],
            [202.2, 4.03],
            [270.15, 6.03],
            [470.15, 6.03],
            [570.15, 6.03],
            [670.15, 6.03]
        ]

        # This represents the temperature difference between the core and the room on the first column
        # and the fraction of air flow relating to the nominal as defined above on the second column
        self.labs_tests: list
        if self.__air_flow_type == AirFlowType.FAN_ASSISTED:
            self.labs_tests = self.labs_tests_400_fan
        elif self.__air_flow_type == AirFlowType.DAMPER_ONLY:
            self.labs_tests = self.labs_tests_400
        else:
            sys.exit('AirFlowType does not have characteristic data for EHS system')

        # Initial conditions
        self.t_core: float = 200.0 # self.__zone.temp_internal_air()
        self.t_wall: float = 50.0 #self.__zone.temp_internal_air()

        self.damper_fraction: float = 1.0
        self.__energy_in: float = 0.0

    def temp_setpnt(self):
        return self.__control.setpnt()

    def in_required_period(self):
        return self.__control.in_required_period()

    def frac_convective(self):
        return self.__frac_convective

    def __convert_to_kwh(self, power: float, timestep: float) -> float:
        """
        Converts power value supplied to the correct energy unit
        Arguments
        power -- Power value in watts
        timestep -- length of the timestep

        returns -- Energy in kWh
        """
        return power / units.W_per_kW * timestep * self.__n_units
    
    def __temp_charge_cut_corr(self) -> float:
        """
        Correct nominal/json temp_charge_cut with monthly table 
        Arguments

        returns -- temp_charge_cut (corrected)
        """
        temp_charge_cut_delta = [ -1.2, -0.6, 0.0, 0.6, 1.2, 1.2, 1.2, 1.2, 0.6, 0.0, -0.6, -1.2 ]
        current_month = self.__simtime.current_month()
        temp_charge_cut = self.__temp_charge_cut + temp_charge_cut_delta[current_month]

        return temp_charge_cut

    def __electric_charge(self, time: float, t_core: float) -> float:
        """
        Calculates power required for unit
        Arguments
        time -- current time period that we are looking at
        t_core -- current temperature of the core

        returns -- Power required in watts
        """

        if self.temp_air >= self.__temp_charge_cut_corr():
            return 0.0
        
        target_charge: float = self.__charge_control.target_charge()
        if self.__charge_control.is_on() and t_core <= self.__t_core_target * target_charge:
            pwr_required: float = (self.__t_core_target * target_charge - t_core) \
                                  * self.__mass * self.__c_pcore / self.__simtime.timestep()
            if pwr_required > self.__pwr_in:
                return self.__pwr_in
            else:
                return pwr_required
        else:
            return 0.0

    def __lab_test_ha(self, t_core_rm_diff: float) -> float:
        # labs_test for electric storage heater
        x: list = [row[0] for row in self.labs_tests]
        y: list = [row[1] for row in self.labs_tests]
        return np.interp(t_core_rm_diff, x, y)

    def __calulate_q_dis(self, time: float, t_core: float, q_out_wall: float, q_dis_modo: str) -> float:
        q_dis: float
        if q_dis_modo == "max":
            q_dis = self.__lab_test_ha(t_core - self.temp_air) * (t_core - self.temp_air)

        elif q_dis_modo == 0:
            q_dis = q_dis_modo

        else:
            # Equation for heat transfer to room from wall/case of heater
            q_dis = q_dis_modo - q_out_wall

        return q_dis

    def __return_q_released(self,
                            new_temp_core_and_wall: list,
                            q_released: float,
                            q_dis: float,
                            q_in: float,
                            timestep: int,
                            time: float,
                            q_instant: float = 0.0) -> float:
        # Setting core and wall temperatures to new values, for next iteration
        self.t_core = new_temp_core_and_wall[0]
        self.t_wall = new_temp_core_and_wall[1]

        # the purpose of this calculation is to calculate fan energy required by the device
        energy_for_fan_kwh: float = 0.0
        power_for_fan: float = 0.0
        if self.__air_flow_type == AirFlowType.FAN_ASSISTED and q_dis > 0:
            # TODO: Calculate fan energy using SFP
            # Air flow rate in kg/s
            # mass_flow_air = q_dis * self.__frac_convective / (self.__c_p * (self.__t_dis_safe - self.temp_air))
            # TODO: Improve mass to volumen conversion for air stream at 60 degC approx
            # Currently assumed air density of 1.04 kg/m3
            # Air flow volumen in m3/s
            # vol_flow_air = mass_flow_air / 1.04
            # print("vol: ", vol_flow_air)
            # power_for_fan: float = vol_flow_air * self.__sfp * units.W_per_kW
            # TODO: consider fan energy as part of q_dis and modify control to reduce q_dis accordingly
            power_for_fan: float = self.__fan_pwr
            energy_for_fan_kwh = self.__convert_to_kwh(power=power_for_fan, timestep=timestep)

        # Convert values to correct kwh unit
        q_in_kwh: float = self.__convert_to_kwh(power=q_in, timestep=timestep)
        q_instant_kwh: float = self.__convert_to_kwh(power=q_instant, timestep=timestep)

        # Save demand energy
        self.__energy_supply_conn.demand_energy(q_in_kwh + energy_for_fan_kwh + q_instant_kwh)


        # Multipy energy released by number of devices installed in the zone
        return self.__convert_to_kwh(power=(q_released + q_instant), timestep=timestep)

    def __heat_balance(self, temp_core_and_wall: list, time: float, q_dis_modo=0) -> tuple:
        """
        Calculates heat balance
        """
        t_core: float = temp_core_and_wall[0]
        t_wall: float = temp_core_and_wall[1]
        # Equation for electric charging
        q_in: float = self.__electric_charge(time, t_core)
        self.__energy_in = q_in

        q_dis: float

        q_out_wall: float
        if t_wall >= self.temp_air:
            q_out_wall = self.__c * (t_wall - self.temp_air) ** self.__n
        else:
            q_out_wall = self.__c * (t_wall - self.temp_air)

        # Equation for calculating q_dis
        q_dis = self.__calulate_q_dis(time=time, t_core=t_core, q_out_wall=q_out_wall, q_dis_modo=q_dis_modo)

        # Calculation of the U value between core and wall/case as
        # U value for the insulation and resistance of the air layer between the insulation and the wall/case
        if q_dis > 0:
            insulation: float = 1 / (1 / self.__Uins + self.__Rair_on)
        else:
            insulation: float = 1 / (1 / self.__Uins + self.__Rair_off)

        # Equation for the heat transfer between the core and the wall/case of the heater
        q_out_ins: float = insulation * self.__A * (t_core - t_wall)

        # Variation of Core temperature as per heat balance inside the heater
        dT_core: float = (1 / (self.__mass * self.__c_pcore)) * (q_in - q_out_ins - q_dis)

        # Variation of Wall/case temperature as per heat balance in surface of the heater
        dT_wall: float = (1 / self.__thermal_mass_wall) * (q_out_ins - q_out_wall)

        q_released: float = q_dis + q_out_wall

        return [dT_core, dT_wall], q_released, q_dis, q_in

    def __func_core_temperature_change_rate(self, q_dis_modo: Union[str, float]) -> types.FunctionType:
        """
        Lambda function for differentiation
        """
        return lambda time, t_core_and_wall: self.__heat_balance(temp_core_and_wall=t_core_and_wall,
                                                                 time=time,
                                                                 q_dis_modo=q_dis_modo)[0]

    def __calculate_sol_and_q_released(self,
                                       time_range: list,
                                       temp_core_and_wall: list,
                                       q_dis_modo: Union[str, float]) -> tuple:

        # first calculate how much the system is leaking without active discharging
        sol: OdeResult = solve_ivp(fun=self.__func_core_temperature_change_rate(q_dis_modo=q_dis_modo),
                                   t_span=time_range,
                                   y0=temp_core_and_wall,
                                   method='BDF')

        new_temp_core_and_wall: list = sol.y[:, -1]

        values: tuple = self.__heat_balance(temp_core_and_wall=new_temp_core_and_wall,
                                            time=time_range[1],
                                            q_dis_modo=q_dis_modo)
        q_released: float = values[1]
        q_dis: float = values[2]
        q_in: float = values[3]

        return new_temp_core_and_wall, q_released, q_dis, q_in

    def demand_energy(self, energy_demand: float) -> float:

        # Initialising Variables
        self.temp_air = self.__zone.temp_internal_air()
        timestep: float = self.__simtime.timestep()
        self.__time_unit: float = 3600 * timestep
        current_hour: int = self.__simtime.current_hour()
        time_range: list = [current_hour * self.__time_unit, (current_hour + 1) * self.__time_unit]
        temp_core_and_wall: list = [self.t_core, self.t_wall]

        # Converting energy_demand from kWh to Wh and distributing it through all units
        energy_demand: float = energy_demand * units.W_per_kW / self.__n_units
        

        #################################################
        # Step 1                                        #
        #################################################
        # first calculate how much the system is leaking without active discharging
        new_temp_core_and_wall: list
        q_released: float
        q_dis: float
        q_in: float
        new_temp_core_and_wall, q_released, q_dis, q_in = \
            self.__calculate_sol_and_q_released(time_range=time_range,
                                                temp_core_and_wall=temp_core_and_wall,
                                                q_dis_modo=0)

        # if Q_released is more than what the zone wants, that's it
        if q_released >= energy_demand / timestep:
            # More energy than needed to be released. End of calculations.
            return self.__return_q_released(new_temp_core_and_wall=new_temp_core_and_wall,
                                            q_released=q_released,
                                            q_dis=q_dis,
                                            q_in=q_in,
                                            timestep=timestep,
                                            time=time_range[1])

        #################################################
        # Step 2                                        #
        #################################################
        # Zone needs more than leaked, let's calculate with max discharging capability
        new_temp_core_and_wall, q_released, q_dis, q_in = \
            self.__calculate_sol_and_q_released(time_range=time_range,
                                                temp_core_and_wall=temp_core_and_wall,
                                                q_dis_modo="max")

        # If Q_released is not sufficient for zone demand, that's it
        # unless there is instnat backup that can top up the energy provided
        if q_released < energy_demand / timestep:
            if self.__pwr_instant > 0:
                power_supplied_instant = min(energy_demand / timestep - q_released, self.__pwr_instant)
            else:
                power_supplied_instant = 0.0

            # The system can only discharge the maximum amount, zone doesn't get everything it needs
            return self.__return_q_released(new_temp_core_and_wall=new_temp_core_and_wall,
                                            q_released=q_released,
                                            q_dis=q_dis,
                                            q_in=q_in,
                                            timestep=timestep,
                                            time=time_range[1],
                                            q_instant=power_supplied_instant)

        #################################################
        # Step 3                                        #
        #################################################
        # Zone actually needs an amount of energy that can be released by the system:
        # Let's call the heat balance forcing that amount (assuming perfect damper or
        # fan assisted control of the unit)
        q_dis = energy_demand / timestep
        new_temp_core_and_wall, q_released, q_dis, q_in = \
            self.__calculate_sol_and_q_released(time_range=time_range,
                                                temp_core_and_wall=temp_core_and_wall,
                                                q_dis_modo=q_dis)

        return self.__return_q_released(new_temp_core_and_wall=new_temp_core_and_wall,
                                        q_released=q_released,
                                        q_dis=q_dis,
                                        q_in=q_in,
                                        timestep=timestep,
                                        time=time_range[1])
