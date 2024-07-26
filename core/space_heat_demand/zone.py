#!/usr/bin/env python3

"""
This module provides objects to represent the thermal zones in the building,
and to calculate the temperatures in the zone and associated building elements.
"""

# Standard library imports
import sys

# Third-party imports
import numpy as np
import cython

# Local imports
import core.units as units
from core.material_properties import AIR
from core.space_heat_demand.building_element import \
    BuildingElementAdjacentZTU_Simple, BuildingElementAdjacentZTC, \
    BuildingElementGround, BuildingElementOpaque, BuildingElementTransparent

# Convective fractions
# (default values from BS EN ISO 52016-1:2017, Table B.11)
f_int_c: cython.double = 0.4 # Can be different for each source of internal gains
f_sol_c: cython.double = 0.1

# Areal thermal capacity of air and furniture
# (default value from BS EN ISO 52016-1:2017, Table B.17)
k_m_int: cython.double = 10000 # J / (m2.K)


@cython.cclass
class Zone:
    """ An object to represent a thermal zone in the building

    Temperatures of nodes in associated building elements are also calculated
    in this object by solving a matrix equation.
    """
    __useful_area: cython.double
    __volume: cython.double
    __building_elements: list
    __vent_elements: list
    __vent_cool_extra: object
    __tb_heat_trans_coeff: cython.double
    __area_el_total: cython.double
    __c_int: cython.double
    __element_positions: object
    __zone_idx: cython.int
    __no_of_temps: cython.int
    __temp_prev: cython.double[:]
    __print_heat_balance: cython.bint
    __use_fast_solver  : cython.bint
    __control_obj : object

    def __init__(
            self,
            area: cython.double,
            volume: cython.double,
            building_elements: list,
            thermal_bridging: list,
            vent_obj : object,
            temp_ext_air_init: cython.double,
            temp_setpnt_init: cython.double,
            control_obj : object = None,
            print_heat_balance: cython.bint=False,
            use_fast_solver: cython.bint=False,
            ):
        """ Construct a Zone object

        Arguments:
        area              -- useful floor area of the zone, in m2
        volume            -- total volume of the zone, m3
        building_elements -- list of BuildingElement objects (walls, floors, windows etc.)
        thermal_bridging  -- Either:
                             - overall heat transfer coefficient for thermal
                               bridges in the zone, in W / K
                             - list of ThermalBridge objects for this zone
        vent_elements     -- list of ventilation elements (infiltration, mech vent etc.)
        vent_cool_extra   -- element providing additional ventilation in response to high
                             internal temperature
        temp_ext_air_init -- external air temperature to use during initialisation, in Celsius
        temp_setpnt_init -- setpoint temperature to use during initialisation, in Celsius
        print_heat_balance-- flag to indicate whether to print the heat balance breakdown
        use_fast_solver   -- flag to indicate whether to use the optimised solver (results
                             may differ slightly due to reordering of floating-point ops)
        vent_obj          -- object reference for ventilation class
        control_obj       -- object reference for SetpointTimeControl class
         
        Other variables:
        area_el_total     -- total area of all building elements associated
                             with this zone, in m2
        c_int             -- internal thermal capacity of the zone, in J / K
        element_positions -- dictionary where key is building element and
                             values are 2-element tuples storing matrix row and
                             column numbers (both same) where the first element
                             of the tuple gives the position of the heat
                             balance eqn (row) and node temperature (column)
                             for the external surface and the second element
                             gives the position for the internal surface.
                             Positions in between will be for the heat balance
                             and temperature of the inside nodes of the
                             building element
        zone_idx          -- matrix row and column number (both same)
                             corresponding to heat balance eqn for zone (row)
                             and temperature of internal air (column)
        no_of_temps       -- number of unknown temperatures (each node in each
                             building element + 1 for internal air) to be
                             solved for
        temp_prev         -- list of temperatures (nodes and internal air) from
                             previous timestep. Positions in list defined in
                             self.__element_positions and self.__zone_idx
        """
        self.__useful_area       = area
        self.__volume            = volume
        self.__building_elements = building_elements
        self.__vent_obj          = vent_obj
        self.__control_obj       = control_obj

        # If thermal_bridging is a list of ThermalBridge objects, calculate the
        # overall heat transfer coefficient for thermal bridges, otherwise just
        # use the coefficient given
        if isinstance(thermal_bridging, list):
            self.__tb_heat_trans_coeff = 0.0
            for tb in thermal_bridging:
                self.__tb_heat_trans_coeff += tb.heat_trans_coeff()
        else:
            self.__tb_heat_trans_coeff = thermal_bridging

        self.__area_el_total = sum([eli.area for eli in self.__building_elements])
        self.__c_int = k_m_int * area

        # Calculate:
        # - size of required matrix/vectors (total number of nodes across all
        #   building elements + 1 for internal air)
        # - positions of heat balance eqns and temperatures in matrix for each node
        self.__element_positions = {}
        n: cython.Py_ssize_t = 0
        start_idx: cython.Py_ssize_t
        end_idx: cython.Py_ssize_t
        eli: object
        for eli in self.__building_elements:
            start_idx = n
            n = n + eli.no_of_nodes()
            end_idx = n - 1
            self.__element_positions[eli] = (start_idx, end_idx)
        self.__zone_idx = n
        self.__no_of_temps = n + 1

        self.__print_heat_balance = print_heat_balance
        self.__use_fast_solver = use_fast_solver

        self.__init_node_temps(temp_ext_air_init, temp_setpnt_init)

    @cython.cfunc
    def __init_node_temps(self, temp_ext_air_init: cython.double, temp_setpnt_init: cython.double) -> cython.void:
        """ Initialise temperatures of heat balance nodes

        Arguments:
        temp_ext_air_init -- external air temperature to use during initialisation, in Celsius
        temp_setpnt_init -- setpoint temperature to use during initialisation, in Celsius
        """
        # Use yearly timestep for warm-up period
        # - solution converges significantly faster with larger timestep
        # - solution is the same to ~5 significant figures for hourly vs. yearly timestep
        delta_t_h: cython.double = 8760
        delta_t = delta_t_h * units.seconds_per_hour

        # Assume default convective fraction for heating/cooling suggested in
        # BS EN ISO 52016-1:2017 Table B.11
        frac_convective: cython.double = 0.4

        # Set starting point for all node temperatures (elements of
        # self.__temp_prev) as average of external air temp and setpoint. This
        # is somewhat arbitrary, but of all options for a uniform initial
        # temperature, this should lead to relatively fast stabilisation of
        # fabric temperatures, which are expected to be close to the external
        # air temperature towards the external surface nodes and close to the
        # setpoint temperature towards the internal surface nodes, and therefore
        # node temperatures on average should be close to the average of the
        # external air and setpoint temperatures.
        temp_start: cython.double = (temp_ext_air_init + temp_setpnt_init) / 2.0
        self.__temp_prev = np.array([temp_start] * self.__no_of_temps)

        # Iterate over space heating calculation and meet all space heating
        # demand until temperatures stabilise, under steady-state conditions
        # using specified constant setpoint and external air temperatures.
        space_heat_demand: cython.double
        space_cool_demand: cython.double
        while True:
            space_heat_demand, space_cool_demand, _, _ = self.space_heat_cool_demand(
                delta_t_h,
                temp_ext_air_init,
                0.0, # Internal gains
                0.0, # Solar gains
                frac_convective,
                frac_convective,
                temp_setpnt_init,
                temp_setpnt_init,
                avg_air_supply_temp=temp_ext_air_init,
                ach_cooling=0.0,
                )

            # Note: space_cool_demand returned by function above is negative,
            # and only one of space_heat_demand and space_cool_demand will be
            # non-zero.
            gains_heat_cool: cython.double = (space_heat_demand + space_cool_demand) * units.W_per_kW / delta_t_h

            temps_updated: cython.double[:]

            temps_updated, _ = self.__calc_temperatures(
                delta_t,
                self.__temp_prev,
                temp_ext_air_init,
                0.0, # Internal gains
                0.0, # Solar gains
                gains_heat_cool,
                frac_convective,
                0.0, # ach
                temp_ext_air_init, # avg_supply_temp
                )

            if not np.isclose(temps_updated, self.__temp_prev, rtol=1e-08).all():
                self.__temp_prev = temps_updated
            else:
                break

    def area(self) -> cython.double:
        return self.__useful_area

    def volume(self) -> cython.double:
        return self.__volume

    def gains_solar(self) -> cython.double:
        """sum solar gains for all elements in the zone
        only transparent elements will have solar gains > 0 """

        solar_gains: cython.double = 0
        eli: object
        for eli in self.__building_elements:
            solar_gains += eli.solar_gains()

        return solar_gains

    def __calc_vent_heat_transfer_coeff(self, air_changes_per_hour):
        """ Calculate ventilation heat transfer co-efficient from air changes per hour """
        q_ve = air_changes_per_hour * self.__volume / units.seconds_per_hour
        h_ve = AIR.density_kg_per_m3() * AIR.specific_heat_capacity() * q_ve
        return h_ve

    def __calc_temperatures(self,
            delta_t: cython.double,
            temp_prev: cython.double[:],
            temp_ext_air: cython.double,
            gains_internal: cython.double,
            gains_solar: cython.double,
            gains_heat_cool: cython.double,
            f_hc_c: cython.double,
            ach: cython.double,
            avg_supply_temp : cython.double,
            print_heat_balance: cython.bint = False,
            ):
        """ Calculate temperatures according to procedure in BS EN ISO 52016-1:2017, section 6.5.6

        Arguments:
        delta_t         -- calculation timestep, in seconds
        temp_prev       -- temperature vector X (see below) from previous timestep
        temp_ext_air    -- temperature of external air, in deg C
        gains_internal  -- total internal heat gains, in W
        gains_solar     -- directly transmitted solar gains, in W
        gains_heat_cool -- gains from heating (positive) or cooling (negative), in W
        f_hc_c          -- convective fraction for heating/cooling
        ach             -- air changes per hour
        print_heat_balance -- flag to record whether to return the heat balance outputs
        avg_supply_temp -- average supply temperature

        Temperatures are calculated by solving (for X) a matrix equation A.X = B, where:
        A is a matrix of known coefficients
        X is a vector of unknown temperatures
        B is a vector of known quantities

        Each row in vector X is a temperature variable - one for each node in each
        building element plus the internal air temperature in the zone.

        Each row of matrix A contains the coefficients from the heat balance equations
        for each of the nodes in each building element, plus one row for the heat
        balance equation of the zone.

        Each column of matrix A contains the coefficients for a particular temperature
        variable (in same order that they appear in vector X). Where the particular
        temperature does not appear in the equation this coefficient will be zero.

        Note that for this implementation, the columns and rows will be in corresponding
        order, so the heat balance equation for node i will be in row i and the
        coefficients in each row for the temperature at node i will be in column i.

        Each row of vector B contains the other quantities (i.e. those that are not
        coefficients of the temperature variables) from the heat balance equations
        for each of the nodes in each building element, plus one row for the heat
        balance equation of the zone, in the same order that the rows appear in matrix
        A.
        """
        h_ve = self.__calc_vent_heat_transfer_coeff(ach)

        # Init matrix with zeroes
        # Number of rows in matrix = number of columns
        # = total number of nodes + 1 for overall zone heat balance (and internal air temp)
        matrix_a: cython.double[:, :] = np.zeros((self.__no_of_temps, self.__no_of_temps))

        # Init vector_b with zeroes (length = number of nodes + 1 for overall zone heat balance)
        vector_b: cython.double[:] = np.zeros(self.__no_of_temps)

        # One term in eqn 39 is sum from k = 1 to n of (A_elk / A_tot). Given
        # that A_tot is defined as the sum of A_elk from k = 1 to n, this term
        # will always evaluate to 1.
        # TODO Check this is correct. It seems a bit pointless if it is but we
        #      should probably retain it as an explicit term anyway to match
        #      the standard.
        sum_area_frac: cython.double = 1.0

        # Node heat balances - loop through building elements and their nodes:
        # - Construct row of matrix_a for each node energy balance eqn
        # - Calculate RHS of node energy balance eqn and add to vector_b
        idx: cython.Py_ssize_t
        i: cython.Py_ssize_t
        h_ci: cython.double
        eli: object
        elk: object
        i_sol_dir: cython.double
        i_sol_dif: cython.double
        f_sh_dir: cython.double
        f_sh_dif: cython.double

        for eli in self.__building_elements:
            # External surface node (eqn 41)
            # Get position (row == column) in matrix previously calculated for the first (external) node
            idx = self.__element_positions[eli][0]
            # Position of first (external) node within element is zero
            i = 0
            # Coeff for temperature of this node
            matrix_a[idx][idx] = (eli.k_pli[i] / delta_t) + eli.h_ce() + eli.h_re() + eli.h_pli[i]
            # Coeff for temperature of next node
            matrix_a[idx][idx + 1] = - eli.h_pli[i]
            # RHS of heat balance eqn for this node
            i_sol_dir, i_sol_dif = eli.i_sol_dir_dif()
            f_sh_dir, f_sh_dif = eli.shading_factors_direct_diffuse()
            vector_b[idx] = (eli.k_pli[i] / delta_t) * temp_prev[idx] \
                          + (eli.h_ce() + eli.h_re()) * eli.temp_ext() \
                          + eli.a_sol * (i_sol_dif * f_sh_dif + i_sol_dir * f_sh_dir) \
                          - eli.therm_rad_to_sky

            # Inside node(s), if any (eqn 40)
            for i in range(1, eli.no_of_inside_nodes() + 1):
                idx = idx + 1
                # Coeff for temperature of prev node
                matrix_a[idx][idx - 1] = - eli.h_pli[i - 1]
                # Coeff for temperature of this node
                matrix_a[idx][idx] = (eli.k_pli[i] / delta_t) + eli.h_pli[i] + eli.h_pli[i - 1]
                # Coeff for temperature of next node
                matrix_a[idx][idx + 1] = - eli.h_pli[i]
                # RHS of heat balance eqn for this node
                vector_b[idx] = (eli.k_pli[i] / delta_t) * temp_prev[idx]

            # Internal surface node (eqn 39)
            idx = idx + 1
            assert idx == self.__element_positions[eli][1]
            i = i + 1
            assert i == eli.no_of_nodes() - 1
            # Get internal convective surface heat transfer coefficient, which
            # depends on direction of heat flow, which depends in temperature of
            # zone and internal surface
            h_ci = eli.h_ci(temp_prev[self.__zone_idx], temp_prev[idx])
            # Coeff for temperature of prev node
            matrix_a[idx][idx - 1] = - eli.h_pli[i - 1]
            # Coeff for temperature of this node
            matrix_a[idx][idx] = (eli.k_pli[i] / delta_t) + h_ci \
                               + eli.h_ri() * sum_area_frac + eli.h_pli[i - 1]
            # Add final sum term for LHS of eqn 39 in loop below.
            # These are coeffs for temperatures of internal surface nodes of
            # all building elements in the zone
            col: cython.Py_ssize_t
            for elk in self.__building_elements:
                col = self.__element_positions[elk][1]
                # The line below must be an adjustment to the existing value
                # to handle the case where col = idx (i.e. where we have
                # already partially set the value of the matrix element above
                # (before this loop) and do not want to overwrite it)
                matrix_a[idx][col] = matrix_a[idx][col] \
                                   - (elk.area / self.__area_el_total) * eli.h_ri()
            # Coeff for temperature of thermal zone
            matrix_a[idx][self.__zone_idx] = - h_ci
            # RHS of heat balance eqn for this node
            vector_b[idx] = (eli.k_pli[i] / delta_t) * temp_prev[idx] \
                          + ( (1.0 - f_int_c) * gains_internal \
                            + (1.0 - f_sol_c) * gains_solar \
                            + (1.0 - f_hc_c) * gains_heat_cool \
                            ) \
                          / self.__area_el_total

        # Zone heat balance:
        # - Construct row of matrix A for zone heat balance eqn
        # - Calculate RHS of zone heat balance eqn and add to vector_b

        # Coeff for temperature of thermal zone
        # TODO Throughput factor only applies to MVHR and WHEV, therefore only
        #      these systems accept throughput_factor as an argument to the h_ve
        #      function, hence the branch on the type in the loop below. This
        #      means that the MVHR and WHEV classes no longer have the same
        #      interface as other ventilation element classes, which could make
        #      future development more difficult. Ideally, we would find a
        #      cleaner way to implement this difference.

        matrix_a[self.__zone_idx][self.__zone_idx] \
            = (self.__c_int / delta_t) \
            + sum([ eli.area
                  * eli.h_ci(
                      temp_prev[self.__zone_idx],
                      temp_prev[self.__element_positions[eli][1]]
                      )
                  for eli in self.__building_elements
                  ]) \
            + h_ve \
            + self.__tb_heat_trans_coeff
        # Add final sum term for LHS of eqn 38 in loop below.
        # These are coeffs for temperatures of internal surface nodes of
        # all building elements in the zone
        for eli in self.__building_elements:
            col = self.__element_positions[eli][1] # Column for internal surface node temperature
            matrix_a[self.__zone_idx][col] \
                = - eli.area \
                * eli.h_ci(temp_prev[self.__zone_idx], temp_prev[self.__element_positions[eli][1]])
        # RHS of heat balance eqn for zone
        vector_b[self.__zone_idx] \
            = (self.__c_int / delta_t) * temp_prev[self.__zone_idx] \
            + h_ve * avg_supply_temp \
            + self.__tb_heat_trans_coeff * temp_ext_air \
            + f_int_c * gains_internal \
            + f_sol_c * gains_solar \
            + f_hc_c * gains_heat_cool

        # Solve matrix eqn A.X = B to calculate vector_x (temperatures)
        if self.__use_fast_solver:
            vector_x = self.__fast_solver(matrix_a, vector_b)
        else:
            vector_x = np.linalg.solve(matrix_a, vector_b)

        heat_balance_dict: dict
        temp_internal: cython.double
        hb_gains_solar: cython.double
        hb_gains_internal: cython.double
        hb_gains_heat_cool: cython.double
        hb_energy_to_change_temp: cython.double
        hb_loss_thermal_bridges: cython.double
        hb_loss_infiltration_ventilation: cython.double
        hb_loss_fabric: cython.double
        fabric_int_sol: cython.double
        fabric_int_int_gains: cython.double
        fabric_int_heat_cool: cython.double
        fabric_int_air_convective: cython.double
        temp_int_surface: cython.double
        air_node_temp: cython.double
        hb_fabric_ext_air_convective: cython.double
        hb_fabric_ext_air_radiative: cython.double
        hb_fabric_ext_sol: cython.double
        hb_fabric_ext_sky: cython.double
        hb_fabric_ext_opaque: cython.double
        hb_fabric_ext_transparent: cython.double
        hb_fabric_ext_ground: cython.double
        hb_fabric_ext_ZTC: cython.double
        hb_fabric_ext_ZTU: cython.double

        if print_heat_balance:
            heat_balance_dict = {}

            # Collect outputs, in W, for heat balance at air node
            temp_internal = vector_x[self.__zone_idx]
            hb_gains_solar = f_sol_c * gains_solar
            hb_gains_internal = f_int_c * gains_internal
            hb_gains_heat_cool = f_hc_c * gains_heat_cool
            hb_energy_to_change_temp = -(self.__c_int / delta_t)*(temp_internal-temp_prev[self.__zone_idx])
            hb_loss_thermal_bridges = self.__tb_heat_trans_coeff*(temp_internal-temp_ext_air)
            hb_loss_infiltration_ventilation = h_ve * (temp_internal- avg_supply_temp)
            hb_loss_fabric = (hb_gains_solar+hb_gains_internal+hb_gains_heat_cool+hb_energy_to_change_temp)\
                            -(hb_loss_thermal_bridges+hb_loss_infiltration_ventilation)
            heat_balance_dict['air_node'] = {
                'solar gains' : hb_gains_solar,
                'internal gains' : hb_gains_internal,
                'heating or cooling system gains' : hb_gains_heat_cool,
                'energy to change internal temperature' : hb_energy_to_change_temp,
                'thermal_bridges' : - hb_loss_thermal_bridges,
                'infiltration_ventilation' : - hb_loss_infiltration_ventilation,
                'fabric' : - hb_loss_fabric
                                }

            # Collect outputs, in W, for heat balance at internal fabric boundary
            fabric_int_sol = (1.0 - f_sol_c) * gains_solar
            fabric_int_int_gains = (1.0 - f_int_c) * gains_internal
            fabric_int_heat_cool = (1.0 - f_hc_c) * gains_heat_cool
            
            fabric_int_air_convective = 0.0
            
            for eli in self.__building_elements:
                idx = self.__element_positions[eli][1]
                temp_int_surface = vector_x[idx]
                air_node_temp = vector_x[self.__zone_idx]
                
                fabric_int_air_convective += eli.area \
                     * ( (eli.h_ci(temp_prev[self.__zone_idx], temp_prev[self.__element_positions[eli][1]])) * (air_node_temp - temp_int_surface))
               
            heat_balance_dict['internal_boundary'] = {
                'fabric_int_air_convective': fabric_int_air_convective,
                'fabric_int_sol': fabric_int_sol,
                'fabric_int_int_gains': fabric_int_int_gains,
                'fabric_int_heat_cool':fabric_int_heat_cool,
                }

            # Collect outputs, in W, for heat balance at external boundary
            hb_fabric_ext_air_convective = 0.0
            hb_fabric_ext_air_radiative = 0.0
            hb_fabric_ext_sol = 0.0
            hb_fabric_ext_sky = 0.0
            hb_fabric_ext_opaque = 0.0
            hb_fabric_ext_transparent = 0.0
            hb_fabric_ext_ground = 0.0
            hb_fabric_ext_ZTC = 0.0
            hb_fabric_ext_ZTU = 0.0
            
            for eli in self.__building_elements:
                # Get position in vector for the first (external) node of the building element
                idx = self.__element_positions[eli][0]
                temp_ext_surface = vector_x[idx]
                i_sol_dir, i_sol_dif = eli.i_sol_dir_dif()
                f_sh_dir, f_sh_dif = eli.shading_factors_direct_diffuse()
                hb_fabric_ext_air_convective += eli.area \
                     * ( (eli.h_ce()) * (eli.temp_ext() - temp_ext_surface))
                hb_fabric_ext_air_radiative += eli.area \
                     * (eli.h_re()) * (eli.temp_ext() - temp_ext_surface)
                hb_fabric_ext_sol += eli.area \
                    * eli.a_sol * (i_sol_dif * f_sh_dif + i_sol_dir * f_sh_dir)
                hb_fabric_ext_sky += eli.area * (- eli.therm_rad_to_sky)
                #fabric heat loss per building element type
                hb_fabric_ext = eli.area \
                     * ( (eli.h_ce()) * (eli.temp_ext() - temp_ext_surface)) \
                     + eli.area \
                     * (eli.h_re()) * (eli.temp_ext() - temp_ext_surface) \
                     + eli.area \
                     * eli.a_sol * (i_sol_dif * f_sh_dif + i_sol_dir * f_sh_dir) \
                     + eli.area * (- eli.therm_rad_to_sky)
                if type(eli) in (BuildingElementOpaque,):
                    hb_fabric_ext_opaque += hb_fabric_ext
                elif type(eli) in (BuildingElementTransparent,):
                    hb_fabric_ext_transparent += hb_fabric_ext
                elif type(eli) in (BuildingElementGround,):
                    hb_fabric_ext_ground += hb_fabric_ext
                elif type(eli) in (BuildingElementAdjacentZTC,):
                    hb_fabric_ext_ZTC += hb_fabric_ext
                elif type(eli) in (BuildingElementAdjacentZTU_Simple,):
                    hb_fabric_ext_ZTU += hb_fabric_ext
            heat_balance_dict['external_boundary'] = {
                'solar gains': gains_solar,
                'internal gains': gains_internal,
                'heating or cooling system gains': gains_heat_cool,
                'thermal_bridges': - hb_loss_thermal_bridges,
                'infiltration_ventilation': - hb_loss_infiltration_ventilation,
                'fabric_ext_air_convective': hb_fabric_ext_air_convective,
                'fabric_ext_air_radiative': hb_fabric_ext_air_radiative,
                'fabric_ext_sol': hb_fabric_ext_sol,
                'fabric_ext_sky': hb_fabric_ext_sky,
                'opaque_fabric_ext': hb_fabric_ext_opaque,
                'transparent_fabric_ext': hb_fabric_ext_transparent,
                'ground_fabric_ext': hb_fabric_ext_ground,
                'ZTC_fabric_ext': hb_fabric_ext_ZTC,
                'ZTU_fabric_ext': hb_fabric_ext_ZTU,
                }
        else:
            heat_balance_dict = None
        return vector_x, heat_balance_dict

    @cython.cfunc
    def __fast_solver(self, coeffs: cython.double[:, :], rhs: cython.double[:]) -> cython.double[:]:
        """ Optimised heat balance solver

        Arguments:
        coeffs -- full matrix of coefficients for the heat balance eqns
        rhs -- full vector of values that are not temperatures or coefficients
               (i.e. terms on right hand side of heat balance eqns)

        The heat balance equations from BS EN ISO 52016-1:2017 are expressed as a matrix equation and
        solved simultaneously. While this provides a generic calculation procedure that works for an
        arbitrary number of nodes (N), it also has a runtime proportional to N^3 which means that more
        complex buildings can take a long time to simulate. However, many of the nodes are known not to
        interact (e.g. the node in the middle of one wall has no heat transfer with the node in another
        wall) and therefore we do not require the full flexibility of the matrix approach to solve for
        every node temperature. The only part of the heat balance calculation where this flexibility is
        needed is in the interaction between internal air and internal surfaces, so the calculation of the
        other node temperatures can be removed from the matrix equation using algebraic substitution.

        Consider generic heat balance eqns for a 4-node element:

        A1a + B1b                               = Z1    # Heat balance at node a (external surface node)
        A2a + B2b + C2c                         = Z2    # Heat balance at node b (inside node)
              B3b + C3c + D3d                   = Z3    # Heat balance at node c (inside node)
                    C4c + D4d + J4j + K4k + Y4y = Z4    # Heat balance at node d (internal surface node)
        where:
        - a, b, c and d are the node temperatures in the building element to be solved for
        - j and k are the node temperatures for the internal surfaces of other elements
        - y is the internal air temperature
        - A1 is the coefficient for temperature a in equation 1, A2 is the coefficient for temperature a in
          equation 2, etc.
        - Z1, Z2, etc. are the terms in equation 1, 2, etc. that are not the temperatures to be solved for
          or their coefficients (i.e. Z1, Z2 etc. are the terms on the RHS of the equations)

        The heat balance equation for node a (external surface node) can be rearranged to solve for a:

        A1a + B1b = Z1
        A1a = Z1 - B1b
        a = (Z1 - B1b) / A1

        Using the rearranged heat balance equation for node a, we can substitute a in the heat balance
        equation for node b (next inside node) to eliminate a as a variable:

        A2a + B2b + C2c = Z2
        A2 * (Z1 - B1b) / A1 + B2b + C2c = Z2

        Rearranging to consolidate the occurances of b gives a new heat balance equation for b:

        A2 * Z1 / A1 + A2 * (- B1b) / A1 + B2b + C2c = Z2
        b * (B2 - A2 * B1 / A1) + A2 * Z1 / A1 + C2c = Z2
        (B2 - A2 * B1 / A1) * b + C2c = Z2 - A2 * Z1 / A1

        This new heat balance equation can then be expressed in terms of modified versions of B2 and Z2:

        B2'b + C2c = Z2'
        where:
        B2' = B2 - B1 * A2 / A1
        Z2' = Z2 - Z1 * A2 / A1

        The process can then be repeated, rearranging this new heat balance equation to solve for b and
        then substituting into the heat balance equation for c:

        b = (Z2' - C2c) / B2'

        C3'c + D3d = Z3'
        where:
        C3' = C3 - C2 * B3 / B2'
        Z3' = Z3 - Z2' * B3 / B2'

        And repeated again to generate a new equation for node d:

        c = (Z3' - D3d) / C3'

        D4'd + J4j + K4k + Y4y = Z4'
        where:
        D4' = D4 - D3 * C4 / C3'
        Z4' = Z4 - Z3' * C4 / C3'

        At this point, we have reached the internal surface and we need the flexibility of the matrix
        approach, but we have reduced the number of nodes to be solved for by 3. The process can be
        repeated for each building element, and once the matrix solver has solved for the internal surface
        temperatures, we can then go back through the other nodes, using the rearranged heat balance
        equations that solve for (in this case) c, b and a.

        In order to deal with different building elements having different numbers of nodes, we can express
        the above relationships generically:

        temperature[i] = (Z_adjusted[i] - coeff[i][i+1] * temperature[i+1]) / coeff_adjusted[i][i]

        coeff_adjusted[i][i] = coeff[i][i] - coeff[i-1][i] * coeff[i][i-1] / coeff_adjusted[i-1][i-1]
        Z_adjusted[i] = Z[i] - Z_adjusted[i-1] * coeff[i][i-1] / coeff_adjusted[i-1][i-1]

        where i is the number of the node and its heat balance equation (counting from the external surface
        to the internal surface), e.g. coeff[i-1][i] would be the coeffient for the temperature of the
        current node in the heat balance equation for the previous node.

        The optimised calculation procedure is therefore:
        - Loop over nodes, from external surface to internal surface, and calculate adjusted coeffs and RHS
          for each heat balance eqn
        - Construct matrix eqn for inside and air nodes only
        - Solve heat balance eqns for inside and air nodes using normal matrix solver
        - Loop over nodes, from internal inside node (i.e. inside node nearest to the internal surface) to
          external surface, and calculate temperatures in sequence
        """
        # Init matrix with zeroes
        # Number of rows in matrix = number of columns
        # = total number of nodes + 1 for overall zone heat balance (and internal air temp)
        coeffs_adj: cython.double[:, :] = np.zeros((self.__no_of_temps, self.__no_of_temps))

        coeffs_adj[0]

        # Init vector_b with zeroes (length = number of nodes + 1 for overall zone heat balance)
        rhs_adj: cython.double[:] = np.zeros(self.__no_of_temps)

        # Init matrix with zeroes
        # Number of rows in matrix = number of columns
        # = total number of internal surface nodes + 1 for internal air node
        num_rows_cols_optimised: cython.int = len(self.__building_elements) + 1
        zone_idx: cython.Py_ssize_t = num_rows_cols_optimised - 1
        matrix_a: cython.double[:, :] = np.zeros((num_rows_cols_optimised, num_rows_cols_optimised))

        # Init vector_b with zeroes (length = number of internal surfaces + 1 for air node)
        vector_b: cython.double[:] = np.zeros(num_rows_cols_optimised)

        # Loop over building elements
        idx_ext_surface: cython.Py_ssize_t
        idx_int_surface: cython.Py_ssize_t
        el_idx: cython.Py_ssize_t
        eli: object
        idx: cython.Py_ssize_t
        for el_idx, eli in enumerate(self.__building_elements):
            idx_ext_surface = self.__element_positions[eli][0]
            idx_int_surface = self.__element_positions[eli][1]

            # No adjusted coeffs and RHS for external surface heat balance eqn
            coeffs_adj[idx_ext_surface][idx_ext_surface] = coeffs[idx_ext_surface][idx_ext_surface]
            rhs_adj[idx_ext_surface] = rhs[idx_ext_surface]

            # Loop over nodes, from inside node adjacent to external surface, to internal surface
            for idx in range(idx_ext_surface + 1, idx_int_surface + 1):
                # Calculate adjusted coeffs and RHS for each heat balance eqn
                coeffs_adj[idx][idx] \
                    = coeffs[idx][idx] \
                    - coeffs[idx - 1][idx] * coeffs[idx][idx - 1] / coeffs_adj[idx - 1][idx - 1]
                rhs_adj[idx] \
                    = rhs[idx] \
                    - rhs_adj[idx - 1] * coeffs[idx][idx - 1] / coeffs_adj[idx - 1][idx - 1]

            # Construct matrix eqn for internal surface nodes only (and air node, after this loop)
            matrix_a[el_idx][el_idx] = coeffs_adj[idx_int_surface][idx_int_surface]
            vector_b[el_idx] = rhs_adj[idx_int_surface]

            # Add coeffs for temperatures other than the int surface temp of this building element
            el_idx_other: cython.Py_ssize_t
            elk: object
            for el_idx_other, elk in enumerate(self.__building_elements):
                # Skip the current building element
                if elk is eli:
                    continue

                idx_other_int_surface = self.__element_positions[elk][1]
                matrix_a[el_idx][el_idx_other] = coeffs[idx_int_surface][idx_other_int_surface]

            # Add coeff for air temperature to this element's internal surface heat balance eqn
            matrix_a[el_idx][zone_idx] = coeffs[idx_int_surface][self.__zone_idx]
            # Add coeff for this element's internal surface temp to the air node heat balance eqn
            matrix_a[zone_idx][el_idx] = coeffs[self.__zone_idx][idx_int_surface]

        # Add rest of air node heat balance eqn to matrix
        # Coeffs for temperatures other than the air temp are added in the loop above
        matrix_a[zone_idx][zone_idx] = coeffs[self.__zone_idx][self.__zone_idx]
        vector_b[zone_idx] = rhs[self.__zone_idx]

        # Solve heat balance eqns for inside and air nodes using normal matrix solver
        # Solve matrix eqn A.X = B to calculate vector_x (temperatures)
        vector_x: cython.double[:] = np.linalg.solve(matrix_a, vector_b)

        # Init vector_x with zeroes (length = number of nodes + 1 for overall zone heat balance)
        temperatures: cython.double[:] = np.zeros(self.__no_of_temps)

        # Populate air node temperature result
        temperatures[self.__zone_idx] = vector_x[zone_idx]

        # Populate node temperature results for each building element
        for el_idx, eli in enumerate(self.__building_elements):
            idx_ext_surface = self.__element_positions[eli][0]
            idx_int_surface = self.__element_positions[eli][1]

            # Populate internal surface temperature result
            temperatures[idx_int_surface] = vector_x[el_idx]

            # Loop over nodes, from internal inside node (i.e. inside node nearest to the
            # internal surface) to external surface, and calculate temperatures in sequence
            for idx in reversed(range(idx_ext_surface, idx_int_surface)):
                temperatures[idx] \
                    = (rhs_adj[idx] - coeffs[idx][idx + 1] * temperatures[idx + 1]) \
                    / coeffs_adj[idx][idx]

        return temperatures

    @cython.cfunc
    def __temp_operative(self, temp_vector: cython.double[:]) -> cython.double:
        """ Calculate the operative temperature, in deg C

        According to the procedure in BS EN ISO 52016-1:2017, section 6.5.5.3.

        Arguments:
        temp_vector -- vector (list) of temperatures calculated from the heat balance equations
        """
        temp_int_air: cython.double = temp_vector[self.__zone_idx]

        # Mean radiant temperature is weighted average of internal surface temperatures
        temp_mean_radiant: cython.double = sum([eli.area * temp_vector[self.__element_positions[eli][1]] \
                                 for eli in self.__building_elements]) \
                          / self.__area_el_total

        return (temp_int_air + temp_mean_radiant) / 2.0

    def temp_operative(self) -> cython.double:
        """ Return operative temperature, in deg C """
        return self.__temp_operative(self.__temp_prev)

    def temp_internal_air(self) -> cython.double:
        """ Return internal air temperature, in deg C """
        return self.__temp_prev[self.__zone_idx]

    def ___ach_req_to_reach_temperature(
            self,
            temp_target: cython.double,
            ach_min: cython.double,
            ach_max: cython.double,
            temp_ach_min: cython.double,
            temp_ach_max: cython.double,
            temp_int_air_min: cython.double,
            temp_int_air_max: cython.double,
            temp_supply: cython.double,
            ):
        """ Calculate the air change rate required to meet a target temperature """
        if temp_ach_max >= temp_ach_min or temp_int_air_max >= temp_int_air_min \
        or temp_int_air_min <= temp_supply:
            return ach_min

        frac_interp = (temp_target - temp_ach_min) / (temp_ach_max - temp_ach_min)
        temp_int_air_req = temp_int_air_min + frac_interp * (temp_int_air_max - temp_int_air_min)

        if temp_int_air_req <= temp_supply:
            return ach_max

        ach_req \
            = ( ach_max * frac_interp \
              * ( (temp_int_air_max - temp_supply)
                / (temp_int_air_min - temp_supply)
                ) \
              + ach_min * (1.0 - frac_interp) \
              ) \
            / ( (temp_int_air_req - temp_supply)
              / (temp_int_air_min - temp_supply)
              )
        return min(max(ach_req, ach_min), ach_max)

    def __calc_cooling_potential_from_ventilation(
            self,
            delta_t,
            temp_ext_air,
            gains_internal,
            gains_solar,
            gains_heat_cool,
            frac_conv_gains_heat_cool,
            temp_setpnt_heat,
            temp_setpnt_cool,
            temp_setpnt_cool_vent,
            temp_operative_free,
            temp_int_air_free,
            ach_cooling,
            ach_windows_open,
            ach_target,
            avg_supply_temp,
            ):
        # If ach required for cooling has not been provided, check if the
        # maximum air changes when window shut and open are different
        if ach_cooling is None and ach_windows_open != ach_target:
            #Calculate node and internal air temperatures with maximum additional ventilation
            temp_vector_vent_max, _ = self.__calc_temperatures(
                delta_t,
                self.__temp_prev,
                temp_ext_air,
                gains_internal,
                gains_solar,
                gains_heat_cool,
                frac_conv_gains_heat_cool,
                ach_windows_open,
                avg_supply_temp,
                )

            # Calculate internal operative temperature with maximum ventilation
            temp_operative_vent_max = self.__temp_operative(temp_vector_vent_max)
            temp_int_air_vent_max = temp_vector_vent_max[self.__zone_idx]

            ach_to_trigger_heating = self.___ach_req_to_reach_temperature(
                temp_setpnt_heat,
                ach_target,
                ach_windows_open,
                temp_operative_free,
                temp_operative_vent_max,
                temp_int_air_free,
                temp_int_air_vent_max,
                avg_supply_temp,
                )

            # If there is cooling potential from additional ventilation, and
            # free-floating temperature exceeds setpoint for additional ventilation
            if temp_operative_vent_max < temp_operative_free \
            and temp_operative_free > temp_setpnt_cool_vent \
            and temp_int_air_free > avg_supply_temp:
                # Calculate ventilation required to reach cooling setpoint for ventilation
                ach_cooling = self.___ach_req_to_reach_temperature(
                    temp_setpnt_cool_vent,
                    ach_target,
                    ach_windows_open,
                    temp_operative_free,
                    temp_operative_vent_max,
                    temp_int_air_free,
                    temp_int_air_vent_max,
                    avg_supply_temp,
                    )

                # Calculate node and internal air temperatures with heating/cooling gains of zero
                temp_vector_no_heat_cool_vent_extra, _ = self.__calc_temperatures(
                    delta_t,
                    self.__temp_prev,
                    temp_ext_air,
                    gains_internal,
                    gains_solar,
                    gains_heat_cool,
                    frac_conv_gains_heat_cool,
                    ach_cooling,
                    avg_supply_temp,
                    )

                # Calculate internal operative temperature at free-floating conditions
                # i.e. with no heating/cooling
                temp_operative_free_vent_extra = self.__temp_operative(temp_vector_no_heat_cool_vent_extra)

                # If temperature achieved by additional ventilation is above setpoint
                # for active cooling, assume cooling system will be used instead of
                # additional ventilation. Otherwise, use resultant operative temperature
                # in calculation of space heating/cooling demand.
                if temp_operative_free_vent_extra > temp_setpnt_cool:
                    ach_cooling = ach_target
                else:
                    temp_operative_free = temp_operative_free_vent_extra
        else:
            ach_to_trigger_heating = None

        if ach_cooling is None:
            ach_cooling = ach_target

        return temp_operative_free, ach_cooling, ach_to_trigger_heating

    def __interp_heat_cool_demand(
            self,
            delta_t_h,
            temp_setpnt,
            heat_cool_load_upper,
            temp_operative_free,
            temp_operative_upper,
            ):
        if (temp_operative_upper - temp_operative_free) == 0.0:
            sys.exit(
                "ERROR: Divide-by-zero in calculation of heating/cooling demand. "
                "This may be caused by the specification of very low overall "
                "areal heat capacity of BuildingElements and/or very high thermal "
                "mass of WetDistribution."
                )

        heat_cool_load_unrestricted \
            = heat_cool_load_upper \
            * (temp_setpnt - temp_operative_free) \
            / (temp_operative_upper - temp_operative_free)
        # Convert from W to kWh
        heat_cool_demand: cython.double \
            = heat_cool_load_unrestricted / units.W_per_kW * delta_t_h
        return heat_cool_demand

    def space_heat_cool_demand(
            self,
            delta_t_h: cython.double,
            temp_ext_air: cython.double,
            gains_internal: cython.double,
            gains_solar: cython.double,
            frac_convective_heat: cython.double,
            frac_convective_cool: cython.double,
            temp_setpnt_heat: cython.double,
            temp_setpnt_cool: cython.double,
            avg_air_supply_temp : cython.double,
            gains_heat_cool_convective: cython.double=0.0,
            gains_heat_cool_radiative: cython.double=0.0,
            ach_windows_open :cython.double = None,
            ach_target : cython.double = None,
            ach_cooling : cython.double = None,
            ) -> tuple:
        """ Calculate heating and cooling demand in the zone for the current timestep

        According to the procedure in BS EN ISO 52016-1:2017, section 6.5.5.2, steps 1 to 4.

        Arguments:
        delta_t_h -- calculation timestep, in hours
        temp_ext_air -- temperature of the external air for the current timestep, in deg C
        gains_internal -- internal gains for the current timestep, in W
        gains_solar -- directly transmitted solar gains, in W
        frac_convective_heat -- convective fraction for heating
        frac_convective_cool -- convective fraction for cooling
        temp_setpnt_heat -- temperature setpoint for heating, in deg C
        temp_setpnt_cool -- temperature setpoint for cooling, in deg C
        ach_windows_open -- air changes per hour when all windows open
        ach_target -- air changes per hour required for ventilation requirement/target
        ach_cooling -- air changes per hour required to meet cooling requirement/target
        """
        if temp_setpnt_cool < temp_setpnt_heat:
            sys.exit('ERROR: Cooling setpoint is below heating setpoint.')

        if self.__control_obj is not None:
            temp_setpnt_cool_vent_response = self.__control_obj.setpnt()
            if temp_setpnt_cool_vent_response is None:
                # Set cooling setpoint to Planck temperature to ensure no cooling demand
                temp_setpnt_cool_vent_response = units.Kelvin2Celcius(1.4e32)
            if temp_setpnt_cool_vent_response < temp_setpnt_heat:
                sys.exit('ERROR: Setpoint for additional ventilation is below heating setpoint.')
            temp_setpnt_cool_vent = temp_setpnt_cool_vent_response
        else:
            temp_setpnt_cool_vent = units.Kelvin2Celcius(1.4e32)
        if temp_setpnt_cool_vent < temp_setpnt_heat:
            sys.exit('ERROR: Cooling ventilation setpoint is below heating setpoint.')

        # Calculate timestep in seconds
        delta_t: cython.double = delta_t_h * units.seconds_per_hour

        # For calculation of demand, set heating/cooling gains to zero
        gains_heat_cool: cython.double = gains_heat_cool_convective + gains_heat_cool_radiative
        if gains_heat_cool == 0.0:
            frac_conv_gains_heat_cool: cython.double = 0.0
        else:
            frac_conv_gains_heat_cool: cython.double \
                = gains_heat_cool_convective / gains_heat_cool

        # Calculate node and internal air temperatures with heating/cooling gains of zero
        temp_vector_no_heat_cool: cython.double[:]
        temp_vector_no_heat_cool, _ = self.__calc_temperatures(
            delta_t,
            self.__temp_prev,
            temp_ext_air,
            gains_internal,
            gains_solar,
            gains_heat_cool,
            frac_conv_gains_heat_cool,
            ach_target if ach_cooling is None else ach_cooling,
            avg_air_supply_temp,
            )

        # Calculate internal operative temperature at free-floating conditions
        # i.e. with no heating/cooling
        temp_operative_free: cython.double = self.__temp_operative(temp_vector_no_heat_cool)
        temp_int_air_free: cython.double = temp_vector_no_heat_cool[self.__zone_idx]

        temp_operative_free, ach_cooling, ach_to_trigger_heating \
            = self.__calc_cooling_potential_from_ventilation(
                delta_t,
                temp_ext_air,
                gains_internal,
                gains_solar,
                gains_heat_cool,
                frac_conv_gains_heat_cool,
                temp_setpnt_heat,
                temp_setpnt_cool,
                temp_setpnt_cool_vent,
                temp_operative_free,
                temp_int_air_free,
                ach_cooling,
                ach_windows_open,
                ach_target,
                avg_air_supply_temp,
                )

        # Determine relevant setpoint (if neither, then return space heating/cooling demand of zero)
        # Determine maximum heating/cooling
        temp_setpnt: cython.double
        heat_cool_load_upper: cython.double
        frac_convective: cython.double
        if temp_operative_free > temp_setpnt_cool:
            # Cooling
            # TODO Implement eqn 26 "if max power available" case rather than just "otherwise" case?
            #      Could max. power be available at this point for all heating/cooling systems?
            temp_setpnt = temp_setpnt_cool
            heat_cool_load_upper = - 10.0 * self.__useful_area
            frac_convective = frac_convective_cool
        elif temp_operative_free < temp_setpnt_heat:
            # Heating
            # TODO Implement eqn 26 "if max power available" case rather than just "otherwise" case?
            #      Could max. power be available at this point for all heating/cooling systems?
            temp_setpnt = temp_setpnt_heat
            heat_cool_load_upper = 10.0 * self.__useful_area
            frac_convective = frac_convective_heat
        else:
            return 0.0, 0.0, ach_cooling, ach_to_trigger_heating # No heating or cooling load

        # Calculate node and internal air temperatures with maximum heating/cooling
        gains_heat_cool_upper: cython.double = gains_heat_cool + heat_cool_load_upper
        frac_convective_heat_cool_upper: cython.double \
            = ( gains_heat_cool * frac_conv_gains_heat_cool 
              + heat_cool_load_upper * frac_convective
              ) \
            / gains_heat_cool_upper
        temp_vector_upper_heat_cool, _ = self.__calc_temperatures(
            delta_t,
            self.__temp_prev,
            temp_ext_air,
            gains_internal,
            gains_solar,
            gains_heat_cool_upper,
            frac_convective_heat_cool_upper,
            ach_cooling,
            avg_air_supply_temp,
            )

        # Calculate internal operative temperature with maximum heating/cooling
        temp_operative_upper: cython.double = self.__temp_operative(temp_vector_upper_heat_cool)

        # Calculate heating (positive) or cooling (negative) required to reach setpoint
        heat_cool_demand = self.__interp_heat_cool_demand(
            delta_t_h,
            temp_setpnt,
            heat_cool_load_upper,
            temp_operative_free,
            temp_operative_upper,
            )

        space_heat_demand: cython.double = 0.0 # in kWh
        space_cool_demand: cython.double = 0.0 # in kWh
        if heat_cool_demand < 0.0:
            space_cool_demand = heat_cool_demand
        elif heat_cool_demand > 0.0:
            space_heat_demand = heat_cool_demand
        else:
            pass
        return space_heat_demand, space_cool_demand, ach_cooling, ach_to_trigger_heating

    def update_temperatures(self,
            delta_t: cython.double,
            temp_ext_air: cython.double,
            gains_internal: cython.double,
            gains_solar: cython.double,
            gains_heat_cool: cython.double,
            frac_convective: cython.double,
            ach: cython.double,
            avg_supply_temp,
            ) -> object:
        """ Update node and internal air temperatures for calculation of next timestep

        Arguments:
        delta_t         -- calculation timestep, in seconds
        temp_ext_air    -- temperature of external air, in deg C
        gains_internal  -- total internal heat gains, in W
        gains_solar     -- directly transmitted solar gains, in W
        gains_heat_cool -- gains from heating (positive) or cooling (negative), in W
        frac_convective -- convective fraction for heating/cooling (as appropriate)
        ach             -- air changes per hour
        heat_balance_dict -- dictionary to record heat balance outputs
        """

        # Calculate node and internal air temperatures with calculated heating/cooling gains.
        # Save as "previous" temperatures for next timestep
        self.__temp_prev, heat_balance_dict = self.__calc_temperatures(
            delta_t,
            self.__temp_prev,
            temp_ext_air,
            gains_internal,
            gains_solar,
            gains_heat_cool,
            frac_convective,
            ach,
            avg_supply_temp,
            print_heat_balance = self.__print_heat_balance,
            )
        return heat_balance_dict

    def total_fabric_heat_loss(self) -> cython.double:
        """ Return the total fabric heat loss from all
        building elements in a zone, in W / K """
        total_fabric_heat_loss: cython.double = 0
        for be in self.__building_elements:
            total_fabric_heat_loss += be.fabric_heat_loss()
        return total_fabric_heat_loss

    def total_heat_loss_area(self) -> cython.double:
        """ return the total heat loss area, in m2"""
        total_heat_loss_area: cython.double = 0
        for be in self.__building_elements:
            if type(be) in (BuildingElementOpaque, BuildingElementTransparent,
             BuildingElementGround, BuildingElementAdjacentZTU_Simple):
                total_heat_loss_area += be.area
            elif type(be) in (BuildingElementAdjacentZTC,):
                continue
            else:
                sys.exit(f"Building element {type(be)} is not recognised. \
                    Please update zone.py function total_heat_loss_area.")
        return total_heat_loss_area

    def total_heat_capacity(self) -> cython.double:
        """ Return the total heat capacity from all building elements in a zone
        excluding ground and transparent elements, in kJ / K """
        # TODO Exclude solid door (opaque building element), or define convention
        #      that heat capacity of solid doors can be entered as zero
        total_heat_capacity: cython.double = 0
        for be in self.__building_elements:
            total_heat_capacity += be.heat_capacity()
        return total_heat_capacity

    def total_thermal_bridges(self) -> cython.double:
        """ Return the total heat transfer coefficient for all
        thermal bridges in a zone, in W / K """
        return self.__tb_heat_trans_coeff

    def total_vent_heat_loss(self) -> cython.double:
        """ Return the ventilation heat loss from all ventilation elements, in W / K """
        total_vent_heat_loss: cython.double = 0
        # TODO This doesn't work yet
        # if self.__vent_obj is not None:
        #     total_vent_heat_loss = self.__vent_obj.h_ve_average(self.__volume)
        return total_vent_heat_loss