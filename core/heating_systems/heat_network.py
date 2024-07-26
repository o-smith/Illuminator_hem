#!/usr/bin/env python3

"""
This module provides objects to model heat networks
"""

# Standard library imports
import sys

# Local imports
from core.material_properties import WATER
from core.units import hours_per_day, W_per_kW


class HeatNetworkService:
    """ A base class for objects representing services (e.g. water heating) provided by a heat network.

    This object encapsulates the name of the service, meaning that the system
    consuming the energy does not have to specify this on every call, and
    helping to enforce that each service has a unique name.

    Derived objects provide a place to handle parts of the calculation (e.g.
    distribution flow temperature) that may differ for different services.

    Separate subclasses need to be implemented for different types of service
    (e.g. HW and space heating). These should implement the following functions:
    - demand_energy(self, energy_demand)
    """

    def __init__(self, heat_network, service_name, control=None):
        """ Construct a HeatNetworkService object

        Arguments:
        heat_network -- reference to the HeatNetwork object providing the service
        service_name -- name of the service demanding energy from the boiler
        control -- reference to a control object which must implement is_on() func
        """
        self._heat_network = heat_network
        self._service_name = service_name
        self.__control = control

    def is_on(self):
        if self.__control is not None:
            service_on = self.__control.is_on()
        else:
            service_on = True
        return service_on


class HeatNetworkServiceWaterDirect(HeatNetworkService):
    """ An object to represent a water heating service provided by a heat network.

    This object contains the parts of the heat network calculation that are
    specific to providing hot water directly to the dwelling.
    """

    def __init__(self,
                 heat_network,
                 service_name,
                 temp_hot_water,
                 cold_feed,
                 simulation_time
                ):
        """ Construct a HeatNetworkWater object

        Arguments:
        heat_network       -- reference to the HeatNetwork object providing the service
        service_name       -- name of the service demanding energy from the heat network
        temp_hot_water     -- temperature of the hot water to be provided, in deg C
        cold_feed          -- reference to ColdWaterSource object
        simulation_time    -- reference to SimulationTime object
        """
        super().__init__(heat_network, service_name)

        self.__temp_hot_water = temp_hot_water
        self.__cold_feed = cold_feed
        self.__service_name = service_name
        self.__simulation_time = simulation_time

    def get_cold_water_source(self):
        return self.__cold_feed

    def get_temp_hot_water(self):
        return self.__temp_hot_water

    def demand_hot_water(self, volume_demanded_target):
        """ Demand energy for hot water (in kWh) from the heat network """
        # Calculate energy needed to meet hot water demand
        if 'temp_hot_water' in volume_demanded_target:
            volume_demanded = volume_demanded_target['temp_hot_water']['warm_vol']
        else:
            volume_demanded = 0.0

        
        energy_content_kWh_per_litre = WATER.volumetric_energy_content_kWh_per_litre(
            self.__temp_hot_water,
            self.__cold_feed.temperature()
            )
        energy_demand = volume_demanded * energy_content_kWh_per_litre 

        return self._heat_network._HeatNetwork__demand_energy(
            self.__service_name,
            energy_demand,
            )

class HeatNetworkServiceWaterStorage(HeatNetworkService):
    """ An object to represent a water heating service provided by a heat network.

    This object contains the parts of the heat network calculation that are
    specific to providing hot water to the dwelling via a hot water cylinder.
    """

    def __init__(
            self,
            heat_network,
            service_name,
            temp_hot_water,
            control = None,
            ):
        """ Construct a HeatNetworkWaterStorage object

        Arguments:
        heat_network -- reference to the HeatNetwork object providing the service
        service_name -- name of the service demanding energy from the heat network
        temp_hot_water -- temperature of the hot water to be provided, in deg C
        control -- reference to a control object which must implement is_on() func
        """
        super().__init__(heat_network, service_name, control)
        self.__temp_hot_water = temp_hot_water

        self.__service_name = service_name

    def demand_energy(self, energy_demand, temp_return):
        """ Demand energy (in kWh) from the heat network """
        if not self.is_on():
            return 0.0

        # Calculate energy needed to cover losses
        return self._heat_network._HeatNetwork__demand_energy(
            self.__service_name,
            energy_demand,
            )

    def energy_output_max(self, temp_return):
        """ Calculate the maximum energy output of the heat network"""
        if not self.is_on():
            return 0.0

        return self._heat_network._HeatNetwork__energy_output_max(self.__temp_hot_water)


class HeatNetworkServiceSpace(HeatNetworkService):
    """ An object to represent a space heating service provided by a heat network.

    This object contains the parts of the heat network calculation that are
    specific to providing space heating-.
    """
    def __init__(self, heat_network, service_name, control):
        """ Construct a HeatNetworkSpace object

        Arguments:
        heat_network -- reference to the HeatNetwork object providing the service
        service_name -- name of the service demanding energy from the heat network
        control -- reference to a control object which must implement is_on() and setpnt() funcs
        """
        super().__init__(heat_network, service_name, control)

        self.__service_name = service_name
        self.__control = control

    def demand_energy(self, energy_demand, temp_flow, temp_return):
        """ Demand energy (in kWh) from the heat network """
        if not self.is_on():
            return 0.0

        return self._heat_network._HeatNetwork__demand_energy(
            self.__service_name,
            energy_demand,
            )

    def energy_output_max(self, temp_output, temp_return_feed):
        """ Calculate the maximum energy output of the heat network"""
        if not self.is_on():
            return 0.0

        return self._heat_network._HeatNetwork__energy_output_max(temp_output)

    def temp_setpnt(self):
        return self.__control.setpnt()

    def in_required_period(self):
        return self.__control.in_required_period()


class HeatNetwork:
    """ An object to represent a heat network """

    def __init__(
            self, 
            power_max,
            daily_loss,
            building_level_distribution_losses,
            energy_supply,
            energy_supply_conn_name_auxiliary,
            energy_supply_conn_name_building_level_distribution_losses,
            simulation_time,
            ):
        """ Construct a HeatNetwork object

        Arguments:
        power_max -- maximum power output of HIU, in kW
        daily_loss -- daily loss from the HIU, in kWh/day
        building_level_distribution_losses -- building level distribution losses in Watts
        energy_supply       -- reference to EnergySupply object
        energy_supply_conn_name_auxiliary -- name to use for reporting auxiliary energy use
        energy_supply_conn_name_building_level_distribution_losses 
            -- name to use for reporting building level distribution losses energy use
        simulation_time     -- reference to SimulationTime object

        Other variables:
        energy_supply_connections -- dictionary with service name strings as keys and corresponding
                                     EnergySupplyConnection objects as values
        temp_hot_water            -- temperature of the hot water to be provided, in deg C
        cold_feed                 -- reference to ColdWaterSource object
        """
        self.__power_max = power_max
        self.__daily_loss = daily_loss
        self.__building_level_distribution_losses = building_level_distribution_losses
        self.__energy_supply = energy_supply
        self.__simulation_time = simulation_time
        self.__energy_supply_connections = {}
        self.__energy_supply_connection_aux \
            = self.__energy_supply.connection(energy_supply_conn_name_auxiliary)
        self.__energy_supply_connection_building_level_distribution_losses \
            = self.__energy_supply.connection(
                energy_supply_conn_name_building_level_distribution_losses
                )
        self.__total_time_running_current_timestep = 0.0

    def __create_service_connection(self, service_name):
        """ Create an EnergySupplyConnection for the service name given """
        # Check that service_name is not already registered
        if service_name in self.__energy_supply_connections.keys():
            sys.exit("Error: Service name already used: "+service_name)
            # TODO Exit just the current case instead of whole program entirely?

        # Set up EnergySupplyConnection for this service
        self.__energy_supply_connections[service_name] = \
            self.__energy_supply.connection(service_name)

    def create_service_hot_water_direct(
            self,
            service_name,
            temp_hot_water,
            cold_feed,
            ):
        """ Return a HeatNetworkSeriviceWaterDirect object and create an EnergySupplyConnection for it
        
        Arguments:
        service_name      -- name of the service demanding energy from the heat network
        temp_hot_water    -- temperature of the hot water to be provided, in deg C
        cold_feed         -- reference to ColdWaterSource object
        """
        self.__create_service_connection(service_name)

        return HeatNetworkServiceWaterDirect(
            self,
            service_name,
            temp_hot_water,
            cold_feed,
            self.__simulation_time
            )

    def create_service_hot_water_storage(self, service_name, temp_hot_water, control=None):
        """ Return a HeatNetworkSeriviceWaterStorage object and create an EnergySupplyConnection for it

        Arguments:
        service_name -- name of the service demanding energy from the heat network
        temp_hot_water -- temperature of the hot water to be provided, in deg C
        control -- reference to a control object which must implement is_on() func
        """
        self.__create_service_connection(service_name)

        return HeatNetworkServiceWaterStorage(self, service_name, temp_hot_water, control)

    def create_service_space_heating(self, service_name, control):
        """ Return a HeatNetworkServiceSpace object and create an EnergySupplyConnection for it

        Arguments:
        service_name -- name of the service demanding energy from the heat network
        control -- reference to a control object which must implement is_on() and setpnt() funcs
        """
        self.__create_service_connection(service_name)

        return HeatNetworkServiceSpace(self, service_name, control)

    def __energy_output_max(self, temp_output):
        """ Calculate the maximum energy output of the heat network, accounting
            for time spent on higher-priority services.

        Note: Call via a HeatNetworkService object, not directly.
        """
        timestep = self.__simulation_time.timestep()
        time_available = timestep - self.__total_time_running_current_timestep
        return self.__power_max * time_available

    def __demand_energy(
            self,
            service_name,
            energy_output_required,
            ):
        """ Calculate energy required by heat network to satisfy demand for the service indicated."""
        energy_output_max = self.__energy_output_max(None)
        if energy_output_max == 0.0:
            return 0.0
        energy_output_provided = max(0.0, min(energy_output_required, energy_output_max))
        self.__energy_supply_connections[service_name].demand_energy(energy_output_provided)

        time_available \
            = self.__simulation_time.timestep() - self.__total_time_running_current_timestep
        self.__total_time_running_current_timestep \
            += (energy_output_provided / energy_output_max) * time_available

        return energy_output_provided

    def timestep_end(self):
        """ Calculations to be done at the end of each timestep """
        # Energy required to overcome losses
        self.__energy_supply_connection_aux.demand_energy(self.HIU_loss())
        self.__energy_supply_connection_building_level_distribution_losses.demand_energy(
            self.building_level_loss()
            )

        #Variables below need to be reset at the end of each timestep
        self.__total_time_running_current_timestep = 0.0

    def HIU_loss(self):
        """ Standing heat loss from the HIU (heat interface unit) in kWh """
        # daily_loss to be sourced from the PCDB, in kWh/day
        return self.__daily_loss / hours_per_day * self.__simulation_time.timestep()

    def building_level_loss(self):
        """ Converts building level distribution loss from watts to kWh """
        return self.__building_level_distribution_losses \
             / W_per_kW * self.__simulation_time.timestep()

