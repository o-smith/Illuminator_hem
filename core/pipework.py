#!/usr/bin/env python3

"""
This module provides object(s) to represent pipework
"""

# Standard library imports
from math import pi, log
import sys
from enum import Enum, auto

# Local imports
import core.units as units
import core.material_properties as material_properties


# Set default values for the heat transfer coefficients inside the pipe, in W / m^2 K
INTERNAL_HTC_AIR = 15.5 # CIBSE Guide C, Table 3.25, air flow rate approx 3 m/s
INTERNAL_HTC_WATER = 1500.0 # CIBSE Guide C, Table 3.32 #Note, consider changing to 1478.4

# Set default values for the heat transfer coefficient at the outer surface, in W / m^2 K
EXTERNAL_REFLECTIVE_HTC = 5.7 # low emissivity reflective surface, CIBSE Guide C, Table 3.25
EXTERNAL_NONREFLECTIVE_HTC = 10.0 # high emissivity non-reflective surface, CIBSE Guide C, Table 3.25

class Location(Enum):
    EXTERNAL = auto()
    INTERNAL = auto()

    #there's a newer way of doing this
    @classmethod
    def from_string(cls, strval):
        if strval == 'external':
            return cls.EXTERNAL
        elif strval == 'internal':
            return cls.INTERNAL
        else:
            sys.exit('Location (' + str(strval) + ') not valid.')

class Pipework_Simple:
    """ An object to represent heat loss from pipework after flow has stopped """

    def __init__(self, location: Enum, internal_diameter: float, length: float, contents: str):
        """Construct a Pipework_Simple object

        Arguments:
        location              -- location of the pipework
        internal_diameter     -- internal diameter of the pipe, in m
        length                -- length of pipe, in m
        contents              -- whether the pipe is carrying air or water
        """
        self.__location: Enum = location
        self.__length: float = length
        self.__internal_diameter: float = internal_diameter
        self.__volume_litres: float \
            = pi * (self.__internal_diameter/2) * (self.__internal_diameter/2) \
            * self.__length * units.litres_per_cubic_metre

        if contents == "water":
            self.__contents_properties: material_properties.MaterialProperties \
                = material_properties.WATER
        else:
            sys.exit("No properties available for specified pipe content (" + str(contents) + ")")

    def get_location(self):
        return self.__location

    def volume_litres(self):
        return self.__volume_litres

    def cool_down_loss(self, inside_temp, outside_temp):
        """Calculates the total heat loss from a full pipe from demand temp to ambient
        temp in kWh

        Arguments:
        inside_temp   -- temperature of water (or air) inside the pipe, in degrees C
        outside_temp  -- temperature outside the pipe, in degrees C
        """
        cool_down_loss: float \
            = self.__contents_properties.volumetric_energy_content_kWh_per_litre(
                inside_temp,
                outside_temp,
                ) \
            * self.__volume_litres

        return(cool_down_loss) # returns kWh


class Pipework(Pipework_Simple):
    """ An object to represent steady state heat transfer in a hollow cyclinder (pipe)
    with radial heat flow. Method taken from 2021 ASHRAE Handbook, Section 4.4.2 """

    def __init__(self, location, internal_diameter, external_diameter, length, k_insulation, thickness_insulation, reflective, contents):
        """Construct a Pipework object

        Arguments:
        location              -- location of the pipework
        internal_diameter     -- internal diameter of the pipe, in m
        external_diameter     -- external diameter of the pipe, in m
        length                -- length of pipe, in m
        k_insulation          -- thermal conductivity of the insulation, in W / m K
        thickness_insulation  -- thickness of the pipe insulation, in m
        reflective            -- whether the surface is reflective or not (boolean input)
        contents              -- whether the pipe is carrying air or water
        """
        if external_diameter <= internal_diameter:
            sys.exit("Pipework: external diameter must be greater than internal diameter")

        super().__init__(location, internal_diameter, length, contents)

        """ Set the heat transfer coefficient inside the pipe, in W / m^2 K """
        if contents == 'air':
            internal_htc = INTERNAL_HTC_AIR
        elif contents == 'water':
            internal_htc = INTERNAL_HTC_WATER
        else:
            sys.exit('Contents of pipe not valid.')
                # TODO Exit just the current case instead of whole program entirely?
                # TODO Add code to log error

        """ Set the heat transfer coefficient at the outer surface, in W / m^2 K """
        if reflective:
            external_htc = EXTERNAL_REFLECTIVE_HTC
        else:
            external_htc = EXTERNAL_NONREFLECTIVE_HTC

        """ Calculate the diameter of the pipe including the insulation (D_insulation), in m"""
        self.__D_insulation = external_diameter + (2.0 * thickness_insulation)

        """ Calculate the interior surface resistance, in K m / W  """
        self.__interior_surface_resistance = 1.0 / (internal_htc * pi * internal_diameter)

        """ Calculate the insulation resistance, in K m / W  """
        self.__insulation_resistance = log(self.__D_insulation / internal_diameter) / (2.0 * pi * k_insulation)

        """ Calculate the external surface resistance, in K m / W  """
        self.__external_surface_resistance = 1.0 / (external_htc * pi * self.__D_insulation)

    def heat_loss(self, inside_temp, outside_temp):
        """" Return the heat loss from the pipe for the current timestep

        Arguments:
        inside_temp    -- temperature of water (or air) inside the pipe, in degrees C
        outside_temp   -- temperature outside the pipe, in degrees C
        """
        # Calculate total thermal resistance
        total_resistance = self.__interior_surface_resistance + self.__insulation_resistance + self.__external_surface_resistance

        # Calculate the heat loss for the current timestep, in W
        heat_loss \
            = (inside_temp - outside_temp) / (total_resistance) * self._Pipework_Simple__length

        return heat_loss

    def temperature_drop(self, inside_temp, outside_temp):
        """ Calculates by how much the temperature of water in a full pipe will fall
        over the timestep.

        Arguments:
        inside_temp   -- temperature of water (or air) inside the pipe, in degrees C
        outside_temp  -- temperature outside the pipe, in degrees C
        """
        heat_loss_kWh = (units.seconds_per_hour * self.heat_loss(inside_temp, outside_temp)) / units.W_per_kW # heat loss for the one hour timestep in kWh

        temp_drop: float \
            = min(
                ( ( heat_loss_kWh * units.J_per_kWh ) \
                / ( self._Pipework_Simple__contents_properties.volumetric_heat_capacity() 
                  * self.volume_litres()
                  )
                ),
                inside_temp - outside_temp,
                )  # Q = C m ∆t
        # temperature cannot drop below outside temperature

        return(temp_drop) # returns DegC

