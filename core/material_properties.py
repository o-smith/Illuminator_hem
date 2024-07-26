#!/usr/bin/env python3

"""
This module contains data on the properties of materials, and classes to
organise this data.
"""

# Local imports
import core.units as units


class MaterialProperties:
    """ An object to hold properties of materials, both specified and derived """
    # TODO Consider splitting this into different kinds of materials (e.g. heat
    #      transfer fluid, heat conductor etc.) so that not everything has to
    #      be defined (or made optional, adding complexity to the code)

    def __init__(self, density=None, specific_heat_capacity=None):
        """ Construct a MaterialProperties object

        Arguments:
        density                  -- in kg / litre
        specific_heat_capacity   -- in J / (kg.K)

        Derived values:
        volumetric_heat_capacity -- in J / (litre.K)
        """
        self.__density                  = density
        self.__specific_heat_capacity   = specific_heat_capacity
        self.__volumetric_heat_capacity = specific_heat_capacity * density

    def density(self):
        """ Return density of the material, in kg / litre """
        return self.__density

    def density_kg_per_m3(self):
        return self.__density * units.litres_per_cubic_metre

    def specific_heat_capacity(self):
        """ Return specific heat capacity, in J / (kg.K) """
        return self.__specific_heat_capacity

    def specific_heat_capacity_kWh(self):
        """ Return specific heat capacity, in kWh / (kg.K) """
        return self.__specific_heat_capacity / units.J_per_kWh

    def volumetric_heat_capacity(self):
        """ Return volumetric heat capacity, in J / (litre.K) """
        return self.__volumetric_heat_capacity

    def volumetric_energy_content_J_per_litre(self, temp_high, temp_base):
        """ Return energy content of material, in J / litre

        Arguments:
        temp_high -- temperature for which energy content should be calculated, in deg C or K
        temp_base -- temperature which defines "zero energy", in same units as temp_high
        """
        return (temp_high - temp_base) * self.__volumetric_heat_capacity

    def volumetric_energy_content_kWh_per_litre(self, temp_high, temp_base):
        """ Return energy content of material, in kWh / litre

        Arguments:
        temp_high -- temperature for which energy content should be calculated, in deg C or K
        temp_base -- temperature which defines "zero energy", in same units as temp_high
        """
        return self.volumetric_energy_content_J_per_litre(temp_high, temp_base) / units.J_per_kWh


# Define materials
WATER = MaterialProperties(density=1.0, specific_heat_capacity=4184)
# TODO Check the figures above. Use more precise ones?
AIR = MaterialProperties(density=0.001204, specific_heat_capacity=1006.0)
