#!/usr/bin/env python3

"""
This module provides object(s) to represent ductwork
"""

# Standard library imports
import sys
from math import pi, log
from enum import Enum, auto

# Set default value for the heat transfer coefficient inside the duct, in W / m^2 K 
INTERNAL_HTC = 15.5 # CIBSE Guide C, Table 3.25, air flow rate approx 3 m/s

# Set default values for the heat transfer coefficient at the outer surface, in W / m^2 K
EXTERNAL_REFLECTIVE_HTC = 5.7 # low emissivity reflective surface, CIBSE Guide C, Table 3.25
EXTERNAL_NONREFLECTIVE_HTC = 10.0 # high emissivity non-reflective surface, CIBSE Guide C, Table 3.25

class DuctType(Enum):
    INTAKE = auto()
    SUPPLY = auto()
    EXTRACT = auto()
    EXHAUST = auto()

    @classmethod
    def from_string(cls, strval):
        if strval == "intake":
            return cls.INTAKE
        elif strval == "supply":
            return cls.SUPPLY
        elif strval == "extract":
            return cls.EXTRACT
        elif strval == "exhaust":
            return cls.EXHAUST
        else:
            sys.exit("DuctType (" + str(strval) + ") not valid.")


class Ductwork:
    """ An object to represent ductwork for mechanical ventilation with heat recovery 
    (MVHR), assuming steady state heat transfer in, 1. A hollow cyclinder (duct)
    with radial heat flow and 2. A rectangular cross-section. ISO 12241:2022 """

    def __init__(self,
                 cross_section_shape,
                 duct_perimeter,
                 internal_diameter,
                 external_diameter,
                 length,
                 k_insulation,
                 thickness_insulation,
                 reflective,
                 duct_type,
                 MVHR_location,
                 mvhr_eff,
                 ):
        """Construct a ductwork object
        Arguments:
        cross_section_shape  -- whether cross-section of duct is circular or rectangular(sqaure)
        duct_perimeter       -- if ductwork is rectangular(sqaure) enter perimeter, in m
        internal_diameter    -- internal diameter of the duct, in m
        external_diameter    -- external diameter of the duct, in m
        length               -- length of duct, in m
        k_insulation         -- thermal conductivity of the insulation, in W / m K
        thickness_insulation -- thickness of the duct insulation, in m
        reflective           -- whether the outer surface of the duct is reflective (true) or not (false) (boolean input)
        duct_type            -- intake, supply, extract or exhaust
        MVHR_location        -- location of the MVHR unit (inside or outside the thermal envelope) 
        mvhr_eff             -- heat recovery efficiency of MVHR unit (0 to 1)
        """
        self.__length = length
        self.__MVHR_location = MVHR_location 
        self.__mvhr_eff = mvhr_eff
        self.__duct_type = DuctType.from_string(duct_type)

        """ Select the correct heat transfer coefficient for the outer surface, in W / m^2 K """
        if reflective:
            external_htc = EXTERNAL_REFLECTIVE_HTC
        else:
            external_htc = EXTERNAL_NONREFLECTIVE_HTC

        if cross_section_shape == "circular":
            """ Calculate the diameter of the duct including the insulation (D_ins), in m"""
            self.__D_ins = external_diameter + (2.0 * thickness_insulation)
    
            """ Calculate the interior linear surface resistance, in K m / W  """
            self.__internal_surface_resistance = 1.0 / (INTERNAL_HTC * pi * internal_diameter)
    
            """ Calculate the insulation linear thermal resistance, in K m / W  """
            self.__insulation_resistance = log(self.__D_ins / internal_diameter) / (2.0 * pi * k_insulation)

            """ Calculate the exterior linear surface resistance, in K m / W  """
            self.__external_surface_resistance = 1.0 / (external_htc * pi * self.__D_ins)

        elif cross_section_shape == "rectangular":
            """ Calculate the perimeter of the duct including the insulation, in m"""
            # the value 8 is specified in the standard ISO 12241:2022 and not assigned a description
            duct_perimeter_external = duct_perimeter + (8 * thickness_insulation)
    
            """ Calculate the interior linear surface resistance, in K m / W  """
            self.__internal_surface_resistance = 1.0 / (INTERNAL_HTC * duct_perimeter)
    
            """ Calculate the insulation linear thermal resistance, in K m / W  """
            self.__insulation_resistance = (2.0 * thickness_insulation) \
                                           / (k_insulation * (duct_perimeter + duct_perimeter_external))
    
            """ Calculate the exterior linear surface resistance, in K m / W  """
            self.__external_surface_resistance = 1.0 / (external_htc * duct_perimeter_external)
        else:
            sys.exit("Duct shape not valid")

    def get_duct_type(self):
        return self.__duct_type

    def duct_heat_loss(self, inside_temp, outside_temp):
        """" Return the heat loss for air inside the duct for the current timestep
        Arguments:
        inside_temp    -- temperature of air inside the duct, in degrees C
        outside_temp   -- temperature outside the duct, in degrees C
        """
        # Calculate total thermal resistance
        total_resistance = self.__internal_surface_resistance + self.__insulation_resistance + self.__external_surface_resistance

        # Calculate the heat loss, in W
        duct_heat_loss = (inside_temp - outside_temp) / (total_resistance) * self.__length

        return duct_heat_loss

    def total_duct_heat_loss(self, temp_indoor_air, temp_outdoor_air):
        """" Return the heat loss for air inside the duct for the current timestep
        Arguments:
        outside_temp       -- temperature outside the duct, in degrees C
        supply_duct_temp   -- temperature of air inside the supply duct, in degrees C
        extract_duct_temp  -- temperature of air inside the extract duct, in degrees C
        intake_duct_temp   -- temperature of air inside the intake duct, in degrees C
        exhaust_duct_temp  -- temperature of air inside the exhaust duct, in degrees C
        efficiency         -- heat recovery efficiency of MVHR
        """
        temp_diff = temp_indoor_air - temp_outdoor_air
        # Outside location
        # Air inside the duct loses heat, external environment gains heat
        # Loses energy to outside in extract duct - losses must be X by the efficiency of heat recovery
        # Loses energy to outside in supply duct - lose all because after MVHR unit
        if self.__MVHR_location == 'outside':
            if self.__duct_type in (DuctType.INTAKE, DuctType.EXHAUST):
                return 0.0
            elif self.__duct_type == DuctType.SUPPLY:
                supply_duct_temp = temp_outdoor_air + (self.__mvhr_eff * temp_diff)
                return self.duct_heat_loss(supply_duct_temp, temp_outdoor_air)
            elif self.__duct_type == DuctType.EXTRACT:
                return self.duct_heat_loss(temp_indoor_air, temp_outdoor_air) * self.__mvhr_eff
            else:
                sys.exit("Heat loss for DuctType (" + str(self.__duct_type) + ") not defined.")
        # Inside location
        # This will be a negative heat loss i.e. air inside the duct gains heat, dwelling loses heat
        # Gains energy from zone in intake duct - benefit of gain must be X by the efficiency of heat recovery
        # Gains energy from zone in exhaust duct
        elif self.__MVHR_location == 'inside':
            if self.__duct_type in (DuctType.SUPPLY, DuctType.EXTRACT):
                return 0.0
            elif self.__duct_type == DuctType.INTAKE:
                return self.duct_heat_loss(temp_outdoor_air, temp_indoor_air) * self.__mvhr_eff
            elif self.__duct_type == DuctType.EXHAUST:
                exhaust_duct_temp = temp_outdoor_air + ((1.0 - self.__mvhr_eff) * temp_diff)
                return self.duct_heat_loss(exhaust_duct_temp, temp_indoor_air)
            else:
                sys.exit("Heat loss for DuctType (" + str(self.__duct_type) + ") not defined.")
        else:
            sys.exit('MVHR location not valid.')
            # TODO Exit just the current case instead of whole program entirely?
            # TODO Add code to log an error
