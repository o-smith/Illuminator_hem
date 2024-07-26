#!/usr/bin/env python3

"""
This module provides objects to represent thermal bridges.

All ThermalBridge objects must provide a function heat_trans_coeff() which
returns the heat transfer coefficient of the thermal bridge
"""

class ThermalBridgeLinear:
    """ A class to represent linear thermal bridges """

    def __init__(self, linear_therm_trans, length):
        """ Construct a ThermalBridgeLinear object.

        Arguments:
        linear_therm_trans -- linear thermal transmittance of the thermal bridge, in W / (m.K)
        length             -- length of the thermal bridge over which the
                              linear thermal transmittance applies, in m
        """
        self.__linear_therm_trans = linear_therm_trans
        self.__length = length

    def heat_trans_coeff(self):
        """ Return the heat transfer coefficient for this thermal bridge. """
        return self.__linear_therm_trans * self.__length


class ThermalBridgePoint:
    """ A class to represent point thermal bridges """

    def __init__(self, heat_transfer_coeff):
        """ Construct a ThermalBridgePoint object.

        Arguments:
        heat_transfer_coeff -- heat transfer coefficient of the thermal bridge, in W / K
        """
        self.__heat_trans_coeff = heat_transfer_coeff

    def heat_trans_coeff(self):
        """ Return the heat transfer coefficient for this thermal bridge. """
        return self.__heat_trans_coeff