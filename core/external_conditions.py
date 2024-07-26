#!/usr/bin/env python3

"""
This module provides object(s) to store and look up data on external conditions
(e.g. external air temperature)

Calculation of solar radiation on a surface of a given orientation and tilt is
based on BS EN ISO 52010-1:2017.
"""

# Standard library imports
import sys
from math import cos, sin, tan, pi, asin, acos, radians, degrees, exp, sqrt, floor
from itertools import product
from copy import deepcopy

# Local imports
import core.units as units

class ExternalConditions:
    """ An object to store and look up data on external conditions """

    def __init__(self,
            simulation_time,
            air_temps,
            wind_speeds,
            wind_directions,
            diffuse_horizontal_radiation,
            direct_beam_radiation,
            solar_reflectivity_of_ground,
            latitude,
            longitude,
            timezone,
            start_day,
            end_day,
            time_series_step,
            january_first,
            daylight_savings,
            leap_day_included,
            direct_beam_conversion_needed,
            shading_segments,
            ):
        """ Construct an ExternalConditions object

        Arguments:
        simulation_time -- reference to SimulationTime object
        air_temps       -- list of external air temperatures, in deg C (one entry per hour)
        wind_speeds     -- list of wind speeds, in m/s (one entry per hour)
        wind_directions -- list of wind directions in degrees where North=0, East=90, 
                            South=180, West=270. Values range: 0 to 360.
                            Wind direction is reported by the direction from which it originates.
                            E.g, a southernly (180 degree) wind blows from the south to the north.
        diffuse_horizontal_radiation    -- list of diffuse horizontal radiation values, in W/m2 (one entry per hour)
        direct_beam_radiation           -- list of direct beam radiation values, in W/m2 (one entry per hour)
        solar_reflectivity_of_ground    -- list of ground reflectivity values, 0 to 1 (one entry per hour)
        latitude        -- latitude of weather station, angle from south, in degrees (single value)
        longitude       -- longitude of weather station, easterly +ve westerly -ve, in degrees (single value)
        timezone        -- timezone of weather station, -12 to 12 (single value)
        start_day       -- first day of the time series, day of the year, 0 to 365 (single value)
        end_day         -- last day of the time series, day of the year, 0 to 365 (single value)
        time_series_step -- timestep of the time series data, in hours
        january_first   -- day of the week for January 1st, monday to sunday, 1 to 7 (single value)
        daylight_savings    -- handling of daylight savings time, (single value)
                            e.g. applicable and taken into account, 
                            applicable but not taken into account, 
                            not applicable
        leap_day_included   -- whether climate data includes a leap day, true or false (single value)
        direct_beam_conversion_needed -- A flag to indicate whether direct beam radiation from climate data needs to be 
                                        converted from horizontal to normal incidence. If normal direct beam radiation 
                                        values are provided then no conversion is needed.
        shading_segments -- data splitting the ground plane into segments (8-36) and giving height
                            and distance to shading objects surrounding the building
        """     

        self.__simulation_time  = simulation_time
        self.__air_temps        = air_temps
        self.__wind_speeds      = wind_speeds
        self.__wind_directions  = wind_directions
        self.__solar_reflectivity_of_ground = solar_reflectivity_of_ground
        self.__latitude = latitude # practical  range -90 to +90
        self.__longitude = longitude # practical range -180 to +180
        self.__timezone = timezone
        self.__start_day = start_day
        self.__end_day = end_day
        self.__january_first = january_first
        self.__daylight_savings = daylight_savings
        self.__leap_day_included = leap_day_included
        self.__direct_beam_conversion_needed = direct_beam_conversion_needed
        self.__shading_segments = shading_segments
        self.__time_series_step = time_series_step

        # Initialise results cache (to improve performance)
        self.__cached_results = {}
        self.__cached_timestep = None

        days_in_year = 366 if leap_day_included else 365
        hours_in_year = days_in_year * 24
        time_shift = self.__init_time_shift()

        # Calculate earth orbit deviation for each day of year
        earth_orbit_deviation = [
            self.__init_earth_orbit_deviation(current_day)
            for current_day in range(0, days_in_year)
            ]
        # Calculate extra terrestrial radiation
        self.__extra_terrestrial_radiation = [
            self.__init_extra_terrestrial_radiation(earth_orbit_deviation[current_day])
            for current_day in range(0, days_in_year)
            ]
        # Calculate solar declination for each day of year
        self.__solar_declination = [
            self.__init_solar_declination(earth_orbit_deviation[current_day])
            for current_day in range(0, days_in_year)
            ]
        # Calculate equation of time for each day of year
        equation_of_time = [
            self.__init_equation_of_time(current_day)
            for current_day in range(0, days_in_year)
            ]
        # Calculate solar time for each hour of year
        self.__solar_time = [
            self.__init_solar_time(
                floor(current_hour % 24),
                equation_of_time[floor(current_hour / 24)],
                time_shift,
                )
            for current_hour in range(0, hours_in_year)
            ]
        # Calculate solar hour angle for each hour of year
        self.__solar_hour_angle = [
            self.__init_solar_hour_angle(self.__solar_time[current_hour])
            for current_hour in range(0, hours_in_year)
            ]
        # Calculate solar altitude for each hour of year
        self.__solar_altitude = [
            self.__init_solar_altitude(
                self.__solar_declination[floor(current_hour / 24)],
                self.__solar_hour_angle[current_hour],
                )
            for current_hour in range(0, hours_in_year)
            ]
        # Calculate solar zenith angle for each hour of year
        self.__solar_zenith_angle = [
            self.__init_solar_zenith_angle(self.__solar_altitude[current_hour])
            for current_hour in range(0, hours_in_year)
            ]
        # Calculate solar azimuth angle for each hour of year
        self.__solar_azimuth_angle = [
            self.__init_solar_azimuth_angle(
                self.__solar_declination[floor(current_hour / 24)],
                self.__solar_hour_angle[current_hour],
                self.__solar_altitude[current_hour],
                )
            for current_hour in range(0, hours_in_year)
            ]
        # Calculate air mass for each hour of year
        self.__air_mass = [
            self.__init_air_mass(self.__solar_altitude[current_hour])
            for current_hour in range(0, hours_in_year)
            ]

        # Calculate direct beam radiation for each timestep
        simtime = deepcopy(self.__simulation_time)
        self.__direct_beam_radiation = [
            self.__init_direct_beam_radiation(
                direct_beam_radiation[simtime.time_series_idx(self.__start_day, self.__time_series_step)],
                self.__solar_altitude[simtime.current_hour()]
                )
            for _, _, _ in simtime
            ]
        # Calculate diffuse horizontal radiation for each timestep
        simtime = deepcopy(self.__simulation_time)
        self.__diffuse_horizontal_radiation = [
            diffuse_horizontal_radiation[
                simtime.time_series_idx(self.__start_day, self.__time_series_step)
                ]
            for _, _, _ in simtime
            ]
        # Calculate dimensionless clearness parameter for each timestep
        simtime = deepcopy(self.__simulation_time)
        dimensionless_clearness_parameter = [
            self.__init_dimensionless_clearness_parameter(
                self.__diffuse_horizontal_radiation[t_idx],
                self.__direct_beam_radiation[t_idx],
                self.__solar_altitude[simtime.current_hour()],
                )
            for t_idx, _, _ in simtime
            ]
        # Calculate dimensionless sky brightness parameter for each timestep
        simtime = deepcopy(self.__simulation_time)
        dimensionless_sky_brightness_parameter = [
            self.__init_dimensionless_sky_brightness_parameter(
                self.__air_mass[simtime.current_hour()],
                self.__diffuse_horizontal_radiation[t_idx],
                self.__extra_terrestrial_radiation[simtime.current_day()],
                )
            for t_idx, _, _ in simtime
            ]
        # Calculate circumsolar brightness coefficient, F1 for each timestep
        simtime = deepcopy(self.__simulation_time)
        self.__F1 = [
            self.__init_F1(
                dimensionless_clearness_parameter[t_idx],
                dimensionless_sky_brightness_parameter[t_idx],
                self.__solar_zenith_angle[simtime.current_hour()],
                )
            for t_idx, _, _ in simtime
            ]
        # Calculate horizontal brightness coefficient, F2 for each timestep
        simtime = deepcopy(self.__simulation_time)
        self.__F2 = [
            self.__init_F2(
                dimensionless_clearness_parameter[t_idx],
                dimensionless_sky_brightness_parameter[t_idx],
                self.__solar_zenith_angle[simtime.current_hour()],
                )
            for t_idx, _, _ in simtime
            ]

    def testoutput_setup(self,tilt,orientation):
        """ print output to a file for analysis """

        #call this function once at the start of the calculation to test outputs

        import time
        readable = time.ctime()
        with open("test_sunpath.txt", "a") as o:
            o.write("\n")
            o.write("\n")
            o.write("*****************")
            o.write("\n")
            o.write(readable)
            o.write("\n")
            o.write("latitude " + str(self.latitude()))
            o.write("\n")
            o.write("longitude " + str(self.longitude()))
            o.write("\n")
            o.write("day of year " + str(self.start_day()))
            o.write("\n")
            o.write("surface tilt " + str(tilt))
            o.write("\n")
            o.write("surface orientation " + str(orientation))
            o.write("\n")
            o.write("sim hour,solar time,s declination,s hour angle,s altitude,s azimuth,air mass,sun surface azimuth,direct irad,solar angle of incidence,ET irad,F1,F2,E,delta,a over b,diffuse irad,ground reflect irad,circumsolar,final diffuse,final direct")


    def testoutput(self,tilt,orientation):
        """ print output to a file for analysis """

        #call this function once during every timestep to test outputs
        current_timestep = self.__simulation_time.index()
        current_hour = self.__simulation_time.current_hour()
        current_day = self.__simulation_time.current_day()

        #write headers
        with open("test_sunpath.txt", "a") as o:
            o.write("\n")
            o.write(str(self.__simulation_time.hour_of_day()))
            o.write(",")
            o.write(str(self.__solar_time[current_hour]))
            o.write(",")
            o.write(str(self.__solar_declination[current_day]))
            o.write(",")
            o.write(str(self.__solar_hour_angle[current_hour]))
            o.write(",")
            o.write(str(self.__solar_altitude[current_hour]))
            o.write(",")
            o.write(str(self.__solar_azimuth_angle[current_hour]))
            o.write(",")
            o.write(str(self.__air_mass[current_hour]))
            o.write(",")
            o.write(str(self.sun_surface_azimuth(orientation)))
            o.write(",")
            o.write(str(self.direct_irradiance(tilt, orientation)))
            o.write(",")
            o.write(str(self.solar_angle_of_incidence(tilt, orientation)))
            o.write(",")
            o.write(str(self.__extra_terrestrial_radiation[current_day]))
            o.write(",")
            o.write(str(self.__F1[current_timestep]))
            o.write(",")
            o.write(str(self.__F2[current_timestep]))
            o.write(",")
            o.write(str(self.__dimensionless_clearness_parameter[current_timestep]))
            o.write(",")
            o.write(str(self.__dimensionless_sky_brightness_parameter[current_timestep]))
            o.write(",")
            o.write(str(self.a_over_b(tilt, orientation)))
            o.write(",")
            o.write(str(self.diffuse_irradiance(tilt, orientation)[0]))
            o.write(",")
            o.write(str(self.ground_reflection_irradiance(tilt)))
            o.write(",")
            o.write(str(self.circumsolar_irradiance(tilt, orientation)))
            o.write(",")
            o.write(str(self.calculated_diffuse_irradiance(tilt, orientation)))
            o.write(",")
            o.write(str(self.calculated_direct_irradiance(tilt, orientation)))

    def air_temp(self):
        """ Return the external air temperature for the current timestep """
        return self.__air_temps[self.__simulation_time.time_series_idx(self.__start_day, self.__time_series_step)]

    def air_temp_annual(self):
        """ Return the average air temperature for the year """
        assert len(self.__air_temps) == 8760 # Only works if data for whole year has been provided
        return sum(self.__air_temps) / len(self.__air_temps)

    def air_temp_monthly(self):
        """ Return the average air temperature for the current month """
        # Get start and end hours for current month
        idx_start, idx_end = self.__simulation_time.current_month_start_end_hour()
        # Get air temperatures for the current month
        air_temps_month = self.__air_temps[idx_start:idx_end]

        return sum(air_temps_month) / len(air_temps_month)

    def ground_temp(self):
        """ Return the external ground temperature for the current timestep """
        return self.__ground_temps[self.__simulation_time.time_series_idx(self.__start_day, self.__time_series_step)]

    def wind_speed(self):
        """ Return the wind speed for the current timestep """
        return self.__wind_speeds[self.__simulation_time.time_series_idx(self.__start_day, self.__time_series_step)]

    def wind_speed_annual(self):
        """ Return the average wind speed for the year """
        # Only works if data for whole year has been provided, so assert this is true
        assert len(self.__wind_speeds) \
            == units.hours_per_day * units.days_per_year / self.__time_series_step
        return sum(self.__wind_speeds) / len(self.__wind_speeds)

    def wind_direction(self):
        """ Return the wind direction for the current timestep """
        return self.__wind_directions[self.__simulation_time.time_series_idx(self.__start_day, self.__time_series_step)]

    def diffuse_horizontal_radiation(self):
        """ Return the diffuse_horizontal_radiation for the current timestep """
        return self.__diffuse_horizontal_radiation[self.__simulation_time.index()]

    def direct_beam_radiation(self):
        """ Return the direct_beam_radiation for the current timestep """
        return self.__direct_beam_radiation[self.__simulation_time.index()]
    
    def __init_direct_beam_radiation(self, raw_value, solar_altitude):
        # if the climate data to only provide direct horizontal (rather than normal:
        # If only direct (beam) solar irradiance at horizontal plane is available in the climatic data set,
        # it shall be converted to normal incidence by dividing the value by the sine of the solar altitude.
        """ ISO 52010 section 6.4.2
        TODO investigate the impact of these notes further. Applicable for weather from CIBSE file. 
        NOTE 1 If the solar altitude angle is low, this conversion is very sensative for tiny
        errors in the calculation of the solar altitude. Such tiny errors are feasible given the
        sensitivity for the parameters needed to calculate the solar angle and given the atmospheric
        refraction of solar radiation near the ground. there fore the value at normal incidence is
        preferred.
        NOTE 2 method 1 proved to be most effective in mid-latitude climates 
        other models might be more suitable for tropical climates.
        NOTE 3 if the solar altitude angle is low, the conversion from direct horizontal to direct
        normal beam irradiance is very sensative for tiny errors in the calculation of the 
        solar altitude."""
        if self.__direct_beam_conversion_needed:
            sin_asol = sin(radians(solar_altitude))
            #prevent division by zero error. if sin_asol = 0 then the sun is lower than the
            #horizon and there will be no direct radiation to convert
            if sin_asol > 0:
                Gsol_b = raw_value / sin_asol
            else:
                Gsol_b = raw_value # TODO should this be zero?
        else:
            Gsol_b = raw_value

        return Gsol_b

    def solar_reflectivity_of_ground(self):
        """ Return the solar_reflectivity_of_ground for the current timestep """
        return self.__solar_reflectivity_of_ground[self.__simulation_time.time_series_idx(self.__start_day, self.__time_series_step)]

    def latitude(self):
        """ Return the latitude """
        return self.__latitude

    def longitude(self):
        """ Return the longitude """
        return self.__longitude

    def timezone(self):
        """ Return the timezone """
        return self.__timezone

    def start_day(self):
        """ Return the start_day """
        return self.__start_day
        # TODO possibly this input sits better within simulation_time
        # but included here until final form decided
        # currently used as the current day value

    def end_day(self):
        """ Return the end_day """
        return self.__end_day
        # TODO possibly this input sits better within simulation_time
        # but included here until final form decided
        # not current used

    def january_first(self):
        """ Return the january_first """
        return self.__january_first
        # TODO possibly this input sits better within simulation_time
        # but included here until final form decided
        # not current used

    def daylight_savings(self):
        """ Return the daylight_savings """
        return self.__daylight_savings
        # TODO possibly this input sits better within simulation_time
        # but included here until final form decided
        # currently unclear whether this is reffering to a choice by the user
        # or a statement of the contents of the weather data file
        # not current used

    def leap_day_included(self):
        """ Return the leap_day_included """
        return self.__leap_day_included
        # TODO possibly this input sits better within simulation_time
        # but included here until final form decided
        # currently unclear whether this is reffering to a choice by the user
        # or a statement of the contents of the weather data file
        # not current used

    def __init_earth_orbit_deviation(self, current_day):
        """ Calculate Rdc, the earth orbit deviation, as a function of the day, in degrees """

        nday = current_day + 1
        # nday is the day of the year, from 1 to 365 or 366 (leap year)
        # Note that current_day function returns days numbered 0 to 364 or 365,
        # so we need to add 1 above

        Rdc = (360 / 365) * nday

        return Rdc

    def __init_solar_declination(self, earth_orbit_deviation):
        """ Calculate solar declination in degrees """

        Rdc = radians(earth_orbit_deviation)
        # note we convert to radians for the python cos & sin inputs in formula below

        solar_declination = 0.33281 - 22.984 * cos(Rdc) - 0.3499 * cos(2 * Rdc) \
                          - 0.1398 * cos(3 * Rdc) + 3.7872 * sin(Rdc) + 0.03205 \
                          * sin(2 * Rdc) + 0.07187 * sin(3 * Rdc)

        return solar_declination

    def __init_equation_of_time(self, current_day):
        """ Calculate the equation of time """

        """ 
        teq is the equation of time, in minutes;
        nday is the day of the year, from 1 to 365 or 366 (leap year)
        """    

        nday = current_day + 1
        # nday is the day of the year, from 1 to 365 or 366 (leap year)
        # Note that current_day function returns days numbered 0 to 364 or 365,
        # so we need to add 1 here

        # note we convert the values inside the cos() to radians for the python function
        # even though the 180 / pi is converting from radians into degrees
        # this way the formula remains consistent with as written in the ISO document

        if   nday < 21:
            teq = 2.6 + 0.44 * nday
        elif nday < 136:
            teq = 5.2 + 9.0 * cos((nday - 43) * 0.0357)
        elif nday < 241:
            teq = 1.4 - 5.0 * cos((nday - 135) * 0.0449)
        elif nday < 336:
            teq = -6.3 - 10.0 * cos((nday - 306) * 0.036)
        elif nday <= 366:
            teq = 0.45 * (nday - 359)
        else:
            sys.exit("Day of the year ("+str(nday)+") not valid")

        return teq

    def __init_time_shift(self):
        """ Calculate the time shift, in hours, resulting from the fact that the 
        longitude and the path of the sun are not equal

        NOTE Daylight saving time is disregarded in tshift which is time independent
        """

        tshift = self.timezone() - self.longitude() / 15
        return tshift

    def __init_solar_time(self, hour_of_day, equation_of_time, time_shift):
        """ Calculate the solar time, tsol, as a function of the equation of time, 
        the time shift and the hour of the day """

        """ 
        tsol is the solar time, in h
        nhour is the actual (clock) time for the location, the hour of the day, in h
        """
        nhour = hour_of_day + 1
        #note we +1 here because the simulation hour of day starts at 0
        #while the sun path standard hour of day starts at 1 (hour 0 to 1)
        tsol = nhour - (equation_of_time / 60) - time_shift

        return tsol

    def __init_solar_hour_angle(self, solar_time):
        """ Calculate the solar hour angle, in the middle of the 
        current hour as a function of the solar time """

        # TODO How is this to be adjusted for timesteps that are not hourly?
        # would allowing solar_time to be a decimal be all that is needed?

        """ 
        w is the solar hour angle, in degrees

        Notes from ISO 52020 6.4.1.5
        NOTE 1 The limitation of angles ranging between -180 and +180 degrees is 
        needed to determine which shading objects are in the direction of the sun; 
        see also the calculation of the azimuth angle of the sun in 6.4.1.7.
        NOTE 2 Explanation of "12.5": The hour numbers are actually hour sections: 
        the first hour section of a day runs from 0h to 1h. So, the average position 
        of the sun for the solar radiation measured during (solar) hour section N is 
        at (solar) time = (N -0,5) h of the (solar) day.
        """

        w = (180 / 12) * (12.5 - solar_time)

        if w > 180:
            w = w - 360
        elif w < -180:
            w = w + 360

        return w

    def __init_solar_altitude(self, solar_declination, solar_hour_angle):
        """  the angle between the solar beam and the horizontal surface, determined 
             in the middle of the current hour as a function of the solar hour angle, 
             the solar declination and the latitude """

        # TODO How is this to be adjusted for timesteps that are not hourly?
        # would allowing solar_time to be a decimal be all that is needed?

        """ 
        asol is the solar altitude angle, the angle between the solar beam 
        and the horizontal surface, in degrees;
        """

        # note that we convert to radians for the sin & cos python functions and then
        # we need to convert the result back to degrees after the arcsin transformation

        asol = asin(
            sin(radians(solar_declination)) * sin(radians(self.latitude())) \
          + cos(radians(solar_declination)) * cos(radians(self.latitude())) * cos(radians(solar_hour_angle))
          )

        if degrees(asol) < 0.0001:
            return (0)

        return degrees(asol)

    def __init_solar_zenith_angle(self, solar_altitude):
        """  the complementary angle of the solar altitude """

        zenith = 90 - solar_altitude

        return zenith

    def __init_solar_azimuth_angle(self, solar_declination, solar_hour_angle, solar_altitude):
        """  calculates the solar azimuth angle,
        angle from South, eastwards positive, westwards negative, in degrees """

        """
        NOTE The azimuth angles range between −180 and +180 degrees; this is needed to determine which shading 
        objects are in the direction of the sun
        """

        sin_aux1_numerator \
            = cos(radians(solar_declination)) \
            * sin(radians(180 - solar_hour_angle))

        cos_aux1_numerator = cos(radians(self.latitude())) \
                    * sin(radians(solar_declination)) + sin(radians(self.latitude())) \
                    * cos(radians(solar_declination)) \
                    * cos(radians(180 - solar_hour_angle))

        denominator = cos(asin(sin(radians(solar_altitude))))

        sin_aux1 = sin_aux1_numerator / denominator            
        cos_aux1 = cos_aux1_numerator / denominator 
        aux2 = degrees(asin(sin_aux1_numerator) / denominator)

        # BS EN ISO 52010-1:2017. Formula 16
        if (sin_aux1 >= 0 and cos_aux1 > 0):
            solar_azimuth = 180 - aux2
            if solar_azimuth < 0:
                solar_azimuth = -solar_azimuth
        elif cos_aux1 < 0:
            solar_azimuth = aux2
        else:
            solar_azimuth = -(180 + aux2)

        return solar_azimuth

    def __init_air_mass(self, solar_altitude):
        """  calculates the air mass, m, the distance the solar beam travels through the earth atmosphere.
        The air mass is determined as a function of the sine of the solar altitude angle """

        sa = solar_altitude

        if sa >= 10:
            m = 1 / sin(radians(sa))
        else:
            m = 1 / (sin(radians(sa)) + 0.15 * (sa + 3.885)**-1.253)

        return m

    def solar_angle_of_incidence(self, tilt, orientation):
        """  calculates the solar angle of incidence, which is the angle of incidence of the 
        solar beam on an inclined surface and is determined as function of the solar hour angle 
        and solar declination 
        
        Arguments:
        tilt           -- is the tilt angle of the inclined surface from horizontal, measured 
                          upwards facing, 0 to 180, in degrees;
        orientation    -- is the orientation angle of the inclined surface, expressed as the 
                          geographical azimuth angle of the horizontal projection of the inclined 
                          surface normal, -180 to 180, in degrees;
        """

        # set up the parameters first just to make the very long equation slightly more readable
        current_day = self.__simulation_time.current_day()
        solar_declination = self.__solar_declination[current_day]
        sin_dec = sin(radians(solar_declination))
        cos_dec = cos(radians(solar_declination))
        sin_lat = sin(radians(self.latitude()))
        cos_lat = cos(radians(self.latitude()))
        sin_t = sin(radians(tilt))
        cos_t = cos(radians(tilt))
        sin_o = sin(radians(orientation))
        cos_o = cos(radians(orientation))
        current_hour = self.__simulation_time.current_hour()
        solar_hour_angle = self.__solar_hour_angle[current_hour]
        sin_sha = sin(radians(solar_hour_angle))
        cos_sha = cos(radians(solar_hour_angle))

        solar_angle_of_incidence = acos( \
                                   sin_dec * sin_lat * cos_t \
                                 - sin_dec * cos_lat * sin_t * cos_o \
                                 + cos_dec * cos_lat * cos_t * cos_sha \
                                 + cos_dec * sin_lat * sin_t * cos_o * cos_sha \
                                 + cos_dec * sin_t * sin_o * sin_sha \
                                 )

        return degrees(solar_angle_of_incidence)

    def sun_surface_azimuth(self, orientation):
        """  calculates the azimuth angle between sun and the inclined surface,
        needed as input for the calculation of the irradiance in case of solar shading by objects 

        Arguments:

        orientation    -- is the orientation angle of the inclined surface, expressed as the 
                          geographical azimuth angle of the horizontal projection of the inclined 
                          surface normal, -180 to 180, in degrees;

        """

        test_angle = self.solar_hour_angle() - orientation

        if test_angle > 180:
            azimuth = -360 + test_angle
        elif test_angle < -180:
            azimuth = 360 + test_angle
        else:
            azimuth = test_angle

        return azimuth

    def sun_surface_tilt(self, tilt):
        """  calculates the tilt angle between sun and the inclined surface,
        needed as input for the calculation of the irradiance in case of solar shading by objects 

        Arguments:

        tilt           -- is the tilt angle of the inclined surface from horizontal, measured 
                          upwards facing, 0 to 180, in degrees;

        """

        test_angle = tilt - self.solar_zenith_angle()

        if test_angle > 180:
            sun_surface_tilt = -360 + test_angle
        elif test_angle < -180:
            sun_surface_tilt = 360 + test_angle
        else:
            sun_surface_tilt = test_angle

        return sun_surface_tilt

    # TODO section 6.4.2 of ISO 52010 is not implemented here as it relates to methods of
    # obtaiing the needed irradiance values if they are not available from the climatic dataset.
    # If we decide this might be useful later then my suggestion would be to implement this as a
    # preprocessing step so that the core calculation always recieves the 'correct' climate data.

    # TODO solar reflectivity of the ground is expected to be initially fixed at 0.2, as per the 
    # default listed in ISO 52010 Annex B. However, the implementation here allows for one value 
    # per time step so alternative methods can be used in the future. options include taking values
    # from a climatic dataset that contains them, basing the values on ground surface material,
    # or ground cover such as snow. This can be implemented in a preprocess rather than the core.

    def direct_irradiance(self, tilt, orientation):
        """  calculates the direct irradiance on the inclined surface, determined as function 
        of cosine of the solar angle of incidence and the direct normal (beam) solar irradiance
        NOTE The solar beam irradiance is defined as falling on an surface normal to the solar beam. 
        This is not the same as direct horizontal radiation.

        Arguments:
        tilt           -- is the tilt angle of the inclined surface from horizontal, measured 
                          upwards facing, 0 to 180, in degrees;
        orientation    -- is the orientation angle of the inclined surface, expressed as the 
                          geographical azimuth angle of the horizontal projection of the inclined 
                          surface normal, -180 to 180, in degrees;
                          
        """

        direct_irradiance = max(0, self.direct_beam_radiation() * cos(radians(self.solar_angle_of_incidence(tilt, orientation))))

        return direct_irradiance

    def __init_extra_terrestrial_radiation(self, earth_orbit_deviation):
        """  calculates the extra terrestrial radiation, the normal irradiance out of the atmosphere 
        as a function of the day

        Arguments:
        tilt           -- is the tilt angle of the inclined surface from horizontal, measured 
                          upwards facing, 0 to 180, in degrees;
        orientation    -- is the orientation angle of the inclined surface, expressed as the 
                          geographical azimuth angle of the horizontal projection of the inclined 
                          surface normal, -180 to 180, in degrees;
                          
        """

        #NOTE the ISO 52010 has an error in this formula.
        #it lists Gsol,c as the solar angle of incidence on the inclined surface
        #when it should be the solar constant, given elsewhere as 1367
        #we use the correct version of the formula here

        extra_terrestrial_radiation = 1367 * (1 + 0.033 * cos(radians(earth_orbit_deviation)))

        return extra_terrestrial_radiation

    def brightness_coefficient(self, E, Fij):
        """ returns brightness coefficient as a look up from Table 8 in ISO 52010 

        Arguments:
        E    -- dimensionless clearness parameter
        Fij  -- the coefficient to be returned. e.g. f12 or f23
        """

        #TODO I've not had a need for the clearness index parameters contained in this table yet, 
        #if they are needed as input or output later then this function can be reworked

        if E < 1.065:
            # overcast
            index = 1
        elif E < 1.23:
            index = 2
        elif E < 1.5:
            index = 3
        elif E < 1.95:
            index = 4
        elif E < 2.8:
            index = 5
        elif E < 4.5:
            index = 6
        elif E < 6.2:
            index = 7
        else:
            #clear
            index = 8

        brightness_coeff_dict = {
            1: {'f11': -0.008, 'f12': 0.588, 'f13': -0.062, 'f21': -0.06, 'f22': 0.072, 'f23': -0.022},
            2: {'f11': 0.13, 'f12': 0.683, 'f13': -0.151, 'f21': -0.019, 'f22': 0.066, 'f23': -0.029},
            3: {'f11': 0.33, 'f12': 0.487, 'f13': -0.221, 'f21': 0.055, 'f22': -0.064, 'f23': -0.026},
            4: {'f11': 0.568, 'f12': 0.187, 'f13': -0.295, 'f21': 0.109, 'f22': -0.152, 'f23': -0.014},
            5: {'f11': 0.873, 'f12': -0.392, 'f13': -0.362, 'f21': 0.226, 'f22': -0.462, 'f23': 0.001},
            6: {'f11': 1.132, 'f12': -1.237, 'f13': -0.412, 'f21': 0.288, 'f22': -0.823, 'f23': 0.056},
            7: {'f11': 1.06, 'f12': -1.6, 'f13': -0.359, 'f21': 0.264, 'f22': -1.127, 'f23': 0.131},
            8: {'f11': 0.678, 'f12': -0.327, 'f13': -0.25, 'f21': 0.156, 'f22': -1.377, 'f23': 0.251}
            }

        return brightness_coeff_dict[index][Fij]

    def __init_F1(self, E, delta, solar_zenith_angle):
        """ returns the circumsolar brightness coefficient, F1

        Arguments:
        E -- dimensionless clearness parameter for the current timestep
        delta -- dimensionless sky brightness parameter for the current timestep
        solar_zenith_angle -- solar zenith angle for the current hour
        """

        #brightness coeffs
        f11 = self.brightness_coefficient(E, 'f11')
        f12 = self.brightness_coefficient(E, 'f12')
        f13 = self.brightness_coefficient(E, 'f13')
        #The formulation of F1 is made so as to avoid non-physical negative values 
        #that may occur and result in unacceptable distortions if the model is used 
        #for very low solar elevation angles
        F1 = max(0, f11 + f12 * delta + f13 * (pi * solar_zenith_angle / 180))

        return F1

    def __init_F2(self, E, delta, solar_zenith_angle):
        """ returns the horizontal brightness coefficient, F2

        Arguments:
        E -- dimensionless clearness parameter
        delta -- dimensionless sky brightness parameter
        solar_zenith_angle -- solar zenith angle for the current hour
        """

        #horizontal brightness coefficient, F2
        f21 = self.brightness_coefficient(E, 'f21')
        f22 = self.brightness_coefficient(E, 'f22')
        f23 = self.brightness_coefficient(E, 'f23')
        #F2 does not have the same restriction of max 0 as F1
        #from the EnergyPlus Engineering Reference:
        #The horizon brightening is assumed to be a linear source at the horizon 
        #and to be independent of azimuth. In actuality, for clear skies, the 
        #horizon brightening is highest at the horizon and decreases in intensity 
        #away from the horizon. For overcast skies the horizon brightening has a 
        #negative value since for such skies the sky radiance increases rather than 
        #decreases away from the horizon.
        F2 = f21 + f22 * delta + f23 * (pi * solar_zenith_angle / 180)

        return F2

    def __init_dimensionless_clearness_parameter(self, Gsol_d, Gsol_b, asol):
        """ returns the dimensionless clearness parameter, E, anisotropic sky conditions (Perez model)
        
        Arguments:
        Gsol_d -- diffuse horizontal radiation
        Gsol_b -- direct beam radiation
        asol -- solar altitude for the current hour
        """

        #constant parameter for the clearness formula, K, in rad^-3 from table 9 of ISO 52010
        K = 1.014

        if Gsol_d == 0:
            E = 999
        else:
            E = (((Gsol_d + Gsol_b) / Gsol_d) + K * (pi / 180 * asol)**3) \
              / (1 + K * (pi / 180 * asol)**3)

        return E

    def __init_dimensionless_sky_brightness_parameter(
            self,
            air_mass,
            diffuse_horizontal_radiation,
            extra_terrestrial_radiation,
            ):
        """  calculates the dimensionless sky brightness parameter, delta

        Arguments:
        air_mass -- air mass for the current hour
        diffuse_horizontal_radiation -- diffuse horizontal radiation for the current timestep
        extra_terrestrial_radiation -- extra-terrestrial radiation for the current day
        """

        delta = air_mass * diffuse_horizontal_radiation \
              / extra_terrestrial_radiation

        return delta

    def a_over_b(self, tilt, orientation):
        """  calculates the ratio of the parameters a and b

        Arguments:
        tilt           -- is the tilt angle of the inclined surface from horizontal, measured 
                          upwards facing, 0 to 180, in degrees;
        orientation    -- is the orientation angle of the inclined surface, expressed as the 
                          geographical azimuth angle of the horizontal projection of the inclined 
                          surface normal, -180 to 180, in degrees;
        """

        current_hour = self.__simulation_time.current_hour()
        #dimensionless parameters a & b
        #describing the incidence-weighted solid angle sustained by the circumsolar region as seen 
        #respectively by the tilted surface and the horizontal. 
        a = max(0, cos(radians(self.solar_angle_of_incidence(tilt, orientation))))
        b = max(cos(radians(85)), cos(radians(self.__solar_zenith_angle[current_hour])))

        return a / b


    def diffuse_irradiance(self, tilt, orientation):
        """  calculates the diffuse part of the irradiance on the surface (without ground reflection)

        Arguments:
        tilt           -- is the tilt angle of the inclined surface from horizontal, measured 
                          upwards facing, 0 to 180, in degrees;
        orientation    -- is the orientation angle of the inclined surface, expressed as the 
                          geographical azimuth angle of the horizontal projection of the inclined 
                          surface normal, -180 to 180, in degrees;
        """

        #first set up parameters needed for the calculation
        Gsol_d = self.diffuse_horizontal_radiation()
        F1 = self.__F1[self.__simulation_time.index()]
        F2 = self.__F2[self.__simulation_time.index()]

        # Calculate components of diffuse radiation
        diffuse_irr_sky = Gsol_d * (1 - F1) * ((1 + cos(radians(tilt))) / 2)
        diffuse_irr_circumsolar = self.circumsolar_irradiance(tilt, orientation)
        diffuse_irr_horiz = Gsol_d * F2 * sin(radians(tilt))

        diffuse_irr_total = diffuse_irr_sky + diffuse_irr_circumsolar + diffuse_irr_horiz

        return diffuse_irr_total, diffuse_irr_sky, diffuse_irr_circumsolar, diffuse_irr_horiz

    def ground_reflection_irradiance(self, tilt):
        """  calculates the contribution of the ground reflection to the irradiance on the inclined surface, 
        determined as function of global horizontal irradiance, which in this case is calculated from the solar 
        altitude, diffuse and beam solar irradiance and the solar reflectivity of the ground

        Arguments:
        tilt           -- is the tilt angle of the inclined surface from horizontal, measured 
                          upwards facing, 0 to 180, in degrees;
        """

        #first set up parameters needed for the calculation
        current_hour = self.__simulation_time.current_hour()
        Gsol_d = self.diffuse_horizontal_radiation()
        Gsol_b = self.direct_beam_radiation()
        asol = radians(self.__solar_altitude[current_hour])

        ground_reflection_irradiance = (Gsol_d + Gsol_b * sin(asol)) * self.solar_reflectivity_of_ground() \
                                     * ((1 - cos(radians(tilt))) / 2)

        return ground_reflection_irradiance

    def circumsolar_irradiance(self, tilt, orientation):
        """  calculates the circumsolar_irradiance

        Arguments:
        tilt           -- is the tilt angle of the inclined surface from horizontal, measured 
                          upwards facing, 0 to 180, in degrees;
        orientation    -- is the orientation angle of the inclined surface, expressed as the 
                          geographical azimuth angle of the horizontal projection of the inclined 
                          surface normal, -180 to 180, in degrees;
        """

        Gsol_d = self.diffuse_horizontal_radiation()
        F1 = self.__F1[self.__simulation_time.index()]
        a_over_b = self.a_over_b(tilt, orientation)

        circumsolar_irradiance = Gsol_d * F1 * a_over_b

        return circumsolar_irradiance

    def calculated_direct_irradiance(self, tilt, orientation):
        """  calculates the total direct irradiance on an inclined surface including circumsolar

        Arguments:
        tilt           -- is the tilt angle of the inclined surface from horizontal, measured 
                          upwards facing, 0 to 180, in degrees;
        orientation    -- is the orientation angle of the inclined surface, expressed as the 
                          geographical azimuth angle of the horizontal projection of the inclined 
                          surface normal, -180 to 180, in degrees;
        """

        calculated_direct = self.direct_irradiance(tilt, orientation) + self.circumsolar_irradiance(tilt, orientation)

        return calculated_direct

    def calculated_diffuse_irradiance(self, tilt, orientation):
        """  calculates the total diffuse irradiance on an inclined surface excluding circumsolar
        and including ground reflected irradiance

        Arguments:
        tilt           -- is the tilt angle of the inclined surface from horizontal, measured 
                          upwards facing, 0 to 180, in degrees;
        orientation    -- is the orientation angle of the inclined surface, expressed as the 
                          geographical azimuth angle of the horizontal projection of the inclined 
                          surface normal, -180 to 180, in degrees;
        """

        diffuse_irr_total, _, diffuse_irr_circumsolar, _ \
            = self.diffuse_irradiance(tilt, orientation)

        calculated_diffuse = diffuse_irr_total \
                           - diffuse_irr_circumsolar \
                           + self.ground_reflection_irradiance(tilt)

        return calculated_diffuse

    def calculated_total_solar_irradiance(self, tilt, orientation):
        """  calculates the hemispherical or total solar irradiance on the inclined surface 
        without the effect of shading

        Arguments:
        tilt           -- is the tilt angle of the inclined surface from horizontal, measured 
                          upwards facing, 0 to 180, in degrees;
        orientation    -- is the orientation angle of the inclined surface, expressed as the 
                          geographical azimuth angle of the horizontal projection of the inclined 
                          surface normal, -180 to 180, in degrees;
                          
        """

        total_irradiance = self.calculated_direct_irradiance(tilt, orientation) \
                         + self.calculated_diffuse_irradiance(tilt, orientation)

        return total_irradiance

    def calculated_direct_diffuse_total_irradiance(self, tilt, orientation, diffuse_breakdown=False):
        t_idx = self.__simulation_time.index()
        if t_idx != self.__cached_timestep:
            # If we have moved on to a new timestep, then clear the cached results
            self.__cached_results = {}
            self.__cached_timestep = t_idx

        if (tilt, orientation) in self.__cached_results.keys():
            # Look up cached values if this tilt and orientation has already been calculated
            calculated_direct = self.__cached_results[(tilt, orientation)]['calculated_direct']
            calculated_diffuse = self.__cached_results[(tilt, orientation)]['calculated_diffuse']
            total_irradiance = self.__cached_results[(tilt, orientation)]['total_irradiance']
            diffuse_res_breakdown = self.__cached_results[(tilt, orientation)]['diffuse_res_breakdown']
        else:
            # Calculate results if this tilt and orientation has not already been calculated
            diffuse_irr_total, diffuse_irr_sky, diffuse_irr_circumsolar, diffuse_irr_horiz \
                = self.diffuse_irradiance(tilt, orientation)
            ground_refl_irr = self.ground_reflection_irradiance(tilt)

            calculated_direct = self.direct_irradiance(tilt, orientation) + diffuse_irr_circumsolar
            calculated_diffuse = diffuse_irr_total \
                               - diffuse_irr_circumsolar \
                               + ground_refl_irr
            total_irradiance = calculated_direct \
                             + calculated_diffuse

            diffuse_res_breakdown = {
                'sky': diffuse_irr_sky,
                'circumsolar': diffuse_irr_circumsolar,
                'horiz': diffuse_irr_horiz,
                'ground_refl': ground_refl_irr,
                }

            # Cache calculated results
            self.__cached_results[(tilt, orientation)] = {}
            self.__cached_results[(tilt, orientation)]['calculated_direct'] = calculated_direct
            self.__cached_results[(tilt, orientation)]['calculated_diffuse'] = calculated_diffuse
            self.__cached_results[(tilt, orientation)]['total_irradiance'] = total_irradiance
            self.__cached_results[(tilt, orientation)]['diffuse_res_breakdown'] = diffuse_res_breakdown

        if diffuse_breakdown:
            return calculated_direct, calculated_diffuse, total_irradiance, diffuse_res_breakdown
        else:
            return calculated_direct, calculated_diffuse, total_irradiance

    # end of sun path calculations from ISO 52010
    # below are overshading calculations from ISO 52016

    def outside_solar_beam(self, tilt, orientation):
        """ checks if the shaded surface is in the view of the solar beam.
        if not, then shading is complete, total direct rad = 0 and no further
        shading calculation needed for this object for this time step. returns
        a flag for whether the surface is outside solar beam

        Arguments:
        tilt           -- is the tilt angle of the inclined surface from horizontal, measured 
                          upwards facing, 0 to 180, in degrees;
        orientation    -- is the orientation angle of the inclined surface, expressed as the 
                          geographical azimuth angle of the horizontal projection of the 
                          inclined surface normal, -180 to 180, in degrees;
                          
        """

        current_hour = self.__simulation_time.current_hour()
        test1 = orientation - self.__solar_azimuth_angle[current_hour]
        test1 = test1 - 360.0 if test1 > +180.0 \
           else test1 + 360.0 if test1 < -180.0 \
           else test1
        test2 = tilt - self.__solar_altitude[current_hour]

        if (-90 > test1 or test1 > 90):
            # surface outside solar beam
            return 1
        elif (-90 > test2 or test2 > 90):
            # surface outside solar beam
            return 1
        else:
            # surface inside solar beam
            return 0

    def get_segment(self):
        """ for complex (environment) shading objects, we need to know which
        segment the azimuth of the sun occupies at each timestep

        """

        current_hour = self.__simulation_time.current_hour()
        azimuth = self.__solar_azimuth_angle[current_hour]
        
        previous_segment_end = None
        for segment in self.__shading_segments:
            if previous_segment_end != None and previous_segment_end != segment["start"]: 
                sys.exit("Error: No gaps between segments allowed")
            previous_segment_end = segment["end"]
            if segment["end"] > segment["start"]:
                sys.exit("Error: End orientation is less than the start orientation. Check shading inputs")
            if (azimuth < segment["start"] and azimuth > segment["end"]):
                return segment
        #if not exited function yet then segment has not been found and there
        #is some sort of error
        sys.exit("solar segment not found. Check shading inputs")

    def obstacle_shading_height(self, Hkbase, Hobst, Lkobst ):
        """ calculates the height of the shading on the shaded surface (k),
        from the shading obstacle in segment i at time t. Note that "obstacle"
        has a specific meaning in ISO 52016 Annex F

        Arguments:
        Hkbase        -- is the base height of the shaded surface k, in m
        Hobst         -- is the height of the shading obstacle, p, in segment i, in m
        Lkobst        -- is the horizontal distance between the shaded surface k, in m 
                         and the shading obstacle p in segment i, in m
        """

        current_hour = self.__simulation_time.current_hour()
        Hshade = max(0, Hobst - Hkbase - Lkobst * tan(radians(self.__solar_altitude[current_hour])))
        return Hshade

    def overhang_shading_height(self, Hk, Hkbase, Hovh, Lkovh ):
        """ calculates the height of the shading on the shaded surface (k),
        from the shading overhang in segment i at time t. Note that "overhang"
        has a specific meaning in ISO 52016 Annex F

        Arguments:
        Hk            -- is the height of the shaded surface, k, in m
        Hkbase        -- is the base height of the shaded surface k, in m
        Hovh          -- is the lowest height of the overhang q, in segment i, in m
        Lkovh         -- is the horizontal distance between the shaded surface k 
                         and the shading overhang, q, in segment i, in m
        """

        Hshade = max(0, Hk + Hkbase - Hovh + Lkovh * tan(radians(self.solar_altitude())))
        return Hshade

    def direct_shading_reduction_factor(self, base_height, height, width, orientation, window_shading):
        """ calculates the shading factor of direct radiation due to external
        shading objects

        Arguments:
        height         -- is the height of the shaded surface (if surface is tilted then
                          this must be the vertical projection of the height), in m
        base_height    -- is the base height of the shaded surface k, in m
        width          -- is the width of the shaded surface, in m
        orientation    -- is the orientation angle of the inclined surface, expressed as the 
                          geographical azimuth angle of the horizontal projection of the 
                          inclined surface normal, -180 to 180, in degrees;
        window_shading -- data on overhangs and side fins associated to this building element
                          includes the shading object type, depth, anf distance from element
        """

        # start with default assumption of no shading
        Hshade_obst = 0
        Hshade_ovh = 0
        WfinR = 0
        WfinL = 0

        #first process the distant (environment) shading for this building element

        #get the shading segment we are currently in
        segment = self.get_segment()
        #check for any shading objects in this segment
        if "shading" in segment.keys():
            for shade_obj in segment["shading"]:
                if shade_obj["type"] == "obstacle":
                    new_shade_height = self.obstacle_shading_height \
                    (base_height, shade_obj["height"], shade_obj["distance"])

                    Hshade_obst = max(Hshade_obst, new_shade_height)
                elif shade_obj["type"] == "overhang":
                    new_shade_height = self.overhang_shading_height \
                    (height, base_height, shade_obj["height"], shade_obj["distance"])

                    Hshade_ovh = max(Hshade_ovh, new_shade_height)
                else:
                    sys.exit("shading object type" + shade_obj["type"] + "not recognised")

        # then check if there is any simple shading on this building element
        # (note only applicable to transparent building elements so window_shading
        # will always be False for other elements)
        if window_shading:
            current_hour = self.__simulation_time.current_hour()
            altitude = self.__solar_altitude[current_hour]
            azimuth = self.__solar_azimuth_angle[current_hour]
            # if there is then loop through all objects and calc shading heights/widths
            for shade_obj in window_shading:
                depth = shade_obj["depth"]
                distance = shade_obj["distance"]
                if shade_obj["type"] == "overhang":
                    new_shade_height = (depth * tan(radians(altitude)) \
                                    / cos(radians(azimuth - orientation))) \
                                    - distance

                    Hshade_ovh = max(Hshade_ovh, new_shade_height)
                elif shade_obj["type"] == "sidefinright":
                    #check if the sun is in the opposite direction
                    check = azimuth - orientation
                    if check > 0:
                        new_finRshade = 0
                    else:
                        new_finRshade = depth * tan(radians(azimuth - orientation)) \
                                        - distance
                    WfinR = max(WfinR, new_finRshade)
                elif shade_obj["type"] == "sidefinleft":
                    #check if the sun is in the opposite direction
                    check = azimuth - orientation
                    if check < 0:
                        new_finLshade = 0
                    else:
                        new_finLshade = depth * tan(radians(azimuth - orientation)) \
                                        - distance
                    WfinL = max(WfinL, new_finLshade)
                else:
                    sys.exit("shading object type" + shade_obj["type"] + "not recognised")

        # The height of the shade on the shaded surface from all obstacles is the 
        # largest of all, with as maximum value the height of the shaded object
        Hk_obst = min(height, Hshade_obst)

        # The height of the shade on the shaded surface from all overhangs is the 
        # largest of all, with as maximum value the height of the shaded object
        Hk_ovh = min(height, Hshade_ovh)

        # The height of the remaining sunlit area on the shaded surface from 
        # all obstacles and all overhangs
        Hk_sun = max(0, height - (Hk_obst + Hk_ovh))

        # The width of the shade on the shaded surface from all right side fins 
        # is the largest of all, with as maximum value the width of the shaded object
        Wk_finR = min(width, WfinR)

        # The width of the shade on the shaded surface from all left side fins 
        # is the largest of all, with as maximum value the width of the shaded object
        Wk_finL = min(width, WfinL)

        # The width of the remaining sunlit area on the shaded surface from all 
        # right hand side fins and all left hand side fins
        Wk_sun = max(0, width - (Wk_finR + Wk_finL))

        # And then the direct shading reduction factor of the shaded surface for 
        # obstacles, overhangs and side fins
        Fdir = (Hk_sun * Wk_sun) / (height * width)

        return Fdir

    def diffuse_shading_reduction_factor(
            self,
            diffuse_breakdown,
            tilt,
            height,
            width,
            window_shading,
            ):
        """ calculates the shading factor of diffuse radiation due to external
        shading objects

        Arguments:
        height         -- is the height of the shaded surface (if surface is tilted then
                          this must be the vertical projection of the height), in m
        base_height    -- is the base height of the shaded surface k, in m
        width          -- is the width of the shaded surface, in m
        orientation    -- is the orientation angle of the inclined surface, expressed as the 
                          geographical azimuth angle of the horizontal projection of the 
                          inclined surface normal, -180 to 180, in degrees;
        window_shading -- data on overhangs and side fins associated to this building element
                          includes the shading object type, depth, and distance from element
        """
        # Note: Shading factor for circumsolar radiation is same as for direct.
        #       As circumsolar radiation will be subtracted from diffuse and
        #       added to direct later on, we don't need to do anything for
        #       circumsolar radiation here and it is excluded.
        diffuse_irr_sky = diffuse_breakdown['sky']
        diffuse_irr_hor = diffuse_breakdown['horiz']
        diffuse_irr_ref = diffuse_breakdown['ground_refl']
        diffuse_irr_total = diffuse_irr_sky + diffuse_irr_hor + diffuse_irr_ref

        # TODO Calculate shading from distant objects (PD CEN ISO/TR 52016-2:2017 Section F.6.2)
        # TODO Loop over segments
        # TODO     Calculate sky-diffuse shading factor for the segment
        # TODO     Calculate horiz-diffuse shading factor for the segment
        # TODO Sum product of shading factor and irradiance for each element

        # Calculate shading from fins and overhangs (PD CEN ISO/TR 52016-2:2017 Section F.6.3)
        # TODO Alpha is not defined in this standard but is possibly the angular
        #      height of the horizon. This would make some sense as raising the
        #      horizon angle would have essentially the same effect as increasing
        #      the tilt of the building element so it would make sense to add it
        #      to beta when calculating the sky view factor. The overall effect
        #      of changing alpha when there are no fins or overhangs seems to be
        #      to decrease the shading factor (meaning more shading) for diffuse
        #      radiation from sky and horizon, and to increase the shading factor
        #      for radiation reflected from the ground (sometimes to above 1),
        #      which would also seem to make sense as in this case there is less
        #      sky and more ground in view than in the basic assumption of
        #      perfectly flat surroundings.
        # TODO Should angular height of horizon be user input or derived from
        #      calculation of distant shading objects above? Set to zero for now
        angular_height_of_horizon = 0.0
        alpha = radians(angular_height_of_horizon)
        # TODO Beta is not defined in this standard but is used for tilt of the
        #      building element in BS EN ISO 52016-1:2017, so assuming the same.
        #      This seems to give sensible numbers for the sky view factor
        #      F_w_sky when alpha = 0
        beta = radians(tilt)

        #create lists of diffuse shading factors to keep the largest one
        #in case there are multiple shading objects
        Fdiff_list = []

        # Unpack window shading details
        ovh_D_L_ls = [[0.0,1.0]] # [D,L] - L cannot be zero as this leads to divide-by-zero later on
        finR_D_L_ls = [[0.0,1.0]] # [D,L] - L cannot be zero as this leads to divide-by-zero later on
        finL_D_L_ls = [[0.0,1.0]] # [D,L] - L cannot be zero as this leads to divide-by-zero later on

        if window_shading:
            for shade_obj in window_shading:
                if shade_obj["type"] == "overhang":
                    ovh_D_L_ls.append([shade_obj["depth"],shade_obj["distance"]])
                elif shade_obj["type"] == "sidefinright":
                    finR_D_L_ls.append([shade_obj["depth"],shade_obj["distance"]])
                elif shade_obj["type"] == "sidefinleft":
                    finL_D_L_ls.append([shade_obj["depth"],shade_obj["distance"]])
                else:
                    sys.exit("shading object type" + shade_obj["type"] + "not recognised")
        
        #the default values should not be used if shading is specified
        if len(ovh_D_L_ls) >= 2:
            ovh_D_L_ls.pop(0)
        if len(finR_D_L_ls) >= 2:
            finR_D_L_ls.pop(0)
        if len(finL_D_L_ls) >= 2:
            finL_D_L_ls.pop(0)

        #perform the diff shading calculation for each comination of overhangs and fins
        for ovh_D_L,finR_D_L,finL_D_L in product(ovh_D_L_ls,finR_D_L_ls,finL_D_L_ls):
            D_ovh = ovh_D_L[0]
            L_ovh = ovh_D_L[1]
            D_finL = finL_D_L[0]
            L_finL = finL_D_L[1]
            D_finR = finR_D_L[0]
            L_finR = finR_D_L[1]
            # Calculate required geometric ratios
            # Note: PD CEN ISO/TR 52016-2:2017 Section F.6.3 refers to ISO 52016-1:2017
            #       Section F.5.5.1.6 for the definition of P1 and P2. However, this
            #       section does not exist. Therefore, these definitions have been
            #       taken from Section F.3.5.1.2 instead, also supported by Table F.6
            #       in PD CEN ISO/TR 52016-2:2017. These sources define P1 and P2
            #       differently for fins and for overhangs so it is assumed that
            #       should also apply here.
            P1_ovh = D_ovh / height
            P2_ovh = L_ovh / height
            P1_finL = D_finL / width
            P2_finL = L_finL / width
            P1_finR = D_finR / width
            P2_finR = L_finR / width

            # Calculate view factors (eqns F.15 to F.18) required for eqns F.9 to F.14
            # Note: The equations in the standard refer to P1 and P2, but as per the
            #       comment above, there are different definitions of these for fins
            #       and for overhangs. The decision on which ones to use for each of
            #       the equations below has been made depending on which of the
            #       subsequent equations the resulting variables are used in (e.g.
            #       F_w_s is used to calculate F_sh_dif_fins so we use P1 and P2 for
            #       fins).
            # Note: For F_w_r, we could set P1 equal to P1 for fins and P2 equal to
            #       P1 (not P2) for overhangs, as this appears to be consistent with
            #       example in Table F.6
            # F_w_r = 1 - exp(-0.8632 * (P1_fin + P1_ovh))
            # Note: Formula in standard for view factor to fins seems to assume that
            #       fins are the same on each side. Therefore, here we take the
            #       average of this view factor calculated with the dimensions of
            #       each fin.
            F_w_s \
                = ( 0.6514 * (1 - (P2_finL / sqrt(P1_finL * P1_finL + P2_finL * P2_finL)))
                  + 0.6514 * (1 - (P2_finR / sqrt(P1_finR * P1_finR + P2_finR * P2_finR)))
                  ) \
                / 2
            F_w_o = 0.3282 * (1 - (P2_ovh / sqrt(P1_ovh * P1_ovh + P2_ovh * P2_ovh)))
            F_w_sky = (1 - sin(alpha + beta - radians(90))) / 2

            # Calculate denominators of eqns F.9 to F.14
            view_factor_sky_no_obstacles = (1 + cos(beta)) / 2
            view_factor_ground_no_obstacles = (1 - cos(beta)) / 2

            # Setback and remote obstacles (eqns F.9 and F.10): Top half of each eqn
            # is view factor to sky (F.9) or ground (F.10) with setback and distant
            # obstacles
            # TODO Uncomment these lines when definitions of P1 and P2 in formula
            #      for F_w_r have been confirmed.
            # if view_factor_sky_no_obstacles == 0:
            #     # Shading makes no difference if sky not visible (avoid divide-by-zero)
            #     F_sh_dif_setback = 1.0
            # else:
            #     F_sh_dif_setback = (1 - F_w_r) * F_w_sky \
            #                      / view_factor_sky_no_obstacles
            # if view_factor_ground_no_obstacles == 0:
            #     # Shading makes no difference if ground not visible (avoid divide-by-zero)
            #     F_sh_ref_setback = 1.0
            # else:
            #     F_sh_ref_setback = (1 - F_w_r) * (1 - F_w_sky) \
            #                      / view_factor_ground_no_obstacles

            # Fins and remote obstacles (eqns F.11 and F.12): Top half of each eqn
            # is view factor to sky (F.11) or ground (F.12) with fins and distant
            # obstacles
            if view_factor_sky_no_obstacles == 0:
                # Shading makes no difference if sky not visible (avoid divide-by-zero)
                F_sh_dif_fins = 1.0
            else:
                F_sh_dif_fins = (1 - F_w_s) * F_w_sky \
                              / view_factor_sky_no_obstacles
            if view_factor_ground_no_obstacles == 0:
                # Shading makes no difference if ground not visible (avoid divide-by-zero)
                F_sh_ref_fins = 1.0
            else:
                F_sh_ref_fins = (1 - F_w_s) * (1 - F_w_sky) \
                              / view_factor_ground_no_obstacles

            # Overhangs and remote obstacles (eqns F.13 and F.14)
            # Top half of eqn F.13 is view factor to sky with overhangs
            if view_factor_sky_no_obstacles == 0:
                # Shading makes no difference if sky not visible (avoid divide-by-zero)
                F_sh_dif_overhangs = 1.0
            else:
                F_sh_dif_overhangs = (F_w_sky - F_w_o) \
                                   / view_factor_sky_no_obstacles
            # Top half of eqn F.14 is view factor to ground with distant obstacles,
            # but does not account for overhangs blocking any part of the view of
            # the ground, presumably because this will not happen in the vast
            # majority of cases
            if view_factor_ground_no_obstacles == 0:
                # Shading makes no difference if ground not visible (avoid divide-by-zero)
                F_sh_ref_overhangs = 1.0
            else:
                F_sh_ref_overhangs = (1 - F_w_sky) \
                                   / view_factor_ground_no_obstacles

            # Keep the smallest of the three shading reduction factors as the
            # diffuse or reflected shading factor. Also enforce that these cannot be
            # negative (which may happen with some extreme tilt values)
            # TODO Add setback shading factors to the arguments to min function when
            #      definitions of P1 and P2 in formula for F_w_r have been confirmed.
            # F_sh_dif = max(0.0, min(F_sh_dif_setback, F_sh_dif_fins, F_sh_dif_overhangs))
            # F_sh_ref = max(0.0, min(F_sh_ref_setback, F_sh_ref_fins, F_sh_ref_overhangs))
            F_sh_dif = max(0.0, min(F_sh_dif_fins, F_sh_dif_overhangs))
            F_sh_ref = max(0.0, min(F_sh_ref_fins, F_sh_ref_overhangs))

            if diffuse_irr_total == 0:
                sys.exit('Error: Zero diffuse radiation with non-zero direct radiation.')
            Fdiff = ( F_sh_dif * (diffuse_irr_sky + diffuse_irr_hor)
                    + F_sh_ref * diffuse_irr_ref
                    ) \
                  / diffuse_irr_total
            Fdiff_list.append(Fdiff)

        Fdiff = min(Fdiff_list)

        return Fdiff

    def shading_reduction_factor_direct_diffuse(
            self,
            base_height,
            height,
            width,
            tilt,
            orientation,
            window_shading,
            ):
        """ calculates the direct and diffuse shading factors due to external
        shading objects

        Arguments:
        height         -- is the height of the shaded surface (if surface is tilted then
                          this must be the vertical projection of the height), in m
        base_height    -- is the base height of the shaded surface k, in m
        width          -- is the width of the shaded surface, in m
        orientation    -- is the orientation angle of the inclined surface, expressed as the 
                          geographical azimuth angle of the horizontal projection of the 
                          inclined surface normal, -180 to 180, in degrees;
        tilt           -- is the tilt angle of the inclined surface from horizontal, measured 
                          upwards facing, 0 to 180, in degrees;
        window_shading -- data on overhangs and side fins associated to this building element
                          includes the shading object type, depth, anf distance from element
        """

        # first chceck if there is any radiation. This is needed to prevent a potential 
        # divide by zero error in the final step, but also, if there is no radiation 
        # then shading is irrelevant and we can skip the whole calculation
        direct, diffuse, _, diffuse_breakdown \
            = self.calculated_direct_diffuse_total_irradiance(tilt, orientation, True)
        if direct + diffuse == 0:
            return 0.0, 0.0

        window_shading_expanded = []
        if window_shading:
            for shading_obj in window_shading:
                if shading_obj["type"] == "reveal":
                    window_shading_expanded.append(
                        {
                            "type": "overhang",
                            "depth": shading_obj["depth"],
                            "distance": shading_obj["distance"],
                        }
                    )
                    window_shading_expanded.append(
                        {
                            "type": "sidefinleft",
                            "depth": shading_obj["depth"],
                            "distance": shading_obj["distance"]
                        }
                    )
                    window_shading_expanded.append(
                        {
                            "type": "sidefinright",
                            "depth": shading_obj["depth"],
                            "distance": shading_obj["distance"]
                        }
                    )
                else:
                    window_shading_expanded.append(shading_obj)
        # first check if the surface is outside the solar beam
        # if so then direct shading is complete and we don't need to
        # calculate shading from objects
        # TODO The outside solar beam condition is based on a vertical projection
        #      of the surface and does not account for the condition where a
        #      surface that is only slightly pitched is exposed to direct solar
        #      radiation when the sun is high (e.g. a surface pitched slightly
        #      to the north will be exposed to direct solar radiation when the
        #      sun is high in the southern sky). As the solar radiation
        #      calculation already accounts for the situation where the sun is
        #      actually behind the surface (accounting for the combination of
        #      orientation and pitch), there is no need to zero it using the
        #      shading factor. For now, we set the shading factor to 1 and ignore
        #      shading from objects on the other side of the building (which if
        #      significantly pitched would have to be relatively tall and/or
        #      very close to cast a shadow on the surface in question anyway),
        #      so that results in the unshaded case will be correct.
        if self.outside_solar_beam(tilt, orientation):
            Fdir = 1.0
        else:
            Fdir = self.direct_shading_reduction_factor \
                    (base_height, height, width, orientation, window_shading_expanded)

        Fdiff = self.diffuse_shading_reduction_factor(
            diffuse_breakdown,
            tilt,
            height,
            width,
            window_shading_expanded,
            )

        return Fdir, Fdiff

        # TODO suspected bug identified in ISO 52016 as it conflicts with ISO 52010:
        #ISO 52010 states that (6.4.5.2.1) the total irradiance on the inclined surface is
        #Itotal = Fdir * Idirect + Idiffuse
        #This is how the shading factor is used with the solar gains calculation.
        #However, ISO 52016 takes Fdir and performs the calculation below to give a "final"
        #shading factor. This does not make sense to be applied solely to the direct radiation
        #when calculating solar gains. Therefore we return Fdir here.

        #Fshade = (Fdir * direct + diffuse) / (direct + diffuse)

        #return Fshade
