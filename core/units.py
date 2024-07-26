#!/usr/bin/env python3

"""
This module contains common unit conversions for use by other modules.
"""

J_per_kWh = 3600000
kJ_per_kWh = 3600
J_per_kJ = 1000
W_per_kW = 1000
litres_per_cubic_metre = 1000
minutes_per_hour = 60
seconds_per_minute = 60
seconds_per_hour = 3600
hours_per_day = 24
days_per_year = 365
days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
knots_per_m_per_sec = 1 / (1852 / 3600)
mm_per_m = 1000

def average_monthly_to_annual(list_monthly_averages):
    assert len(list_monthly_averages) == 12
    return sum([
        month_ave * days_in_month[month_idx]
        for month_idx, month_ave in enumerate(list_monthly_averages)
        ]) \
        / sum(days_in_month)

def Celcius2Kelvin(temp_C):
    assert temp_C >= -273.15
    return temp_C + 273.15

def Kelvin2Celcius(temp_K):
    assert temp_K >= 0
    return temp_K - 273.15

def convert_profile_to_daily(original_profile, timestep):
    """ Convert profile from per-timestep figures to daily figures """
    total_steps = len(original_profile)
    steps_per_day = int(hours_per_day / timestep)
    daily_profile = [
        sum(original_profile[i : i + steps_per_day])
        for i in range(0, total_steps, steps_per_day)
    ]
    return daily_profile
