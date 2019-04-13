"""
Contains the following atmospheric functions:

 - density = atm_density(alt, mach)
 - mach = atm_mach(alt, velocity)
 - velocity = atm_velocity(alt, mach)
 - pressure = atm_pressure(alt)
 - temperature = atm_temperature(alt)
 - sos = atm_speed_of_sound(alt)
 - mu = atm_dynamic_viscosity_mu(alt)
 - nu = atm_kinematic_viscosity_nu(alt)
 - eas = atm_equivalent_airspeed(alt, mach)

All the default units are in English units because the source equations
are in English units.
"""
from __future__ import print_function, absolute_import
import sys
import numpy as np

from .unit_conversion import _altitude_factor, _temperature_factor, _pressure_factor

def speed_of_sound(T, R=1716., gamma=1.4):
    """
    Calculates the speed of sound without units

    Parameters
    ----------
    T : float, np.ndarray
        the temperature
    R : float; default=1716.0
        1716.59, dir air, R=287.04 J/kg*K
    gamma : float; default=1.4
        the ratio of Cp/Cv
    """
    a = (gamma * R * T) ** 0.5
    return a

