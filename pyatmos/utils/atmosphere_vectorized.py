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

from .atmosphere import atm_temperature, atm_pressure, _log_pressure
from .unitless import speed_of_sound
from .unit_conversion import (
    convert_altitude, convert_velocity, convert_density,
    _altitude_factor, _temperature_factor, _pressure_factor, _velocity_factor,
    _feet_to_alt_units)


def atm_pressure_array(alt, alt_units='ft', pressure_units='psf'):
    """Gets the pressure as a numpy array"""
    alt_ft = alt * _altitude_factor(alt_units, 'ft')
    ln_pressure = np.array([_log_pressure(alti) for alti in alt_ft])
    press_psf = np.exp(ln_pressure)
    return press_psf * _pressure_factor('psf', pressure_units)


def atm_temperature_array(alt, alt_units='ft', temperature_units='R'):
    """Gets the temperature as a numpy array"""
    temp_rankine = np.array([atm_temperature(alti, alt_units=alt_units, temperature_units='R')
                             for alti in alt])
    return temp_rankine * _temperature_factor('R', temperature_units)


def atm_speed_of_sound_array(alt, alt_units='ft', velocity_units='ft/s', gamma=1.4):
    """Gets the speed of sound as a numpy array"""
    T = atm_temperature_array(alt, alt_units=alt_units, temperature_units='R')
    a = speed_of_sound(T, R=1716., gamma=gamma)

    factor = _velocity_factor('ft/s', velocity_units) # ft/s to m/s
    a2 = a * factor
    return a2


def atm_density_array(alt, R=1716., alt_units='ft', density_units='slug/ft^3'):
    """Gets the density as a numpy array"""
    z = convert_altitude(alt, alt_units, 'ft')
    p = atm_pressure_array(z)
    T = atm_temperature_array(z)

    rho = p / (R * T)
    rho2 = convert_density(rho, 'slug/ft^3', density_units)
    return rho2


def atm_dynamic_viscosity_mu_array(alt, alt_units='ft', visc_units='(lbf*s)/ft^2'):
    # type : (Any, str, str) -> Any
    """Gets the dynamic viscosity as a numpy array"""
    z = alt * _altitude_factor(alt_units, 'ft')
    T = atm_temperature_array(z)
    mu = sutherland_viscoscity_array(T)  # (lbf*s)/ft^2

    # same units as pressure, except multiplied by seconds
    if visc_units == '(lbf*s)/ft^2':
        factor = 1.
    elif visc_units in ['(N*s)/m^2', 'Pa*s']:
        factor = 47.88026
    else:
        raise NotImplementedError('visc_units=%r; not in (lbf*s)/ft^2 or (N*s)/m^2 or Pa*s')
    return mu * factor


def atm_kinematic_viscosity_nu_array(alt, alt_units='ft', visc_units='ft^2/s'):
    # type : (Any, str, str) -> Any
    """Gets the kinematic viscosity as a numpy array"""
    z = alt * _altitude_factor(alt_units, 'ft')
    rho = atm_density_array(z)
    mu = atm_dynamic_viscosity_mu_array(z)
    nu = mu / rho

    if visc_units == 'ft^2/s':
        factor = 1.
    elif visc_units == 'm^2/s':
        factor = _feet_to_alt_units(alt_units) ** 2
    else:
        raise NotImplementedError('visc_units=%r' % visc_units)
    return nu * factor


def sutherland_viscoscity_array(T):
    """Gets the Sutherland viscosity as a numpy array"""
    viscosity = 2.27E-8 * (T ** 1.5) / (T + 198.6)
    ilow = np.where(T < 225.)[0]
    if len(ilow):
        viscosity[ilow] = 8.0382436E-10 * T

    ihigh = np.where(T > 5400)[0]
    if len(ihigh):
        sys.stderr.write('WARNING:  viscosity - Temperature is too large '
                         '(T>5400 R) Tmax=%s\n' % T[ihigh].max())
    return viscosity


def atm_equivalent_airspeed_array(alt, mach, alt_units='ft', eas_units='ft/s'):
    # type : (Any, Any, str, str) -> Any
    """Gets the equivalent airspeed as a numpy array"""
    z = convert_altitude(alt, alt_units, 'ft')
    z0 = 0.
    T0 = atm_temperature(z0)
    p0 = atm_pressure(z0)
    p = atm_pressure_array(z)

    gamma = 1.4
    R = 1716.
    #eas = a * mach * sqrt((p * T0) / (T * p0))
    #    = sqrt(gamma * R * T) * mach * sqrt(T0 / p0) * sqrt(p / T)
    #    = sqrt(gamma * R) * mach * sqrt(T0 / p0) * sqrt(T) * sqrt(p / T)
    #    = sqrt(gamma * R * T0 / p0) * mach * sqrt(p)
    #    = k * sqrt(p)
    # rho0 = p0 / (R * T0)
    # k = sqrt(gamma / rho0) * mach
    eas = np.sqrt(gamma * R * T0 / p0) * mach * p ** 0.5
    eas2 = convert_velocity(eas, 'ft/s', eas_units)
    return eas2
