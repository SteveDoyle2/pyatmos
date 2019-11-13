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
import sys
from typing import Union
import numpy as np

from .atmosphere import atm_temperature, _log_pressure, _equivalent_airspeed, atm_pressure
from .unitless import speed_of_sound
from .unit_conversion import (
    convert_velocity, convert_density,
    _rankine_to_temperature_units, _psfs_to_dvisc_units, _ft2s_to_kvisc_units,
    _altitude_factor, _pressure_factor, _velocity_factor,
)

def atm_dynamic_pressure_array(alt: np.array, mach: np.array,
                               alt_units: str='ft', pressure_units: str='psf') -> np.array:
    r"""
    Freestream Dynamic Pressure  \f$ q_{\infty} \f$

    Parameters
    ----------
    alt : np.ndarray
        Altitude in alt_units
    mach : np.ndarray
        Mach Number \f$ M \f$
    alt_units : str; default='ft'
        the altitude units; ft, kft, m
    pressure_units : str; default='psf'
        the pressure units; psf, psi, Pa, kPa, MPa

    Returns
    -------
    dynamic_pressure : np.ndarray
        Returns dynamic pressure in pressure_units

    The common method that requires many calculations...
    \f[  \large q = \frac{1}{2} \rho V^2  \f]
    \f[  \large p = \rho R T  \f]
    \f[  \large M = \frac{V}{a}  \f]
    \f[  \large a = \sqrt{\gamma R T}  \f]
    so...
    \f[  \large q = \frac{\gamma}{2} p M^2  \f]
    """
    alt_shape = alt.shape
    z = alt * _altitude_factor(alt_units, 'ft')
    p = atm_pressure_array(z)
    q = 0.7 * p * mach ** 2

    factor = _pressure_factor('psf', pressure_units)
    q2 = q * factor
    return q2

def atm_pressure_array(alt: np.array,
                       alt_units: str='ft', pressure_units: str='psf') -> np.array:
    """Gets the pressure as a numpy array"""
    alt_ft = alt * _altitude_factor(alt_units, 'ft')
    ln_pressure = np.array([_log_pressure(alti) for alti in alt_ft.ravel()])
    press_psf = np.exp(ln_pressure).reshape(alt.shape)
    return press_psf * _pressure_factor('psf', pressure_units)


def atm_temperature_array(alt: np.array,
                          alt_units: str='ft', temperature_units: str='R') -> np.array:
    """Gets the temperature as a numpy array"""
    temp_rankine = np.array([atm_temperature(alti, alt_units=alt_units, temperature_units='R')
                             for alti in alt.ravel()]).reshape(alt.shape)
    return temp_rankine * _rankine_to_temperature_units(temperature_units)


def atm_speed_of_sound_array(alt: np.array,
                             alt_units: str='ft', velocity_units: str='ft/s', gamma: float=1.4) -> np.array:
    """Gets the speed of sound as a numpy array"""
    T = atm_temperature_array(alt, alt_units=alt_units, temperature_units='R')
    a = speed_of_sound(T, R=1716., gamma=gamma)

    factor = _velocity_factor('ft/s', velocity_units) # ft/s to m/s
    a2 = a * factor
    return a2


def atm_density_array(alt: np.array, R: float=1716.,
                      alt_units: str='ft', density_units: str='slug/ft^3') -> np.array:
    """Gets the density as a numpy array"""
    alt = np.asarray(alt)
    alt_ft = alt * _altitude_factor(alt_units, 'ft')
    p = atm_pressure_array(alt_ft)
    T = atm_temperature_array(alt_ft)

    rho = p / (R * T)
    rho2 = convert_density(rho, 'slug/ft^3', density_units)
    return rho2


def atm_dynamic_viscosity_mu_array(alt: np.array,
                                   alt_units: str='ft', visc_units: str='(lbf*s)/ft^2') -> np.array:
    """Gets the dynamic viscosity as a numpy array"""
    alt = np.asarray(alt)
    alt_ft = alt * _altitude_factor(alt_units, 'ft')
    T = atm_temperature_array(alt_ft)  # R
    mu = sutherland_viscoscity_array(T)  # (lbf*s)/ft^2
    factor = _psfs_to_dvisc_units(visc_units)
    return mu * factor

def atm_kinematic_viscosity_nu_array(alt: np.array,
                                     alt_units: str='ft', visc_units: str='ft^2/s') -> np.array:
    """Gets the kinematic viscosity as a numpy array"""
    alt = np.asarray(alt)
    alt_ft = alt * _altitude_factor(alt_units, 'ft')
    rho = atm_density_array(alt_ft)
    mu = atm_dynamic_viscosity_mu_array(alt_ft)
    nu = mu / rho  # ft^2/s
    factor = _ft2s_to_kvisc_units(alt_units, visc_units)
    return nu * factor


def sutherland_viscoscity_array(T: np.array) -> np.array:
    """Gets the Sutherland viscosity as a numpy array"""
    viscosity = 2.27E-8 * (T ** 1.5) / (T + 198.6)
    ilow = np.where(T < 225.)[0]
    if len(ilow):
        viscosity[ilow] = 8.0382436E-10 * T[ilow]

    ihigh = np.where(T > 5400)[0]
    if len(ihigh):
        sys.stderr.write('WARNING:  viscosity - Temperature is too large '
                         '(T>5400 R) Tmax=%s\n' % T[ihigh].max())
    return viscosity


def atm_equivalent_airspeed_array(alt: Union[float, np.array], mach: np.array,
                                  alt_units: str='ft', eas_units: str='ft/s') -> np.array:
    """Gets the equivalent airspeed as a numpy array"""
    if isinstance(alt, float):
        pressi = atm_pressure(alt, alt_units=alt_units)
        p_psf = pressi * np.ones(len(mach), dtype=mach.dtype)
    else:
        alt = np.asarray(alt)
        alt_ft = alt * _altitude_factor(alt_units, 'ft')
        p_psf = atm_pressure_array(alt_ft)
    eas_fts = _equivalent_airspeed(mach, p_psf)
    eas2 = convert_velocity(eas_fts, 'ft/s', eas_units)
    return eas2
