"""
Contains the following atmospheric functions:

 - alt, rho, mach, velocity, eas = make_alt_sweep(
       mach, alts, eas_limit=1000.,
       alt_units='m', velocity_units='m/s', density_units='kg/m^3',
       eas_units='m/s')
 - alt, rho, mach, velocity, eas = make_mach_sweep(
       alt, machs, eas_limit=1000.,
       alt_units='m', velocity_units='m/s', density_units='kg/m^3',
       eas_units='m/s')

Sweeps use SI as the default units because English mass units are tricky.
Buyer beware!
"""
from __future__ import print_function, absolute_import
import  numpy as np
from .atmosphere import atm_density, atm_speed_of_sound
from .atmosphere_vectorized import atm_pressure_array, atm_temperature_array
from .unit_conversion import _altitude_factor, _velocity_factor, _density_factor


def make_alt_sweep(mach, alts, eas_limit=1000.,
                   alt_units='m', velocity_units='m/s', density_units='kg/m^3',
                   eas_units='m/s'):
    """
    Makes a sweep across altitude for a constant Mach number.

    Parameters
    ----------
    alt : List[float]
        Altitude in alt_units
    Mach : float
        Mach Number \f$ M \f$
    eas_limit : float
        Equivalent airspeed limiter in eas_units
    alt_units : str; default='m'
        the altitude units; ft, kft, m
    velocity_units : str; default='m/s'
        the velocity units; ft/s, m/s, in/s, knots
    density_units : str; default='kg/m^3'
        the density units; slug/ft^3, slinch/in^3, kg/m^3
    eas_units : str; default='m/s'
        the equivalent airspeed units; ft/s, m/s, in/s, knots
    """
    alt = np.asarray(alts, dtype='float64')
    alt_ft = alt * _altitude_factor(alt_units, 'ft')

    press_psf = atm_pressure_array(alt_ft, alt_units='ft', pressure_units='psf')
    temp_rankine = atm_temperature_array(alt_ft, alt_units='ft', temperature_units='R')

    R = 1716.
    gamma = 1.4

    # p=rho*R*T
    # rho=p/RT
    rho = press_psf / (R * temp_rankine) * _density_factor('slug/ft^3', density_units)
    sos_fts = np.sqrt(gamma * R * temp_rankine)
    sos = sos_fts * _velocity_factor('ft/s', velocity_units)

    mach = np.ones(len(alt)) * mach
    velocity = sos * mach
    alt, rho, mach, velocity, eas = _limit_eas(
        alt, rho, mach, velocity, eas_limit,
        alt_units=alt_units,
        density_units=density_units,
        velocity_units=velocity_units,
        eas_units=eas_units,)
    return alt, rho, mach, velocity, eas

def make_mach_sweep(alt, machs, eas_limit=1000.,
                    alt_units='m', velocity_units='m/s', density_units='kg/m^3',
                    eas_units='m/s'):
    """
    Makes a sweep across Mach number for a constant altitude.

    Parameters
    ----------
    alt : float
        Altitude in alt_units
    Machs : List[float]
        Mach Number \f$ M \f$
    eas_limit : float
        Equivalent airspeed limiter in eas_units
    alt_units : str; default='m'
        the altitude units; ft, kft, m
    velocity_units : str; default='m/s'
        the velocity units; ft/s, m/s, in/s, knots
    density_units : str; default='kg/m^3'
        the density units; slug/ft^3, slinch/in^3, kg/m^3
    eas_units : str; default='m/s'
        the equivalent airspeed units; ft/s, m/s, in/s, knots
    """
    mach = np.asarray(machs)
    rho = np.ones(len(mach)) * atm_density(alt, R=1716., alt_units=alt_units,
                                           density_units=density_units)
    sos = atm_speed_of_sound(alt, alt_units=alt_units,
                             velocity_units=velocity_units)

    alt = np.ones(len(mach)) * alt
    velocity = sos * mach
    alt, rho, mach, velocity, eas = _limit_eas(
        alt, rho, mach, velocity, eas_limit,
        alt_units=alt_units,
        density_units=density_units,
        velocity_units=velocity_units,
        eas_units=eas_units,)
    return alt, rho, mach, velocity, eas

def make_eas_sweep(alt, eass, alt_units='m', velocity_units='m/s', density_units='kg/m^3',
                   eas_units='m/s'):
    """
    Makes a sweep across equivalent airspeed for a constant altitude.

    Parameters
    ----------
    alt : float
        Altitude in alt_units
    eass : List[float]
        Equivalent airspeed in eas_units
    alt_units : str; default='m'
        the altitude units; ft, kft, m
    velocity_units : str; default='m/s'
        the velocity units; ft/s, m/s, in/s, knots
    density_units : str; default='kg/m^3'
        the density units; slug/ft^3, slinch/in^3, kg/m^3
    eas_units : str; default='m/s'
        the equivalent airspeed units; ft/s, m/s, in/s, knots
    """
    # convert eas to output units
    eas = np.atleast_1d(eass) * _velocity_factor(eas_units, velocity_units)

    rho = atm_density(alt, R=1716., alt_units=alt_units,
                      density_units=density_units)
    sos = atm_speed_of_sound(alt, alt_units=alt_units,
                             velocity_units=velocity_units)
    rho0 = atm_density(0., alt_units=alt_units, density_units=density_units)

    alt = np.ones(eass.shape, dtype='float64') * alt
    velocity = eass * np.sqrt(rho0 / rho)
    mach = velocity / sos
    return alt, rho, mach, velocity, eas

def _limit_eas(alt, rho, mach, velocity, eas_limit=1000.,
               alt_units='m', velocity_units='m/s', density_units='kg/m^3',
               eas_units='m/s'):
    """limits the equivalent airspeed"""
    eas_out = None
    if eas_limit:
        rho0 = atm_density(0., alt_units=alt_units, density_units=density_units)

        # eas in velocity units
        eas = velocity * np.sqrt(rho / rho0)
        kvel = _velocity_factor(eas_units, velocity_units)
        eas_limit_in_velocity_units = eas_limit * kvel

        i = np.where(eas < eas_limit_in_velocity_units)
        alt = alt[i]
        rho = rho[i]
        mach = mach[i]
        velocity = velocity[i]
        eas_out = eas[i]

        if len(rho) == 0:
            #print('mach min: %.0f max: %.0f' % (mach.min(), mach.max()))
            #print('vel min: %.0f max: %.0f in/s' % (velocity.min(), velocity.max()))
            #print('EAS min: %.0f max: %.0f in/s' % (eas.min(), eas.max()))
            raise RuntimeError('EAS limit is too struct and has removed all the conditions.\n'
                               'Increase eas_limit or change the mach/altude range\n'
                               '  EAS: min=%.3f max=%.3f limit=%s %s' % (
                                   eas.min() / kvel,
                                   eas.max() / kvel,
                                   eas_limit, eas_units))
    return alt, rho, mach, velocity, eas_out
