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
from math import log, exp
import numpy as np

from .unitless import speed_of_sound, dynamic_pressure_p_mach
from .unit_conversion import  (
    convert_altitude, convert_density, convert_velocity,
    _rankine_to_temperature_units, _psfs_to_dvisc_units, _ft2s_to_kvisc_units,
    _altitude_factor, _pressure_factor, _velocity_factor,
    _reynolds_factor)


def atm_temperature(alt, alt_units='ft', temperature_units='R'):
    # type : (float, str, str) -> float
    r"""
    Freestream Temperature \f$ T_{\infty} \f$

    Parameters
    ----------
    alt : float
        Altitude in alt_units
    alt_units : str; default='ft'
        the altitude units; ft, kft, m
    temperature_units : str; default='R'
        the altitude units; R, K

    Returns
    -------
    T : float
        temperature in degrees Rankine or Kelvin (SI)

    .. note ::
        from BAC-7006-3352-001-V1.pdf\n
        A Manual for Determining Aerodynamic Heating of High Speed Aircraft\n
        page ~236 - Table C.1\n
        These equations were used because they are valid to 300k ft.\n
        Extrapolation is performed above that.
    """
    z = alt * _altitude_factor(alt_units, 'ft')
    if z < 36151.725:
        T = 518.0 - 0.003559996 * z
    elif z < 82344.678:
        T = 389.988
    elif z < 155347.756:
        T = 389.988 + .0016273286 * (z - 82344.678)
    elif z < 175346.171:
        T = 508.788
    elif z < 249000.304:
        T = 508.788 - .0020968273 * (z - 175346.171)
    elif z < 299515.564:
        T = 354.348
    else:
        #print("alt=%i kft > 299.5 kft" % (z / 1000.))
        T = 354.348
        #raise AtmosphereError("altitude is too high")

    factor = _rankine_to_temperature_units(temperature_units)
    T2 = T * factor
    return T2


def atm_pressure(alt, alt_units='ft', pressure_units='psf'):
    # type : (float, str, str) -> float
    r"""
    Freestream Pressure \f$ p_{\infty} \f$

    Parameters
    ----------
    alt : float
        Altitude in alt_units
    alt_units : str; default='ft'
        the altitude units; ft, kft, m
    pressure_units : str; default='psf'
        the pressure units; psf, psi, Pa, kPa, MPa

    Returns
    -------
    pressure : float
        Returns pressure in pressure_units

    .. note ::
        from BAC-7006-3352-001-V1.pdf\n
        A Manual for Determining Aerodynamic Heating of High Speed Aircraft\n
        page ~236 - Table C.1\n
        These equations were used b/c they are valid to 300k ft.\n
        Extrapolation is performed above that.\n
    """
    alt_ft = convert_altitude(alt, alt_units, 'ft')
    ln_pressure = _log_pressure(alt_ft)
    p = exp(ln_pressure)

    factor = _pressure_factor('psf', pressure_units)
    return p * factor


def _log_pressure(alt_ft):
    """calculates the log(pressure) in psf given altitude in feet"""
    if alt_ft < 36151.725:
        ln_pressure = 7.657389 + 5.2561258 * log(1 - 6.8634634E-6 * alt_ft)
    elif alt_ft < 82344.678:
        ln_pressure = 6.158411 - 4.77916918E-5 * (alt_ft - 36151.725)
    elif alt_ft < 155347.756:
        ln_pressure = 3.950775 - 11.3882724 * log(1.0 + 4.17276598E-6 * (alt_ft - 82344.678))
    elif alt_ft < 175346.171:
        ln_pressure = 0.922461 - 3.62635373E-5*(alt_ft - 155347.756)
    elif alt_ft < 249000.304:
        ln_pressure = 0.197235 + 8.7602095 * log(1.0 - 4.12122002E-6 * (alt_ft - 175346.171))
    elif alt_ft < 299515.564:
        ln_pressure = -2.971785 - 5.1533546650E-5 * (alt_ft - 249000.304)
    else:
        #print("alt=%i kft > 299.5 kft" % (alt_ft / 1000.))
        ln_pressure = -2.971785 - 5.1533546650E-5 * (alt_ft - 249000.304)
    return ln_pressure


def atm_dynamic_pressure(alt, mach, alt_units='ft', pressure_units='psf'):
    # type : (float, float, str, str) -> float
    r"""
    Freestream Dynamic Pressure  \f$ q_{\infty} \f$

    Parameters
    ----------
    alt : float
        Altitude in alt_units
    mach : float
        Mach Number \f$ M \f$
    alt_units : str; default='ft'
        the altitude units; ft, kft, m
    pressure_units : str; default='psf'
        the pressure units; psf, psi, Pa, kPa, MPa

    Returns
    -------
    dynamic_pressure : float
        Returns dynamic pressure in pressure_units

    The common method that requires many calculations...
    \f[  \large q = \frac{1}{2} \rho V^2  \f]
    \f[  \large p = \rho R T  \f]
    \f[  \large M = \frac{V}{a}  \f]
    \f[  \large a = \sqrt{\gamma R T}  \f]
    so...
    \f[  \large q = \frac{\gamma}{2} p M^2  \f]
    """
    z = alt * _altitude_factor(alt_units, 'ft')
    p = atm_pressure(z)
    q = dynamic_pressure_p_mach(p, mach)

    factor = _pressure_factor('psf', pressure_units)
    q2 = q * factor
    return q2


def atm_speed_of_sound(alt, alt_units='ft', velocity_units='ft/s', gamma=1.4):
    # type : (float, str, str, float) -> float
    r"""
    Freestream Speed of Sound  \f$ a_{\infty} \f$

    Parameters
    ----------
    alt : float
        Altitude in alt_units
    alt_units : str; default='ft'
        the altitude units; ft, kft, m
    velocity_units : str; default='ft/s'
        the velocity units; ft/s, m/s, in/s, knots

    Returns
    -------
    speed_of_sound, a : float
        Returns speed of sound in velocity_units

    \f[  \large a = \sqrt{\gamma R T}  \f]

    """
    # converts everything to English units first
    z = alt * _altitude_factor(alt_units, 'ft')
    T = atm_temperature(z)
    a = speed_of_sound(T, R=1716., gamma=gamma)

    factor = _velocity_factor('ft/s', velocity_units) # ft/s to m/s
    a2 = a * factor
    return a2


def atm_velocity(alt, mach, alt_units='ft', velocity_units='ft/s'):
    # type : (float, float, str, str) -> float
    r"""
    Freestream Velocity  \f$ V_{\infty} \f$

    Parameters
    ----------
    alt : float
        altitude in alt_units
    Mach : float
        Mach Number \f$ M \f$
    alt_units : str; default='ft'
        the altitude units; ft, kft, m
    velocity_units : str; default='ft/s'
        the velocity units; ft/s, m/s, in/s, knots

    Returns
    -------
    velocity : float
        Returns velocity in velocity_units

    \f[ \large V = M a \f]
    """
    a = atm_speed_of_sound(alt, alt_units=alt_units, velocity_units=velocity_units)
    V = mach * a # units=ft/s or m/s
    return V


def atm_equivalent_airspeed(alt, mach, alt_units='ft', eas_units='ft/s'):
    # type : (float, float, str, str) -> float
    """
    Freestream equivalent airspeed

    Parameters
    ----------
    alt : float
        altitude in alt_units
    Mach : float
        Mach Number \f$ M \f$
    alt_units : str; default='ft'
        the altitude units; ft, kft, m
    eas_units : str; default='ft/s'
        the equivalent airspeed units; ft/s, m/s, in/s, knots

    Returns
    -------
    eas : float
        equivalent airspeed in eas_units

    EAS = TAS * sqrt(rho/rho0)
    p = rho * R * T
    rho = p/(RT)
    rho/rho0 = p/T * T0/p0
    TAS = a * M
    EAS = a * M * sqrt(p/T * T0/p0)
    EAS = a * M * sqrt(p*T0 / (T*p0))

    """
    z = convert_altitude(alt, alt_units, 'ft')
    p_psf = atm_pressure(z)
    eas_fts = _equivalent_airspeed(mach, p_psf)
    eas2 = convert_velocity(eas_fts, 'ft/s', eas_units)
    return eas2


def _equivalent_airspeed(mach, p_psf):
    """helper method for atm_equivalent_airspeed"""
    z0 = 0.
    T0 = atm_temperature(z0)
    p0 = atm_pressure(z0)
    gamma = 1.4
    R = 1716.
    #eas = a * mach * sqrt((p * T0) / (T * p0))
    #    = sqrt(gamma * R * T) * mach * sqrt(T0 / p0) * sqrt(p / T)
    #    = sqrt(gamma * R) * mach * sqrt(T0 / p0) * sqrt(T) * sqrt(p / T)
    #    = sqrt(gamma * R * T0 / p0) * mach * sqrt(p)
    #    = k * sqrt(p)
    # rho0 = p0 / (R * T0)
    # k = sqrt(gamma / rho0) * mach
    eas = np.sqrt(gamma * R * T0 / p0) * mach * p_psf ** 0.5
    return eas


def atm_mach(alt, V, alt_units='ft', velocity_units='ft/s'):
    # type : (float, float, str, str) -> float
    r"""
    Freestream Mach Number

    Parameters
    ----------
    alt : float
        altitude in alt_units
    V : float
        Velocity in velocity_units
    alt_units : str; default='ft'
        the altitude units; ft, kft, m
    velocity_units : str; default='ft/s'
        the velocity units; ft/s, m/s, in/s, knots

    Returns
    -------
    mach : float
        Mach Number \f$ M \f$

    \f[ \large M = \frac{V}{a} \f]
    """
    a = atm_speed_of_sound(alt, alt_units=alt_units, velocity_units=velocity_units)
    mach = V / a
    return mach


def atm_density(alt, R=1716., alt_units='ft', density_units='slug/ft^3'):
    # type : (float, float, str, str) -> float
    r"""
    Freestream Density   \f$ \rho_{\infty} \f$

    Parameters
    ----------
    alt : float
        altitude in feet or meters
    R : float; default=1716.
        gas constant for air in english units (???)
    alt_units : str; default='ft'
        the altitude units; ft, kft, m
    density_units : str; default='slug/ft^3'
        the density units; slug/ft^3, slinch/in^3, kg/m^3

    Returns
    -------
    rho : float
        density \f$ \rho \f$ in density_units

    Based on the formula P=pRT
    \f[ \large \rho=\frac{p}{R T} \f]
    """
    z = convert_altitude(alt, alt_units, 'ft')
    #z = alt * _altitude_factor(alt_units, 'ft')
    p = atm_pressure(z)
    T = atm_temperature(z)

    rho = p / (R * T)
    rho2 = convert_density(rho, 'slug/ft^3', density_units)
    return rho2


def atm_kinematic_viscosity_nu(alt, alt_units='ft', visc_units='ft^2/s'):
    # type : (float, str, str) -> float
    r"""
    Freestream Kinematic Viscosity \f$ \nu_{\infty} \f$

    Parameters
    ----------
    alt : float
        Altitude in alt_units
    alt_units : str; default='ft'
        the altitude units; ft, kft, m
    visc_units : str; default='slug/ft^3'
        the kinematic viscosity units; ft^2/s, m^2/s

    Returns
    -------
    nu : float
        kinematic viscosity \f$ \nu_{\infty} \f$ in visc_units

    \f[ \large \nu = \frac{\mu}{\rho} \f]

    .. seealso::  sutherland_viscoscity
    .. todo:: better debug

    """
    z = alt * _altitude_factor(alt_units, 'ft')
    rho = atm_density(z)
    mu = atm_dynamic_viscosity_mu(z)
    nu = mu / rho  # ft^2/s
    factor = _ft2s_to_kvisc_units(alt_units, visc_units)
    return nu * factor


def atm_dynamic_viscosity_mu(alt, alt_units='ft', visc_units='(lbf*s)/ft^2'):
    # type : (float, str, str) -> float
    r"""
    Freestream Dynamic Viscosity  \f$ \mu_{\infty} \f$

    Parameters
    ----------
    alt : float
        Altitude in alt_units
    alt_units : str; default='ft'
        the altitude units; ft, kft, m
    visc_units : str; default='(lbf*s)/ft^2'
        the viscosity units; (lbf*s)/ft^2, (N*s)/m^2, Pa*s

    Returns
    -------
    mu : float
        dynamic viscosity  \f$ \mu_{\infty} \f$ in (lbf*s)/ft^2 or (N*s)/m^2 (SI)

    .. seealso::  sutherland_viscoscity

    """
    z = alt * _altitude_factor(alt_units, 'ft')
    T = atm_temperature(z)
    mu = sutherland_viscoscity(T)  # (lbf*s)/ft^2
    factor = _psfs_to_dvisc_units(visc_units)
    return mu * factor


def atm_unit_reynolds_number2(alt, mach, alt_units='ft', reynolds_units='1/ft'):
    # type : (float, float, str, str) -> float
    r"""
    Returns the Reynolds Number per unit length.

    Parameters
    ----------
    alt : float
        Altitude in alt_units
    mach : float
        Mach Number \f$ M \f$
    alt_units : str; default='ft'
        the altitude units; ft, kft, m
    reynolds_units : str; default='1/ft'
        the altitude units; 1/ft, 1/m, 1/in

    Returns
    -------
    ReynoldsNumber/L : float
        the Reynolds Number per unit length

    \f[ \large Re_L = \frac{ \rho V}{\mu} = \frac{p M a}{\mu R T} \f]

    .. note ::
        this version of Reynolds number directly caculates the base quantities, so multiple
        calls to atm_press and atm_temp are not made
    """
    z = alt * _altitude_factor(alt_units, 'ft')
    gamma = 1.4
    R = 1716.
    p = atm_pressure(z)
    T = atm_temperature(z)
    mu = sutherland_viscoscity(T)
    #p = rho * R * T
    #a = (gamma * R * T) ** 0.5
    #
    # ReL = p * a * mach / (mu * R * T)
    #     = p * sqrt(gamma * R * T) * mach / (mu * R * T)
    #     = (p * mach / mu) * sqrt(gamma * R * T) / (R * T)
    #     = (p * mach / mu) * sqrt(gamma / (R * T))
    ReL = (p * mach / mu) * (gamma / (R * T)) ** 0.5
    ReL *= _reynolds_factor('1/ft', reynolds_units)
    return ReL


def atm_unit_reynolds_number(alt, mach, alt_units='ft', reynolds_units='1/ft'):
    # type : (float, float, str, str) -> float
    r"""
    Returns the Reynolds Number per unit length.

    Parameters
    ----------
    alt : float
        Altitude in alt_units
    mach : float
        Mach Number \f$ M \f$
    alt_units : str; default='ft'
        the altitude units; ft, kft, m
    reynolds_units : str; default='1/ft'
        the altitude units; 1/ft, 1/m, 1/in

    Returns
    -------
    ReynoldsNumber/L : float
        Reynolds number per unit length in reynolds_units

    \f[ \large Re   = \frac{ \rho V L}{\mu} \f]
    \f[ \large Re_L = \frac{ \rho V  }{\mu} \f]
    """
    z = alt * _altitude_factor(alt_units, 'ft')
    rho = atm_density(z)
    V = atm_velocity(z, mach)
    mu = atm_dynamic_viscosity_mu(z)

    ReL = (rho * V) / mu
    ReL *= _reynolds_factor('1/ft', reynolds_units)
    return ReL


def sutherland_viscoscity(T):
    # type: (float) -> float
    r"""
    Helper function that calculates the dynamic viscosity \f$ \mu \f$ of air at
    a given temperature.

    Parameters
    ----------
    T : float
        Temperature T is in Rankine

    Returns
    -------
    mu : float
       dynamic viscosity  \f$ \mu \f$ of air in (lbf*s)/ft^2

    .. note ::
        prints a warning if T>5400 deg R

    Sutherland's Equation\n
    From Aerodynamics for Engineers 4th Edition\n
    John J. Bertin 2002\n
    page 6 eq 1.5b\n
    """
    if T < 225.: # Rankine
        viscosity = 8.0382436E-10 * T
    else:
        if T > 5400.:
            sys.stderr.write('WARNING:  viscosity - Temperature is too large '
                             '(T>5400 R) T=%s\n' % T)
        viscosity = 2.27E-8 * (T ** 1.5) / (T + 198.6)
    return viscosity
