from __future__ import print_function, absolute_import

def _feet_to_alt_units(alt_units):
    # type : (str) -> float
    """helper method"""
    if alt_units == 'm':
        factor = 0.3048
    elif alt_units == 'ft':
        factor = 1.
    else:
        raise RuntimeError('alt_units=%r is not valid; use [ft, m]' % alt_units)
    return factor

def convert_altitude(alt, alt_units_in, alt_units_out):
    # type : (float, str, str) -> float
    """nominal unit is ft"""
    if alt_units_in == alt_units_out:
        return alt
    return alt * _altitude_factor(alt_units_in, alt_units_out)

def _altitude_factor(alt_units_in, alt_units_out):
    # type : (str, str) -> float
    """helper method for convert_altitude"""
    factor = 1.0
    # units to feet
    if alt_units_in == 'm':
        factor /= 0.3048
    elif alt_units_in == 'ft':
        pass
    elif alt_units_in == 'kft':
        factor *= 1000.
    else:
        raise RuntimeError('alt_units_in=%r is not valid; use [ft, m, kft]' % alt_units_in)

    # ft to m
    if alt_units_out == 'm':
        factor *= 0.3048
    elif alt_units_out == 'ft':
        pass
    elif alt_units_out == 'kft':
        factor /= 1000.
    else:
        raise RuntimeError('alt_units_out=%r is not valid; use [ft, m, kft]' % alt_units_out)
    return factor

def _reynolds_factor(reynolds_units_in, reynolds_units_out):
    # type : (str, str) -> float
    """helper method"""
    factor = 1.0
    # units to 1/feet
    if reynolds_units_in == '1/m':
        factor *= 0.3048
    elif reynolds_units_in == '1/ft':
        pass
    elif reynolds_units_in == '1/in':
        factor *= 12.
    else:
        msg = 'reynolds_units_in=%r is not valid; use [1/ft, 1/m, 1/in]' % reynolds_units_in
        raise RuntimeError(msg)

    # 1/ft to 1/m
    if reynolds_units_out == '1/m':
        factor /= 0.3048
    elif reynolds_units_out == '1/ft':
        pass
    elif reynolds_units_out == '1/in':
        factor /= 12.
    else:
        msg = 'reynolds_units_out=%r is not valid; use [1/ft, 1/m, 1/in]' % reynolds_units_out
        raise RuntimeError(msg)
    return factor

def convert_velocity(velocity, velocity_units_in, velocity_units_out):
    # type : (float, str, str) -> float
    """nominal unit is ft/s"""
    if velocity_units_in == velocity_units_out:
        return velocity
    return velocity * _velocity_factor(velocity_units_in, velocity_units_out)

def _velocity_factor(velocity_units_in, velocity_units_out):
    # type : (str, str) -> float
    """helper method for convert_velocity"""
    factor = 1.0
    if velocity_units_in == 'm/s':
        factor /= 0.3048
    elif velocity_units_in == 'ft/s':
        pass
    elif velocity_units_in == 'in/s':
        factor /= 12.
    elif velocity_units_in == 'knots':
        factor *= 1.68781
    else:
        msg = 'velocity_units_in=%r is not valid; use [ft/s, m/s, knots]' % velocity_units_in
        raise RuntimeError(msg)

    if velocity_units_out == 'm/s':
        factor *= 0.3048
    elif velocity_units_out == 'ft/s':
        pass
    elif velocity_units_out == 'in/s':
        factor *= 12.
    elif velocity_units_out == 'knots':
        factor /= 1.68781
    else:
        msg = 'velocity_units_out=%r is not valid; use [ft/s, m/s, in/s, knots]' % (
            velocity_units_out)
        raise RuntimeError(msg)
    return factor

def convert_pressure(pressure, pressure_units_in, pressure_units_out):
    # type : (float, str, str) -> float
    """nominal unit is psf"""
    if pressure_units_in == pressure_units_out:
        return pressure
    return pressure * _pressure_factor(pressure_units_in, pressure_units_out)

def _pressure_factor(pressure_units_in, pressure_units_out):
    # type : (str, str) -> float
    """helper method for convert_pressure"""
    factor = 1.0
    if pressure_units_in == 'psf':
        pass
    elif pressure_units_in == 'psi':
        factor *= 144
    elif pressure_units_in == 'Pa':
        factor /= 47.880172
    elif pressure_units_in == 'kPa':
        factor *= 20.88543815038
    elif pressure_units_in == 'MPa':
        factor *= 20885.43815038
    else:
        msg = 'pressure_units_in=%r is not valid; use [psf, psi, Pa, kPa, MPa]' % pressure_units_in
        raise RuntimeError(msg)

    if pressure_units_out == 'psf':
        pass
    elif pressure_units_out == 'psi':
        factor /= 144
    elif pressure_units_out == 'Pa':
        factor *= 47.880172
    elif pressure_units_out == 'kPa':
        factor /= 20.88543815038
    elif pressure_units_out == 'MPa':
        factor /= 20885.43815038
    else:
        raise RuntimeError('pressure_units_out=%r is not valid; use [psf, psi, Pa, kPa, MPa]' % (
            pressure_units_out))
    return factor

def convert_density(density, density_units_in, density_units_out):
    # type : (float, str, str) -> float
    """nominal unit is slug/ft^3"""
    if density_units_in == density_units_out:
        return density
    return density * _density_factor(density_units_in, density_units_out)

def _density_factor(density_units_in, density_units_out):
    # type : (str, str) -> float
    """helper method for convert_density"""
    factor = 1.0
    if density_units_in == 'slug/ft^3':
        pass
    elif density_units_in == 'slinch/in^3':
        factor *= 12**4
    elif density_units_in == 'kg/m^3':
        factor /= 515.378818
    else:
        msg = 'density_units_in=%r is not valid; use [slug/ft^3]' % density_units_in
        raise RuntimeError(msg)

    # data is now in slug/ft^3
    if density_units_out == 'slug/ft^3':
        pass
    elif density_units_out == 'slinch/in^3':
        factor /= 12**4
    elif density_units_out == 'kg/m^3':
        factor *= 515.378818
    else:
        msg = 'density_units_out=%r is not valid; use [slug/ft^3, slinch/in^3]' % density_units_out
        raise RuntimeError(msg)
    return factor

def _temperature_factor(temperature_units_in, temperature_units_out):
    if temperature_units_in == temperature_units_out:
        return 1.
    elif  temperature_units_in == 'R' and temperature_units_out == 'K':
        return 5. / 9.
    elif  temperature_units_in == 'R' and temperature_units_out == 'K':
        return 9. / 5.
    else:
        raise NotImplementedError('temperature_units_in=%r temperature_units_out=%r '
                                  'is not supported' % (
                                      temperature_units_in, temperature_units_out))
