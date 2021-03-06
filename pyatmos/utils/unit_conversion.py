"""defines various internal unit conversio functions"""


def _feet_to_alt_units(alt_units: str) -> float:
    """helper method"""
    if alt_units == 'm':
        factor = 0.3048
    elif alt_units == 'ft':
        factor = 1.
    else:
        raise RuntimeError(f'alt_units={alt_units} is not valid; use [ft, m]')
    return factor


def _rankine_to_temperature_units(temperature_units: str) -> float:
    if temperature_units == 'R':
        factor = 1.
    elif temperature_units == 'K':
        factor = 5. / 9.
    else:
        raise RuntimeError('temperature_units=%r is not valid; use [R, K]' % temperature_units)
    return factor


def _psfs_to_dvisc_units(visc_units: str) -> float:
    """same units as pressure, except multiplied by seconds"""
    if visc_units == '(lbf*s)/ft^2':
        factor = 1.
    elif visc_units in ['(N*s)/m^2', 'Pa*s']:
        factor = 47.88026
    else:
        raise RuntimeError('visc_units=%r; not in (lbf*s)/ft^2, (N*s)/m^2, or Pa*s')
    return factor

def _ft2s_to_kvisc_units(alt_units: str, visc_units: str) -> float:
    if visc_units == 'ft^2/s':
        factor = 1.
    elif visc_units == 'm^2/s':
        factor = _feet_to_alt_units(alt_units) ** 2  # TODO: double check this...
    else:
        raise NotImplementedError('visc_units=%r' % visc_units)
    return factor

def convert_altitude(alt: float, alt_units_in: str, alt_units_out: str) -> float:
    """
    Nominal unit is ft.

    Parameters
    ----------
    alt : float
        The altitude in alt_units_in.
    alt_units_in : str
        Original alt unit.
    alt_units_out : str
        Nominal alt unit.

    Returns
    -------
    alt : float
        Altitude in alt_units_out.

    """
    if alt_units_in == alt_units_out:
        return alt
    return alt * _altitude_factor(alt_units_in, alt_units_out)

def _altitude_factor(alt_units_in: str, alt_units_out: str) -> float:
    """
    Helper method for convert_altitude.

    Parameters
    ----------
    alt_units_in : str
        Original alt unit.
    alt_units_out : str
        Nominal alt unit.

    Returns
    -------
    factor : float
        Conversion factor to nominal altitude unit.

    """
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

def _reynolds_factor(reynolds_units_in: str, reynolds_units_out: str) -> float:
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
        raise RuntimeError('reynolds_units_in=%r is not valid; use '
                           '[1/ft, 1/m, 1/in]' % reynolds_units_in)

    # 1/ft to 1/m
    if reynolds_units_out == '1/m':
        factor /= 0.3048
    elif reynolds_units_out == '1/ft':
        pass
    elif reynolds_units_out == '1/in':
        factor /= 12.
    else:
        raise RuntimeError('reynolds_units_out=%r is not valid; use '
                           '[1/ft, 1/m, 1/in]' % reynolds_units_out)
    return factor

def convert_velocity(velocity: float, velocity_units_in: str, velocity_units_out: str) -> float:
    """nominal unit is ft/s"""
    if velocity_units_in == velocity_units_out:
        return velocity
    return velocity * _velocity_factor(velocity_units_in, velocity_units_out)

def _velocity_factor(velocity_units_in: str, velocity_units_out: str) -> float:
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
        raise RuntimeError(f'velocity_units_in={velocity_units_in} is not valid; use '
                           '[ft/s, m/s, in/s, knots]')

    if velocity_units_out == 'm/s':
        factor *= 0.3048
    elif velocity_units_out == 'ft/s':
        pass
    elif velocity_units_out == 'in/s':
        factor *= 12.
    elif velocity_units_out == 'knots':
        factor /= 1.68781
    else:
        raise RuntimeError(f'velocity_units_out={velocity_units_out!r} is not valid; use '
                           '[ft/s, m/s, in/s, knots]')
    return factor

def convert_pressure(pressure: float, pressure_units_in: str, pressure_units_out: str) -> float:
    """nominal unit is psf"""
    if pressure_units_in == pressure_units_out:
        return pressure
    return pressure * _pressure_factor(pressure_units_in, pressure_units_out)

def _pressure_factor(pressure_units_in: str, pressure_units_out: str) -> float:
    """helper method for convert_pressure"""
    factor = 1.0
    if pressure_units_in in ['psf', 'lb/ft^2']:
        pass
    elif pressure_units_in in ['psi', 'lb/in^2']:
        factor *= 144
    elif pressure_units_in == 'Pa':
        factor /= 47.880172
    elif pressure_units_in == 'kPa':
        factor *= 20.88543815038
    elif pressure_units_in == 'MPa':
        factor *= 20885.43815038
    else:
        raise RuntimeError(f'pressure_units_in={pressure_units_in!r} is not valid; use '
                           '[psf, psi, Pa, kPa, MPa]')

    if pressure_units_out in ['psf', 'lb/ft^2']:
        pass
    elif pressure_units_out in ['psi', 'lb/in^2']:
        factor /= 144
    elif pressure_units_out == 'Pa':
        factor *= 47.880172
    elif pressure_units_out == 'kPa':
        factor /= 20.88543815038
    elif pressure_units_out == 'MPa':
        factor /= 20885.43815038
    else:
        raise RuntimeError(f'pressure_units_out={pressure_units_out!r} is not valid; use '
                           '[psf, psi, Pa, kPa, MPa]')
    return factor

def convert_density(density: float, density_units_in: str, density_units_out: str) -> float:
    """nominal unit is slug/ft^3"""
    if density_units_in == density_units_out:
        return density
    return density * _density_factor(density_units_in, density_units_out)

def _density_factor(density_units_in: str, density_units_out: str) -> float:
    """helper method for convert_density"""
    factor = 1.0
    if density_units_in == 'slug/ft^3':
        pass
    elif density_units_in == 'slinch/in^3':
        factor *= 12**4
    elif density_units_in == 'kg/m^3':
        factor /= 515.378818
    else:
        raise RuntimeError(f'density_units_in={density_units_in!r} is not valid; use '
                           '[slug/ft^3, slinch/in^3, kg/m^3]')

    # data is now in slug/ft^3
    if density_units_out == 'slug/ft^3':
        pass
    elif density_units_out == 'slinch/in^3':
        factor /= 12**4
    elif density_units_out == 'kg/m^3':
        factor *= 515.378818
    else:
        raise RuntimeError(f'density_units_out={density_units_out} is not valid; use '
                           '[slug/ft^3, slinch/in^3, kg/m^3]')
    return factor

def _temperature_factor(temperature_units_in: str, temperature_units_out: str) -> float:
    if temperature_units_in == temperature_units_out:
        factor = 1.
    elif  temperature_units_in == 'R' and temperature_units_out == 'K':
        factor = 5. / 9.
    elif  temperature_units_in == 'K' and temperature_units_out == 'R':
        factor = 9. / 5.
    else:
        raise NotImplementedError('temperature_units_in=%r temperature_units_out=%r '
                                  'is not supported' % (
                                      temperature_units_in, temperature_units_out))
    return factor
