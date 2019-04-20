"""
Contains functions that have both consistent units and are vectorized:

 - q = dynamic_pressure_p_mach(p, mach)
 - q = dynamic_pressure_rho_vel(rho, vel)
 - sos = speed_of_sound(T, R=1716., gamma=1.4)

"""
from __future__ import print_function, absolute_import


def dynamic_pressure_p_mach(p, mach):
    """Calculates dynamic pressure without options for units"""
    q = 0.7 * p * mach ** 2
    return q

def dynamic_pressure_rho_vel(rho, vel):
    """Calculates dynamic pressure without options for units"""
    q = 0.5 * rho * vel ** 2
    return q

def speed_of_sound(T, R=1716., gamma=1.4):
    """
    Calculates the speed of sound without options for units

    Parameters
    ----------
    T : float, np.ndarray
        the temperature
    R : float; default=1716.0
        1716.59, dry air, R=287.04 J/kg*K
    gamma : float; default=1.4
        the ratio of Cp/Cv

    """
    a = (gamma * R * T) ** 0.5
    return a
