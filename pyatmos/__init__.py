"""
defines a colorama log
"""
# coding: utf-8
from __future__ import print_function, unicode_literals, absolute_import
#import sys
#import os

__version__ = '1.0'
__desc__ = 'pyatmos'
__long__ = __desc__
__website__ = 'https://github.com/pyatmos/pyatmos'
__license__ = 'BSD-3'
__author__ = ''
__email__ = ''

from .atmosphere import (
    atm_density, atm_dynamic_pressure, atm_temperature,
    atm_pressure, atm_velocity, atm_mach, atm_equivalent_airspeed,
    atm_dynamic_viscosity_mu, atm_kinematic_viscosity_nu,
    get_alt_for_density, get_alt_for_pressure,
    get_alt_for_q_with_constant_mach,
    get_alt_for_eas_with_constant_mach,
    atm_unit_reynolds_number, atm_unit_reynolds_number2,
    make_flfacts_alt_sweep, make_flfacts_mach_sweep,
    make_flfacts_eas_sweep,
)
