# pyatmos
High altitude atmosphere library

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

https://pyatmos.readthedocs.io/en/latest/

|  Version  | Docs  | Status |
| :--- 	  | :--- 	  | :--- 	  |
|   Master | [![Documentation Status](https://readthedocs.org/projects/pynastran-git/badge/?version=latest)](https://pyatmos.readthedocs.io/en/latest/#) | [![Linux Status](https://img.shields.io/travis/SteveDoyle2/pyatmos/master.svg)](https://travis-ci.org/SteveDoyle2/pyatmos) ![Coverage Status](https://coveralls.io/repos/github/SteveDoyle2/pyatmos/badge.svg?branch=master) |

<!---

|  [![PyPi Version](https://img.shields.io/pypi/v/pyatmos.svg)](https://pypi.python.org/pypi/pyatmos) | [docs](http://pynastran.m4-engineering.com/1.1.0/) | [![Build Status](https://img.shields.io/travis/SteveDoyle2/pyatmos/v1.0.svg)](https://travis-ci.org/SteveDoyle2/pyatmos) [![Coverage Status](https://img.shields.io/coveralls/SteveDoyle2/pyatmos/v1.0.svg)](https://coveralls.io/github/SteveDoyle2/pyatmos?branch=v1.0) |

--->
