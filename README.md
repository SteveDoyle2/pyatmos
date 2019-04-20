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
