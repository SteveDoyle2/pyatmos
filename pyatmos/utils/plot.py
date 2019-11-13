import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')

import matplotlib.pyplot as plt
from pyatmos.utils.atmosphere_vectorized import (
    atm_pressure_array, atm_temperature_array, atm_dynamic_viscosity_mu_array,
    atm_kinematic_viscosity_nu_array, atm_equivalent_airspeed_array,
    sutherland_viscoscity_array,
)

alt = np.linspace(-100000, 300000., num=1000)
plt.figure(1)
plt.plot(atm_pressure_array(alt, alt_units='ft', pressure_units='psf'), alt)
plt.grid(True)
plt.xlabel('pressure (psi)')
plt.ylabel('alt (ft)')

plt.figure(1)
plt.semilogx(atm_pressure_array(alt, alt_units='ft', pressure_units='psf'), alt)
plt.xlabel('log(pressure) (psi)')
plt.ylabel('alt (ft)')
plt.grid(True)


plt.figure(2)
plt.plot(atm_temperature_array(alt, alt_units='ft', temperature_units='R'), alt)
plt.xlabel('temperature (R)')
plt.ylabel('alt (ft)')
plt.grid(True)


plt.figure(3)
plt.plot(atm_dynamic_viscosity_mu_array(alt, alt_units='ft', visc_units='(lbf*s)/ft^2'), alt)
plt.xlabel('\nu (lbf*s/ft^2)')
plt.ylabel('alt (ft)')
plt.grid(True)

plt.figure(4)
plt.semilogx(atm_kinematic_viscosity_nu_array(alt, alt_units='ft', visc_units='ft^2/s'), alt)
plt.xlabel('\nu (lbf*s/ft^2)')
plt.ylabel('alt (ft)')
plt.grid(True)


plt.figure(5)
machs = [1e-5, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 2., 3., 4., 5.]
for mach in machs:
    plt.plot(atm_equivalent_airspeed_array(alt, mach, alt_units='ft', eas_units='ft/s'), alt, label='M=%s' % mach)
plt.xlabel('eas (ft/s)')
plt.ylabel('alt (ft)')
plt.grid(True)
plt.legend()

plt.figure(6)
for mach in machs:
    plt.semilogx(atm_equivalent_airspeed_array(alt, mach, alt_units='ft', eas_units='ft/s'), alt, label='M=%s' % mach)
plt.xlabel('eas (ft/s)')
plt.ylabel('alt (ft)')
plt.grid(True)
plt.legend()

plt.figure(7)
alts2 = np.linspace(-50000., 100000, num=11)
mach = np.linspace(0., 2., num=21)
for alti in alts2:
    plt.plot(atm_equivalent_airspeed_array(alti, mach, alt_units='ft', eas_units='ft/s'), mach, label='alt=%s' % alti)
plt.xlabel('eas (ft/s)')
plt.ylabel('mach')
plt.grid(True)
plt.legend()

plt.figure(8)
alts2 = np.linspace(-50000., 100000, num=11)
mach = np.linspace(0., 2., num=21)
for alti in alts2:
    plt.semilogx(atm_equivalent_airspeed_array(alti, mach, alt_units='ft', eas_units='ft/s'), mach, label='alt=%s' % alti)
plt.xlabel('eas (ft/s)')
plt.ylabel('mach')
plt.grid(True)
plt.legend()


T = np.linspace(0., 1000., num=1000)

plt.figure(9)
plt.plot(sutherland_viscoscity_array(T), alt)
plt.xlabel('visc (ft/s)')
plt.ylabel('alt (ft)')
plt.grid(True)

plt.figure(10)
plt.semilogx(sutherland_viscoscity_array(T), alt)
plt.xlabel('visc (ft/s)')
plt.ylabel('alt (ft)')
plt.grid(True)

plt.show()
