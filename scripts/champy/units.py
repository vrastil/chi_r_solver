import numpy as np

CRITICAL_DENSITY = 8.5E-27 # kg/m3
SUN_MASS = 1.989E30 # kg
KPC = 3.086E19 # m

HALO_MASS_MLT = 1E12
DISTANCE_UNITS_HALO = KPC

M_PI = np.pi

def get_omega_for_nfw(R_vir, M_vir):
    """ From virial radius (in kpc) of the halo and its virial mass (in 1E12 M_Sun)
        get average density, in units of critical density """
    # get values in SI units, kg, m
    M_vir *= HALO_MASS_MLT * SUN_MASS
    R_vir *= DISTANCE_UNITS_HALO

    # get simple average density
    rho = M_vir / (4.0/3.0*M_PI*pow(R_vir, 3))

    # calculate omeg
    omega = rho / CRITICAL_DENSITY
    return omega