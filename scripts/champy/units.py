import numpy as np

CRITICAL_DENSITY = 8.5E-27 # kg/m3
SUN_MASS = 1.989E30 # kg
KPC = 3.086E19 # m
LIGHT_SPEED   = 299792458 # m / s
ELECTRON_C = 1.60217663E-19 # C
KG_TO_EV = LIGHT_SPEED * LIGHT_SPEED / ELECTRON_C

HALO_MASS_MLT = 1E12
DISTANCE_UNITS_HALO = KPC
DENSITY_UNITS_HALO = pow(KPC, 3) / SUN_MASS # KG_TO_EV

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

def density_to_halo_cu(Omega_m):
    density = Omega_m * CRITICAL_DENSITY
    return density * DENSITY_UNITS_HALO


def halo_mass_to_cu(M200_sun):
    mass_kg = M200_sun * HALO_MASS_MLT * SUN_MASS
    return mass_kg * DENSITY_UNITS_HALO

def get_rho_c(Omega_m, c):
    rho_0 = density_to_halo_cu(Omega_m)
    return 200 * rho_0*pow(c, 3) / (3*(np.log(1+c)-c/(1+c)))

def get_r_scale(rho_c, c, M200_sun):
    M200 = halo_mass_to_cu(M200_sun)
    return pow(M200/(4*M_PI*rho_c)/(np.log(1+c)-c/(c+1)) , 1 / 3.0)

def get_halo_parameters(params, Omega_m=1):
    params['rho_0'] = density_to_halo_cu(Omega_m)
    params['rho_c'] = get_rho_c(Omega_m, params['c'])
    params['R_s'] = get_r_scale(params['rho_c'], params['c'], params['M200_sun']) / DISTANCE_UNITS_HALO

def get_halo_mass_function(params, Omega_m=1):
    if 'R_s' not in params or 'rho_c' not in params:
        get_halo_parameters(params, Omega_m=Omega_m)

    rs = params['R_s']
    
    return lambda r: 4*M_PI*params['rho_c']*pow(rs, 3)*(np.log(1+r/rs) - (r/(r+rs))) / HALO_MASS_MLT
