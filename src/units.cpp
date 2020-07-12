#include <string>

#include "units.hpp"
#include "io.hpp"

// critical density of the Universe; [kg / m^3]
constexpr double CRITICAL_DENSITY = 8.5E-27;
// constexpr double CRITICAL_DENSITY = 1;

// Hubble constants; [1/s]
constexpr double H0 = 4.55E17;

// radius of the sun; [m]
constexpr double SUN_RADIUS = 696340000;

// mass of the sun; [kg]
constexpr double SUN_MASS = 1.989E30;

// gravitational constatnt; natural units [1/eV^2]
constexpr double GRAV_CON = 6.71E-57;

// speed of the light; [m/s]
constexpr double LIGHT_SPEED   = 299792458L; // m / s

// electron charge; [C]
constexpr double ELECTRON_C = 1.60217663E-19;

// reduced Planck constant; [ev s]
constexpr double RED_PLANCK_CON = 6.582E-16;

// Planck mass; natural units [eV]
constexpr double M_PL = sqrt(1 / (8 * M_PI*GRAV_CON));
// constexpr double M_PL = 2.435E27;
// constexpr double M_PL = 1;

// conversion from kg*c^2 to eV; [eV / kg]
constexpr double KG_TO_EV = LIGHT_SPEED * LIGHT_SPEED / ELECTRON_C;

// conversion from m to 1/eV; [eV / m]
constexpr double M_TO_EV = RED_PLANCK_CON / (LIGHT_SPEED * M_PL);

// kiloparsec; [m]
constexpr double KPC = 3.086E19; // m

// halo mass multiplier from input value
constexpr double HALO_MASS_MLT = 1E12;

// computational units for mass in case of NFW halo; [eV]
// constexpr double MASS_UNITS_HALO = GRAV_CON * KG_TO_EV / KPC;
// constexpr double MASS_UNITS_HALO = GRAV_CON * KG_TO_EV;
constexpr double MASS_UNITS_HALO = KG_TO_EV;

// computational units for mass in case of star; [eV]
// constexpr double MASS_UNITS_STAR = GRAV_CON * KG_TO_EV / SUN_RADIUS;
// constexpr double MASS_UNITS_STAR = GRAV_CON * KG_TO_EV;
constexpr double MASS_UNITS_STAR = KG_TO_EV;

// computational units for density in case of NFW halo; [eV/kpc^3]
// constexpr double DENSITY_UNITS_HALO = MASS_UNITS_HALO *  (KPC * KPC * KPC);
constexpr double DENSITY_UNITS_HALO = MASS_UNITS_HALO;

// computational units for density in case of star; [eV/R_sun^3]
// constexpr double DENSITY_UNITS_STAR = MASS_UNITS_STAR * (SUN_RADIUS * SUN_RADIUS * SUN_RADIUS);
constexpr double DENSITY_UNITS_STAR = MASS_UNITS_STAR;

double phi_prefactor;
double chi_prefactor;

static double get_chi_0()
{
    return 2 * param.chi_opt.beta*M_PL*param.chi_opt.Ys;
}

/* convertes star radius in units of the sun
   to computing units [R_sun] */
double star_radius_to_cu(double star_rad)
{
    return star_rad * SUN_RADIUS;
}

/* convertes halo mass in units of the sun * 1E12
   to computing units [kpc] */
double halo_mass_to_cu(double mass_halo)
{
    double mass_kg = mass_halo * HALO_MASS_MLT * SUN_MASS; // kg
    return mass_kg * MASS_UNITS_HALO; // eV
}


/* convertes star mass in units of the sun
   to computing units [R_sun] */
double star_mass_to_cu(double mass_halo)
{
    double mass_kg = mass_halo * SUN_MASS; // kg
    return mass_kg * MASS_UNITS_STAR; // eV
}

/* converts density in units of the critical density
   to computing units [kpc2]
*/
double density_to_halo_cu(double omega_m)
{
    double density = omega_m * CRITICAL_DENSITY; // kg / m3
    return density * DENSITY_UNITS_HALO;
}

/* converts density in units of the critical density
   to computing units [R_sun2]
*/
double density_to_sun_cu(double omega_m)
{
    double density = omega_m * CRITICAL_DENSITY; // kg / m3
    return density * DENSITY_UNITS_STAR;
}

/* converts chameleon field in units of chi_0,
   i.e. computing units, to physical units
*/
double chi_cu_to_phys(double chi_cu)
{
    return chi_cu*get_chi_0();
}


/* get prefactors for poisson / chameleon EoM */
void get_phi_prefactor()
{
    double R_units;
    switch (param.integration.mod)
    {
        case MOD_STAR:
        {
           R_units = SUN_RADIUS;
           break;
        }
        case MOD_NFW:
        {
           R_units = KPC;
           break;
        }
    }
    
    R_units = 1;
    phi_prefactor = 3/2.*param.spatial.Omega_m*pow(H0 * R_units / LIGHT_SPEED ,2) * param.spatial.rho_0;
    phi_prefactor /= 4*M_PI;
    phi_prefactor = GRAV_CON;
}

/* get prefactors for poisson / chameleon EoM */
void get_chi_prefactor()
{
    // first initialize phi_prefactor
    if (!phi_prefactor) get_phi_prefactor();

    double R_units;
    switch (param.integration.mod)
    {
        case MOD_STAR:
        {
           R_units = SUN_RADIUS;
           break;
        }
        case MOD_NFW:
        {
           R_units = KPC;
           break;
        }
    }
    
    R_units = 1;
    chi_prefactor = 4*M_PI*phi_prefactor*pow(R_units ,2) / param.chi_opt.Ys; // / param.spatial.rho_0;
}

double get_force_mlt()
{
    return -get_chi_0()*param.chi_opt.beta/M_PL;
}

double get_pot_mlt()
{
    return param.chi_opt.Ys;
}


double get_mass_scale()
{
    switch (param.integration.mod)
    {
        case MOD_STAR:
        {
            return SUN_MASS * MASS_UNITS_STAR;
        }
        case MOD_NFW:
        {
            return SUN_MASS * MASS_UNITS_STAR * HALO_MASS_MLT;
        }
    }
}

double get_radius_scale()
{
    // return 1;
    switch (param.integration.mod)
    {
        case MOD_STAR:
        {
            return SUN_RADIUS;
        }
        case MOD_NFW:
        {
            return KPC;
        }
    }

}