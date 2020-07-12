#include <string>

#include "units.hpp"
#include "io.hpp"

// constants
constexpr double CRITICAL_DENSITY = 8.5E-27; // kg / m3
constexpr double SUN_RADIUS     = 696340000; // m
constexpr double SUN_MASS        = 1.989E30; // kg
constexpr double GRAV_CON     = 6.67384E-11; // m3 / kg s2
constexpr double LIGHT_SPEED   = 299792458L; // m / s
constexpr double RED_PLANCK =     1.054E-34; // Js = kg m2 /s2
constexpr double M_PL            = 4.341E-9; // kg
constexpr double KPC             = 3.086E19; // m
constexpr double HALO_MASS_MLT       = 1E12; // halo mass multiplier

// derived
constexpr double GRAV_CON_C1 = GRAV_CON / (LIGHT_SPEED * LIGHT_SPEED); // m / kg; c = 1
constexpr double MASS_UNITS_HALO = GRAV_CON * KPC; // kpc / kg
constexpr double MASS_UNITS_STAR = GRAV_CON * SUN_RADIUS; // R_sun / kg
constexpr double DENSITY_UNITS_HALO = MASS_UNITS_HALO *  (KPC * KPC * KPC);
constexpr double DENSITY_UNITS_STAR = MASS_UNITS_STAR * (SUN_RADIUS * SUN_RADIUS * SUN_RADIUS);

// unit strings
// constexpr std::string STAR_MASS_UNIT = "M_s";
// constexpr std::string HALO_MASS_UNIT = "E12 M_s";
// constexpr std::string HALO_DENSITY_UNIT = "E12 M_s";

static double get_chi_0()
{
    return 2 * param.chi_opt.beta*M_PL*param.chi_opt.Ys;
}

/* convertes halo mass in units of the sun * 1E12
   to computing units [kpc] */
double halo_mass_to_cu(double mass_halo)
{
    double mass_kg = mass_halo * HALO_MASS_MLT * SUN_MASS; // kg
    return mass_kg * MASS_UNITS_HALO;
}


/* convertes star mass in units of the sun
   to computing units [R_sun] */
double star_mass_to_cu(double mass_halo)
{
    double mass_kg = mass_halo * SUN_MASS; // kg
    return mass_kg * MASS_UNITS_STAR;
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


/* get prefactors for chameleon EoM */
double get_v_eff_der_prefactor()
{
    double prefactor = param.chi_opt.beta / (pow(M_PL * get_chi_0(), 2));
    // this will cancel with units of density
    switch (param.integration.mod)
    {
        case MOD_STAR:
        {
            prefactor *= DENSITY_UNITS_STAR;
            break;
        }
        case MOD_NFW:
        {
            prefactor *= DENSITY_UNITS_HALO;
            break;
        }
    }
    return prefactor;
}

double get_v_eff_2nd_der_prefactor()
{
    // not used KG_M3_TO_R_SUN2 / KG_M3_TO_KPC2 conversion
    // return mass in the right computing units such that [m*R] = 1
    return (1 - param.chi_opt.n) * param.chi_opt.beta / M_PL * param.spatial.rho_0 / get_chi_0();
}

double get_star_2a_prefactor()
{
    switch (param.integration.mod)
    {
        case MOD_STAR:
        {
            // density in 1/R_sun^2
            return 2 * param.chi_opt.Ys * M_PL * M_PL / DENSITY_UNITS_STAR; 
        }
        case MOD_NFW:
        {
            // density in 1/pc^2
            return 2 * param.chi_opt.Ys * M_PL * M_PL / DENSITY_UNITS_HALO;
        }
    }
    
}