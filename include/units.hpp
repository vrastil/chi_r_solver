#pragma once

/* convertes star radius in units of the sun
   to computing units [m] */
double star_radius_to_cu(double star_rad);

/* convertes halo radius in units of kpc
   to computing units [m] */
double halo_radius_to_cu(double halo_rad);

/* convertes halo mass in units of the sun * 1E12
   to computing units */
double halo_mass_to_cu(double mass_halo);

/* convertes star mass in units of the sun
   to computing units */
double star_mass_to_cu(double mass_halo);

/* converts density in units of the critical density
   to computing units
*/
double density_to_halo_cu(double omega_m);

/* converts density in units of the critical density
   to computing units
*/
double density_to_sun_cu(double omega_m);


/*
*/
double get_force_mlt();
double get_pot_mlt();

double get_mass_scale();
double get_radius_scale();

double get_star_2ra_prefactor();

// prefactor for gravitational potential
// G*rho_m,0 + computing units
void get_phi_prefactor();
extern double phi_prefactor;

// dimensionless prefactor to poisson equation for chameleon field
// beta*rho_m,0 / (Mpl*chi_0) + computing units
void get_chi_prefactor();
extern double chi_prefactor;