#pragma once

/* convertes halo mass in units of the sun * 1E12
   to computing units [kpc] */
double halo_mass_to_cu(double mass_halo);

/* convertes star mass in units of the sun
   to computing units [R_sun] */
double star_mass_to_cu(double mass_halo);

/* converts density in units of the critical density
   to computing units [kpc]
*/
double density_to_halo_cu(double omega_m);

/* converts density in units of the critical density
   to computing units [R_sub]
*/
double density_to_sun_cu(double omega_m);

/* converts chameleon field in units planck mass,
   i.e. computing units, to physical units [R_sub]
*/
double chi_cu_to_phys(double chi_cu);

/*
*/
double get_v_eff_der_prefactor();
double get_v_eff_2nd_der_prefactor();
double get_star_2a_prefactor();
