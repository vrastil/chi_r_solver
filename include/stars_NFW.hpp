#pragma once

// SPATIAL FUNCTIONS
double get_R200();
double get_r_eq();
double get_rho_c(double _r_eq);
double get_Ys();

// CHAMELEON FUNCTIONS
double chi_bulk(double rho);
double V_eff_2nd_derr(double chi);

// SOLVER
