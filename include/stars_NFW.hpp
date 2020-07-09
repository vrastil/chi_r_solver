#pragma once

// SPATIAL FUNCTIONS
double get_R200();
double pot_star(double r);
double force_star(double r);
double rho_NFW(double r);
double rho_star(double r);
double rho_r(double r);
double pot_NFW(double r);
double force_NFW(double r);
double pot_new(double r);
double force_new(double r);
double r_eq_star();
double r_eq_NFW();
double get_r_eq();
double get_rho_c(double _r_eq);
double get_Ys();

// CHAMELEON FUNCTIONS
double chi_bulk(double rho);
double chi_bulk_der(double r);
double chi_bulk_laplace(double r);
double chi_bulk_r(double r);
double V_eff_der(double r, double chi);
double V_eff_2nd_derr(double chi);
double chi_0_der(double r, double *chi);
double chi_0_der_lin(double r, double *chi, double a);
double chi_1_der(double r, double *chi);
double chi_0_der_NFW_lin(double r, double *chi);
double chi_1_der_NFW_lin(double r, double *chi);


// SHOOTING METHODS
bool is_integrate(double r, double *y, double r_max);
void fce_min_r2(double s, double &t_0, double *y_0, double eps, double m);
void fce_min_chi_0(double s, double &t_0, double *y_0);
double fce_max(double r, double *y);
double shoot_meth_star(double err, double eps);
double shoot_meth_NFW(double err, double eps);
void slv_Chameleon(double r_max, double chi_pot_0, double chi_der_0, double err, double **chi, double step, int &N_i, int reg);
double get_chi_0_star(double err, double eps);
double get_chi_0_NFW(double err, double eps);
void slv_Chameleon_star(double r_max, double err, double **chi, double &step, int &N_i, int reg);
void slv_Chameleon_star(double r_max, double err, double **chi, double &step, int &N_i, int reg, int *cout_max);
void slv_Chameleon_NFW(double r_max, double err, double **chi, double &step, int &N_i, int reg);
void slv_Chameleon_NFW(double r_max, double err, double **chi, double &step, int &N_i, int reg, int *cout_max);