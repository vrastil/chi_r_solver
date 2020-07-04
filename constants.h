//////// CONSTANTS //////
////////////
#pragma once
extern const double PI;
// extern const double G = 6.67384E-11; // gravitation constatnt
extern const double G; // gravitation constatnt, eV
extern const double Mp;				// Planck mass
extern const double B; // chameleon coupling constant

//////// RADIAL OBJECTS ///////////
extern int mod; // 0 = star, 1 = NFW
extern double n;
extern double chi_0;
extern double Ys;
extern double R; // scale radius
extern double R200; // virial radius
extern double c; // concentration parameter
extern double rho_c;
extern double Ms;
extern double r_eq;
extern double rho_0;
extern double chi_B;
extern double m_inf;

extern int h_N;
extern int h_re;

extern const double K;

int sgn(double x);