#include <math.h>

# include "constants.hpp"

const double PI = 4 * atan(1);
// const double G = 6.67384E-11; // gravitation constatnt
 const double G = 6.71E-57; // gravitation constatnt, eV
const double Mp = sqrt(1 / (8 * PI*G));				// Planck mass
const double B = sqrt(1 / 6.0); // chameleon coupling constant

//////// RADIAL OBJECTS ///////////
int mod; // 0 = star, 1 = NFW
double n;
double chi_0;
double Ys;
double R; // scale radius
double R200; // virial radius
double c = 100; // concentration parameter
double rho_c = 1;
double Ms;
double r_eq;
double rho_0;
double chi_B;
double m_inf;

int h_N;
int h_re;

const double M_L = 2.4e-3;
const double n_pl = 8.0;
const double gamma_chi = 1;
const double rho_pl = 1e-11;
const double phi_min = M_L*sgn(1 - n_pl)*pow(B*rho_pl/(abs(n_pl)*gamma_chi*pow(M_L,3)*Mp),1/(1-n_pl));
const double phi_s = phi_min*(1 - 1 / n_pl);
const double K = sqrt((n_pl - 2)*(n_pl - 2)*gamma_chi*pow(M_L, 4 - n_pl)*pow(abs(phi_s), n_pl - 2) / 2);
const double m_pl = 3;



int sgn(double x){
	if (x == 0) return 0;
	else if (x > 0) return 1;
	else return -1;
}