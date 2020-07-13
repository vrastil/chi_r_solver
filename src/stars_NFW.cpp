#include  <cstdlib>
#include <functional>
#include <iostream>
#include <fstream>
#include <cmath>
#include "limits.h"

#define BOOST_LOG_DYN_LINK 1
#include <boost/log/trivial.hpp>

#include "stars_NFW.hpp"
#include "integrator.hpp"
#include "io.hpp"
#include "units.hpp"

using namespace std;

#define mod param.integration.mod
// #define chi (param.integration.chi)

#define R param.spatial.R
#define R200 param.spatial.R200
#define Omega_m param.spatial.Omega_m
#define rho_0 param.spatial.rho_0
#define rho_c param.spatial.rho_c
#define Ms param.spatial.Ms
#define M200 param.spatial.M200
#define M200_sun param.spatial.M200_sun
#define c param.spatial.c

#define Ys param.chi_opt.Ys
#define chi_B param.chi_opt.chi_B
#define n_chi param.chi_opt.n
#define beta param.chi_opt.beta
#define m_inf param.chi_opt.m_inf
#define R_eq param.spatial.R_eq

// SPATIAL FUNCTIONS
static double pot_star(double r);
static double force_star(double r);
static double rho_NFW(double r);
static double rho_star(double r);
static double rho_r(double r);
static double pot_NFW(double r);
static double force_NFW(double r);
static double pot_new(double r);
static double force_new(double r);
static double r_eq_star();
static double r_eq_NFW();
static double get_R200();

// CHAMELEON FUNCTIONS
static double chi_bulk(double rho);
static double chi_bulk_der(double r);
static double chi_bulk_laplace(double r);
static double chi_bulk_r(double r);
static double V_eff_der(double r, double chi);
static double V_eff_2nd_derr(double chi);
static double chi_mass2(double chi);
static double chi_mass(double chi);
static double chi_0_der(double r, double *chi);
static double chi_0_der_lin(double r, double *chi, double a);
static double chi_1_der(double r, double *chi);
static double chi_0_der_NFW_lin(double r, double *chi);
static double chi_1_der_NFW_lin(double r, double *chi);
static double get_r_eq();
static double get_rho_c(double _r_eq);
static double get_Ys();

// SHOOTING METHODS
static bool is_integrate(double r, double *y, double r_max);
static void fce_min_r2(double s, double &t_0, double *y_0, double eps, double m);
static void fce_min_chi_0(double s, double &t_0, double *y_0);
static double fce_max(double r, double *y);

static double shoot_meth_star(double err, double eps);
static double shoot_meth_NFW(double err, double eps);

static void slv_Chameleon(double r_max, double chi_pot_0, double chi_der_0, double err, chi_t& chi, double step);

static double get_chi_0_star(double err, double eps);
static double get_chi_0_NFW(double err, double eps);

static void slv_Chameleon_star(double r_max, double err, chi_t& chi, double &step);
static void slv_Chameleon_star_cout(double r_max, double err, chi_t& chi, double &step);
static void slv_Chameleon_NFW(double r_max, double err, chi_t& chi, double &step);
static void slv_Chameleon_NFW_cout(double r_max, double err, chi_t& chi, double &step);

void solve()
{
	// extract / define variables
	const double r_max = param.integration.r_max;
	const double err = param.integration.err;
	double step = param.integration.step;
	chi_t& chi = param.integration.chi;

	switch (mod){
		case MOD_STAR:
		{
			BOOST_LOG_TRIVIAL(info) << "Start solving the chameleon field for STAR.";
			slv_Chameleon_star_cout(r_max, err, chi, step);
			break;
		}
		case MOD_NFW:
		{
			BOOST_LOG_TRIVIAL(info) << "Start solving the chameleon field for NFW.";
			slv_Chameleon_NFW_cout(r_max, err, chi, step);
			break;
		}
	}
}

void spatial_init()
{
    switch (mod)
    {
        // STAR: from radius and mass of the star compute its density
        case MOD_STAR:
        {
			// convert radius input to computing units
			R = star_radius_to_cu(R);

            // convert mass input to computing units, for star Ms = M200
            Ms = M200 = star_mass_to_cu(M200_sun);

            // convert density to computing units
            rho_0 = density_to_sun_cu(Omega_m);

            // from mass and radius get density
            rho_c = (3*Ms)/(4*M_PI*pow(R, 3));

            // for star parameter R200 only for maximal integration
            R200 = 10 * R;
            break;
        }
        // NFW HALO: from concentration and mass computes density, R200
        case MOD_NFW:
        {
            // convert mass to computing units
            M200 = halo_mass_to_cu(M200_sun);

            // convert density to computing units
            rho_0 = density_to_halo_cu(Omega_m);
            
            // rho_0 parameter of the NFW halo
            rho_c = 200 * rho_0*c*(1 + c)*(1 + c);

            // from mass, rho_0 and concentration we get the scale radius
            R = pow(M200/(4*M_PI*rho_c)/(log(1+c)-c/(c+1)) , 1 / 3.0);

            // virial radius
            R200 = c*R;

			// Calculate 4 PI rho R3, check
			Ms = 4*M_PI*rho_c*pow(R, 3);
            break;
        }
    }
}

void chi_init()
{
	// first initialize phi_prefactor -- needed for potentials
	get_phi_prefactor();

    if (Ys)
    {
        // from the screening potential get R_eq
        R_eq = get_r_eq();
    }
    else
    {
        // from given R_eq we get the screening potential
        Ys = get_Ys();
    }

	// initialize chi_prefactor -- needd for mass
	get_chi_prefactor();

	chi_B = chi_bulk(rho_c + rho_0);
    m_inf = chi_mass(1);

	// debug
    BOOST_LOG_TRIVIAL(debug) << "Pot(0) = " << pot_new(0);
    BOOST_LOG_TRIVIAL(debug) << "Pot(R) = " << pot_new(R);
    BOOST_LOG_TRIVIAL(debug) << "Pot(R_eq) = " << pot_new(R_eq);
}

double chi_bulk(double rho){
	return pow(rho_0 / rho, 1 / (1 - n_chi));
}

double chi_bulk_der(double r){
	double chi_bulk_der;
	double x = r / R;
	
	switch (mod){
	case MOD_STAR: chi_bulk_der = 0; break;
	case MOD_NFW: chi_bulk_der = 1 / R*pow(rho_0 / rho_c, 1 / (1 - n_chi)) / (1 - n_chi)*pow(x*(1 + x)*(1 + x), n_chi / (1 - n_chi))*((1 + x)*(1 + x) + 2 * x*(1 + x)); break;
	}
	return chi_bulk_der;
}

double chi_bulk_laplace(double r){
	double x = r / R;

	double A = 1 / (R*R) / pow(1 - n_chi, 2);
	double B = pow(rho_0 / rho_c, 1 / (1 - n_chi));
	double C = pow(x, (2 * n_chi - 1) / (1 - n_chi));
	double D = pow(1 + x, 2 * n_chi / (1 - n_chi));
	double E = -n_chi*n_chi*(pow(1 + x, 2)) + 2 * n_chi*(1 + x + x*x) + 2 * x*(4 * x + 3);
	return A*B*C*D*E;
}

double chi_bulk_r(double r){
	return pow(rho_0 / rho_r(r), 1 / (1 - n_chi));
}


double V_eff_der(double r, double chi){
	if (chi <= 0) chi = chi_bulk_r(r);

	double v = chi_prefactor*(rho_r(r) - rho_0*pow(chi, n_chi - 1));

	if ((v != v) || (isinf(v))){
		return 0;
	}
	return v;
}

double V_eff_2nd_derr(double chi){
	// right computing units
	return chi_prefactor*rho_0*(1-n_chi)*pow(chi, n_chi - 2);
}

double chi_mass(double chi)
{
	return sqrt(V_eff_2nd_derr(chi));
}

double chi_mass2(double chi)
{
	return V_eff_2nd_derr(chi);
}

double chi_0_der(double r, double *chi){
	return chi[1];
}

double chi_0_der_lin(double r, double *chi, double a){
	return a*exp(-m_inf*r) / (r*r)*(m_inf*r + 1);
}

double chi_1_der(double r, double *chi){
	if (r == 0) return V_eff_der(r, chi[0]) / 3;
	else return V_eff_der(r, chi[0]) - 2 / r*chi[1];
}

double chi_0_der_NFW_lin(double r, double *chi){
	return chi[1];
}

double chi_1_der_NFW_lin(double r, double *chi){
	double m2 = chi_mass2(chi_bulk_r(r));
	if (m2 != m2) m2 = 0;
	if (r == 0) return (m2*chi[0] - chi_bulk_laplace(r)) / 3;
	return m2*chi[0] - chi_bulk_laplace(r) - 2 / r*chi[1];
}

bool is_integrate(double r, double *y, double r_max){
	if ((r_max != 0) && (r > r_max)) return false;
	if (r < R200) return true;
	return abs(1 - y[0]) > 1e-3;
}

void fce_min_r2(double s, double &t_0, double *y_0, double eps, double m){
	switch (mod){
	case 0:{ // stars
			   t_0 = s;
			   y_0[0] = chi_B*(1 + eps);
			   y_0[1] = chi_B*eps / s*(m*s / tanh(m*s) - 1);
	}
	case 1:{ // NFW halo
			   t_0 = s;
			   m = chi_mass(chi_bulk_r(t_0));
			   y_0[0] = chi_bulk_r(s)*(1 + eps);
			   y_0[1] = chi_bulk_der(s)*(1 + eps) + chi_bulk_r(s)*eps / s*(m*s / tanh(m*s) - 1);

	}
	}
}

void fce_min_chi_0(double s, double &t_0, double *y_0){
	switch (mod){
	case 0:{ // stars
			   t_0 = 0;
			   y_0[0] = s;
			   y_0[1] = 0;
	}
	case 1:{ // NFW halo
			   t_0 = 0;
			   y_0[0] = s;
			   y_0[1] = -(s - 1) / (2 * R)*(c + 1) / c; // Newtonian relationship between potential and derivative for NFW halo
	}
	}
}

double fce_max(double r, double *y){
	double dchi = y[0] - 1;
	if (dchi < 0) return (y[1] * r + dchi*(m_inf*r + 1));
	else return 10 * (y[1] * r + dchi*(m_inf*r + 1));
}

double shoot_meth_star(double err, double eps){
	double m = chi_mass(chi_B);
	double s1;

	if (R_eq > 10*R)
	{
		return R;
	}
	else if (R_eq < R)
	{
		s1 = R_eq;
	}
	else
	{
		s1 = R;
	}

	auto fce_min_star = bind(fce_min_r2, placeholders::_1, placeholders::_2, placeholders::_3, eps, m);
	auto integrate_star = bind(is_integrate, placeholders::_1, placeholders::_2, 1 / m_inf);
	t_function f_diff[] = { chi_0_der, chi_1_der };

	return shoot_meth(s1, 0.7, integrate_star, fce_min_star, fce_max, err, f_diff, 2, eps);
}

double shoot_meth_NFW(double err, double eps){
	double s1 = R_eq;
	auto fce_min_star = bind(fce_min_r2, placeholders::_1, placeholders::_2, placeholders::_3, eps, 0);
	auto integrate_star = bind(is_integrate, placeholders::_1, placeholders::_2, 1 / m_inf);
	t_function f_diff[] = { chi_0_der, chi_1_der };

	return shoot_meth(s1, 0.5, integrate_star, fce_min_star, fce_max, err, f_diff, 2, eps);
}

void slv_Chameleon(double r_max, double chi_pot_0, double chi_der_0, double err, chi_t& chi, double step){

	t_function chi_vec_eq[] = { chi_0_der, chi_1_der };
	double h = step;
	double r = 0;
	double chi_vec[2] = { chi_pot_0, chi_der_0 };
	int i = 0;

	chi.push_back(r, chi_vec[0], chi_vec[1]);
	while (r < r_max){
		Runge_Kutta_adap_step(r, chi_vec, h, err, chi_vec_eq, 2);
		if ((r - chi[0][i]) >= step){
			i++;
			chi.push_back(r, chi_vec[0], chi_vec[1]);
		}
	}
}

double get_chi_0_star(double err, double eps){
	double s1 = 1 + pot_star(0); // guess from the analytical solution
	auto integrate_star = bind(is_integrate, placeholders::_1, placeholders::_2, 1 / m_inf);
	t_function f_diff[] = { chi_0_der, chi_1_der };

	return shoot_meth(s1, 0.95, integrate_star, fce_min_chi_0, fce_max, err, f_diff, 2, 1*eps);
}

double get_chi_0_NFW(double err, double eps){
	double s1 = 1 + pot_NFW(0); // guess from the analytical solution
	auto integrate_star = bind(is_integrate, placeholders::_1, placeholders::_2, 1 / (1 * m_inf));
	t_function f_diff[] = { chi_0_der, chi_1_der };

	return shoot_meth(s1, 0.95, integrate_star, fce_min_chi_0, fce_max, err, f_diff, 2, eps*1*m_inf*1e3);
}

void slv_Chameleon_star(double r_max, double err, chi_t& chi, double &step){
	double h = step;
	double eps = 1e-2;
	t_function chi_vec_eq[] = { chi_0_der, chi_1_der };
	double chi_pot_0;
	if ((pot_new(0) + Ys) <= 0) chi_pot_0 = chi_B; // non-linear regime
	else chi_pot_0 = get_chi_0_star(err, eps); // linear regime
	double chi_der_0 = 0;
	double r = 0;
	int i = 0;
	double chi_vec[2] = { chi_pot_0, chi_der_0 };

	if (abs(chi_pot_0 - chi_B) / chi_B < eps){ // non-linear regime
		BOOST_LOG_TRIVIAL(info) << "Chameleon field in the screened regime. Getting <r2>...";
		double m = chi_mass(chi_B);
		double r2;
		if (R_eq>0){
			r2 = shoot_meth_star(err, eps);
			chi.push_back(0, chi_B*(1 + eps*m*r2 / sinh(m*r2)), 0);
		}
		else{
			r2 = 0;
			chi.push_back(0, chi_B*(1 + eps), 0);
		}

		BOOST_LOG_TRIVIAL(info) << "<r2> calculated as <" << r2 << ">. Filling analytical values from 0 to r2.";
		double shr_d_shr2, chr_d_shr2; // sinh(r) / sinhr(r2), cosh(r) / sinhr(r2)
		for (r = step; r < r2; r += step){
			i++;

			if (m*r>30){
				shr_d_shr2 = exp(m*(r - r2)); // approximation 2*sinh(x)=e^x
				chr_d_shr2 = shr_d_shr2;
			}
			else{
				shr_d_shr2 = sinh(m*r) / sinh(m*r2);
				chr_d_shr2 = cosh(m*r) / sinh(m*r2);
			}
			chi.push_back(r, chi_B*(1 + eps*r2 / r * shr_d_shr2), chi_B*eps * r2 / (r*r)*(m*r*chr_d_shr2 - shr_d_shr2));
		}
		if ((R - r2) / R < 1e-6){ // no-shell solution
			r = R;
			chi_vec[0] = chi_B + (1 - chi_B) / (m*R)*tanh(m*R);
			chi_vec[1] = (1 - chi_B)*(m_inf*R + 1) / (R*R)*(R - tanh(m*R) / (m*R));
		}
		else{
			r = r2;
			chi_vec[0] = chi_B*(1 + eps);
			chi_vec[1] = chi_B*eps / r2*(m*r2 / tanh(m*r2) - 1);
		}
	}
	else
	{
		i--;
		BOOST_LOG_TRIVIAL(info) << "Chameleon field in the linear regime.";
	}
	BOOST_LOG_TRIVIAL(info) << "Starting integration from r<" << r << "> to 1/m_inf<" << 1/m_inf << ">.";
	auto integrate_star = bind(is_integrate, placeholders::_1, placeholders::_2, 1 / m_inf);
	integrate_cout(r, chi_vec, integrate_star, err, chi_vec_eq, 2, chi, step, 1, i);

	BOOST_LOG_TRIVIAL(info) << "Starting integration from r<" << r << "> to r_max<" << r_max << ">.";
	double a = chi_vec[1] * r*r*exp(m_inf*r) / (m_inf*r + 1);
	chi_vec_eq[0] = bind(chi_0_der_lin, placeholders::_1, placeholders::_2, a);
	while (r < r_max){
		Runge_Kutta_adap_step(r, chi_vec, h, err, chi_vec_eq, 1);
		chi_vec[1] = a*exp(-m_inf*r) / (r*r)*(m_inf*r + 1);
		if ((r - chi[0][i]) >= step){
			step *= 2;
			i++;
			chi.push_back(r, chi_vec[0], chi_vec[1]);
		}
	}
}

void slv_Chameleon_star_cout(double r_max, double err, chi_t& chi, double &step){
	slv_Chameleon_star(r_max, err, chi, step);
	
	const size_t size = chi.size();

	// potential
	std::string file_name = param.out_opt.out_dir + "potential.dat";
	Ofstream File(file_name);
	double mlt = get_pot_mlt();
	
	BOOST_LOG_TRIVIAL(debug) << "Writing potential into file " << file_name;
	File << "#Radius\tNewtonian potential\tChameleon potential\n";
	File << "# r/r_s\t-Phi_N\t(phi_inf-phi)/(2*beta*M_PL)\n";
    File << std::scientific;
	for (size_t j = 0; j < size; j++){
		File << chi[0][j] / R << "\t" << -pot_new(chi[0][j]) << "\t" << (1 - chi[1][j])*mlt << endl;
	}
	File.close();

	// forces
	file_name = param.out_opt.out_dir + "forces.dat";
	File.open(file_name);
	mlt = get_force_mlt();
	
	BOOST_LOG_TRIVIAL(debug) << "Writing forces into file " << file_name;
	File << "#Radius\tChameleon force / Newtonian force\n";
	File << "# r/r_s\tF_phi/(2beta^2 F_N)\n";
    File << std::scientific;
	for (size_t j = 1; j < size; j++){
		File << chi[0][j] / R << "	" << mlt*chi[2][j] / (2 * beta*beta*force_new(chi[0][j])) << endl;
	}
}

void slv_Chameleon_NFW(double r_max, double err, chi_t& chi, double &step){
	double h = step;
	double eps = 1e-2;
	t_function chi_vec_eq[] = { chi_0_der, chi_1_der };
	double s;
	if ((pot_NFW(0) + Ys) <= 0) s = 0;
	else s = get_chi_0_NFW(err, eps);
	double r;
	int i = 0;
	double chi_vec[2];
	fce_min_chi_0(s, r, chi_vec);

	if (s == 0){ // screening regime
		BOOST_LOG_TRIVIAL(info) << "Chameleon field in the screened regime. Getting <r2>...";
		double m;
		double r2 = shoot_meth_NFW(err, eps);
		BOOST_LOG_TRIVIAL(info) << "<r2> calculated as <" << r2 << ">. Filling analytical values from 0 to r2.";
		chi.push_back(0, 0, 0);
		double shr_d_shr2, chr_d_shr2;
		for (r = step; r < r2; r += step){
			i++;
			m = chi_mass2(chi_bulk_r(r));

			if (m*r>30){ // approximation for m*r>>1
				shr_d_shr2 = exp(m*(r - r2)); // approximation 2*sinh(x)=e^x
				chr_d_shr2 = shr_d_shr2;
			}
			else{
				shr_d_shr2 = sinh(m*r) / sinh(m*r2);
				chr_d_shr2 = cosh(m*r) / sinh(m*r2);
			}
			
			chi.push_back(r, chi_bulk_r(r)*(1 + eps*r2 / r * shr_d_shr2),
							 (chi_bulk_der(r)*(1 + eps*r2 / r * shr_d_shr2) + 
							  chi_bulk_r(r)*eps * r2 / (r*r)*(m*r*chr_d_shr2 - shr_d_shr2)));
		}
		r = r2;
		m = chi_mass(chi_bulk_r(r));
		chi_vec[0] = chi_bulk_r(r)*(1 + eps);
		chi_vec[1] = chi_bulk_der(r)*(1 + eps) + chi_bulk_r(r)*eps / r*(m*r / tanh(m*r) - 1);
	}
	else
	{
		i--;
		BOOST_LOG_TRIVIAL(info) << "Chameleon field in the linear regime.";
	}

	BOOST_LOG_TRIVIAL(info) << "Starting integration from r<" << r << "> to 1/m_inf<" << 1/m_inf << ">.";
	auto integrate_star = bind(is_integrate, placeholders::_1, placeholders::_2, 1 / (1 * m_inf));
	integrate_cout(r, chi_vec, integrate_star, err, chi_vec_eq, 2, chi, step, 1, i);

	BOOST_LOG_TRIVIAL(info) << "Starting integration from r<" << r << "> to r_max<" << r_max << ">.";
	double a = chi_vec[1] * r*r*exp(m_inf*r) / (m_inf*r + 1);
	chi_vec_eq[0] = bind(chi_0_der_lin, placeholders::_1, placeholders::_2, a);
	while (r < r_max){
		Runge_Kutta_adap_step(r, chi_vec, h, err, chi_vec_eq, 1);
		//chi_vec[0] = 1 - a*exp(-m_inf*r) / r;
		chi_vec[1] = a*exp(-m_inf*r) / (r*r)*(m_inf*r + 1);
		if ((r - chi[0][i]) >= step){
			step *= 2;
			i++;
			chi.push_back(r, chi_vec[0], chi_vec[1]);
		}
	}
}

void slv_Chameleon_NFW_cout(double r_max, double err, chi_t& chi, double &step){
	slv_Chameleon_NFW(r_max, err, chi, step);
	const size_t size = chi.size();
	// potential
	std::string file_name = param.out_opt.out_dir + "potential.dat";
	Ofstream File(file_name);
	
	BOOST_LOG_TRIVIAL(debug) << "Writing potential into file " << file_name;
	File << "#Radius\tNewtonian potential\tChameleon potential\n";
	File << "# r/r_s\t-Phi_N\t(phi_inf-phi)/(2*beta*M_PL)\n";
    File << std::scientific;
	for (size_t j = 0; j < size; j++){
		File << chi[0][j] / R << "\t" << -pot_NFW(chi[0][j]) << "\t" << 1 - chi[1][j] << endl;
	}
	File.close();

	// forces
	file_name = param.out_opt.out_dir + "forces.dat";
	File.open(file_name);
	double mlt = get_force_mlt();
	
	BOOST_LOG_TRIVIAL(debug) << "Writing forces into file " << file_name;
	File << "#Radius\tChameleon force / Newtonian force\n";
	File << "# r/r_s\tF_phi/(2beta^2 F_N)\n";
    File << std::scientific;
	for (size_t j = 1; j < size; j++){
		File << chi[0][j] / R << "\t" << mlt*chi[2][j] / (2*beta*beta*force_NFW(chi[0][j])) << endl;
	}
}

double get_R200(){
	double A = rho_c / (200 * rho_0);
	double B = sqrt(27 * A*A + 4 * A);
	double C = 3 * sqrt(3)*B + 27 * A + 2;
	double D = pow(C / 2, 1 / 3.0);
	return R*(D + 1 / D - 2) / 3;
}

double pot_star(double r){
	if (r<R) return -2 * M_PI*rho_c*phi_prefactor*(R*R - r*r / 3.0);
	else return -4 * M_PI*rho_c*phi_prefactor*R*R*R / (3 * r);
}

double force_star(double r){
	if (r < R) return -4 * M_PI*rho_c*phi_prefactor*r / 3.0;
	else return -4 * M_PI*rho_c*phi_prefactor*R*R*R / (3 * r*r);
}

double rho_NFW(double r){
	double x = r / R;
	return rho_0 + rho_c / (x*pow(1 + x, 2));
}

double rho_star(double r){
	if (r < R) return rho_c;
	else return rho_0;
}

double rho_r(double r){
	double rho_r;
	switch (mod){
		case MOD_STAR: rho_r = rho_star(r); break;
		case MOD_NFW: rho_r = rho_NFW(r); break;
	}
	return rho_r;
}

double pot_NFW(double r){
	if (r == 0) return -Ms * phi_prefactor / R;
	double x = r / R;
	return -Ms * phi_prefactor / R*log(1 + x) / x;
}

double force_NFW(double r){
	double x = r / R;
	if (x < 1e-10) return -Ms * phi_prefactor / (R*R)*(1 / 2.0 - 2 * x / 3 + 3 * x*x / 4); // Taylor expansion to prevent great numerical errors
	return Ms * phi_prefactor / (R*R) * (x - (1 + x)*log(1 + x)) / (x*x*(1 + x));
}

double pot_new(double r){
	double pot_new;
	switch (mod){
		case MOD_STAR: pot_new = pot_star(r); break;
		case MOD_NFW: pot_new = pot_NFW(r); break;
	}
	return pot_new;
}

double force_new(double r){
	double force_new;
	switch (mod){
		case MOD_STAR: force_new = force_star(r); break;
		case MOD_NFW: force_new = force_NFW(r); break;
	}
	return force_new;
}

double r_eq_star(){
	double _r_eq;
	// linear case
	if (pot_star(0) + Ys > 0){
		double r = Ys / abs(pot_star(0)) - 1;
		_r_eq = -r * R;
	}
	// non-linear case, inside the star
	else if (pot_star(R) + Ys > 0){
		_r_eq = sqrt(3*(R*R - Ys / (2 * M_PI*phi_prefactor*rho_c)));
	}
	// non-linear case, outside the star
	else
	{
		_r_eq = 4 * M_PI*rho_c*phi_prefactor / 3 * pow(R, 3) / Ys;
	}
	return _r_eq;
}

double r_eq_NFW(){
	double _r_eq;
	// linear case
	if (pot_NFW(0) + Ys > 0){
		double r = Ys / abs(pot_NFW(0)) - 1;
		_r_eq = r / R;
	}
	else
	{
		double r1 = R * 2 * (1 - R*Ys / (Ms * phi_prefactor));
		double r2;

		get_x1_x2(r1, r2, pot_NFW, -Ys, 1.5);
		root_finder(r1, r2, pot_NFW, -Ys, Ys);

		_r_eq = r2;
	}
	return _r_eq;
}

double get_r_eq(){
	double _r_eq;
	switch (mod){
		case MOD_STAR: _r_eq = r_eq_star(); break;
		case MOD_NFW: _r_eq = r_eq_NFW(); break;
	}
	return _r_eq;
}

double get_rho_c(double _r_eq){
	double _rho_c;
	if (_r_eq >= 0)
	{
		_rho_c = -rho_c*Ys / pot_new(_r_eq);
	}
	else
	{
		double x = -_r_eq / R;
		_rho_c = get_rho_c(0) / (1 + x);
	}
	return _rho_c;
}

double get_Ys(){
	double _Ys;
	if (R_eq < 0) _Ys = (1 - R_eq / R)*abs(pot_new(abs(0)));
	else _Ys = abs(pot_new(R_eq));
	return _Ys;
}