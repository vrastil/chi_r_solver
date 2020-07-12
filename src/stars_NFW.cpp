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

static double V_EFF_DER_PREFACTOR;
static double V_EFF_2ND_DER_PREFACTOR;

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
            // convert mass to computing units, for star Ms = M200
            Ms = M200 = star_mass_to_cu(M200_sun);

            // convert density to computing units
            rho_0 = density_to_sun_cu(rho_0);

            // from mass and radius (in star units) get density
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
            rho_0 = density_to_halo_cu(rho_0);
            
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

	V_EFF_DER_PREFACTOR = get_v_eff_der_prefactor();
	V_EFF_2ND_DER_PREFACTOR = get_v_eff_2nd_der_prefactor();

	chi_B = chi_bulk(rho_c + rho_0);
    m_inf = chi_mass(1);
}

double chi_bulk(double rho){
	return pow(rho_0 / rho, 1 / (1 - n_chi));
}

double chi_bulk_der(double r){
	double x = r / R;
	switch (mod){
	case MOD_STAR:return 0;
	case MOD_NFW:return 1 / R*pow(rho_0 / rho_c, 1 / (1 - n_chi)) / (1 - n_chi)*pow(x*(1 + x)*(1 + x), n_chi / (1 - n_chi))*((1 + x)*(1 + x) + 2 * x*(1 + x));
	}

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

	double v = V_EFF_DER_PREFACTOR*(rho_r(r) - rho_0*pow(chi, n_chi - 1));

	if ((v != v) || (isinf(v))){
		return 0;
	}
	return v;
}

double V_eff_2nd_derr(double chi){
	// right computing units -- 1/kpc2 or 1/R_sun2
	return V_EFF_2ND_DER_PREFACTOR*pow(chi, n_chi - 2);
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
	/*
	double m2 = chi_mass2(chi[0]);
	if (m2 != m2) return 0;
	*/
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
	//	double mlt_err=1;
	//	if (dchi > 0) mlt_err *= 1000;
	if (dchi < 0) return (y[1] * r + dchi*(m_inf*r + 1));
	else return 10 * (y[1] * r + dchi*(m_inf*r + 1));
}

double shoot_meth_star(double err, double eps){
	double m = chi_mass(chi_B);
	double h = 1;
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
	// double m = chi_mass(chi_B);
	double h = 1;
	double s1 = R_eq;
	//	s1 = R200;
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
	double h = 1;
	double s1 = 1 + pot_star(0); // guess from the analytical solution
	//	auto fce_min_star = bind(fce_min_chi_0, placeholders::_1, placeholders::_2, placeholders::_3);
	auto integrate_star = bind(is_integrate, placeholders::_1, placeholders::_2, 1 / m_inf);
	t_function f_diff[] = { chi_0_der, chi_1_der };

	return shoot_meth(s1, 0.95, integrate_star, fce_min_chi_0, fce_max, err, f_diff, 2, 1*eps);

	/*
	double f1, f2;
	double sh, fh;
	double r;
	double chi_vec[2];
	//	double chi_B = 1*pow(rho_0 / (rho + rho_0), 1 / (1 - n_chi));
	double chi_00, mlt;
	if ((2*pot_new(0) + Ys) > 0){
	chi_00 = 1;
	mlt = 0.8;
	}
	else{
	chi_00 = chi_B;
	mlt = 2;
	}

	double a;
	r = 0;
	chi_vec[0] = chi_00*s1;
	chi_vec[1] = chi_der_min;


	while (r < 1 / m_inf){
	Runge_Kutta_adap_step(r, chi_vec, h, err, chi_vec_eq, 2);
	}

	f1 = (m_inf*r + 1) / r*(chi_vec[0] - 1) + chi_vec[1];
	//	f1 = chi_vec[0] - chi_max;
	f2 = f1;
	do{ // second guess of opposite sign
	r = 0;
	s1 = s2;
	f1 = f2;
	s2 *= mlt;
	chi_vec[0] = chi_00*s2;
	chi_vec[1] = chi_der_min;

	while (r < 1 / m_inf){
	Runge_Kutta_adap_step(r, chi_vec, h, err, chi_vec_eq, 2);
	if (chi_vec[0] <= 0){
	break;
	}
	}


	f2 = (m_inf*r + 1) / r*(chi_vec[0] - 1) + chi_vec[1];
	//	f2 = chi_vec[0] - chi_max;
	} while (f1*f2 > 0);

	return shoot_meth(0, r_max, chi_der_min, chi_max, s1*chi_00, s2*chi_00, err, chi_vec_eq, 1);
	*/
}

double get_chi_0_NFW(double err, double eps){
	double h = 1;
	double s1 = 1 + pot_NFW(0); // guess from the analytical solution
	//	auto fce_min_NFW = bind(fce_min_chi_0, placeholders::_1, placeholders::_2, placeholders::_3);
	auto integrate_star = bind(is_integrate, placeholders::_1, placeholders::_2, 1 / (1 * m_inf));
	t_function f_diff[] = { chi_0_der, chi_1_der };

	return shoot_meth(s1, 0.95, integrate_star, fce_min_chi_0, fce_max, err, f_diff, 2, eps*1*m_inf*1e3);
	/*
	double h = step;
	t_function chi_vec_eq[] = { chi_0_der, chi_1_der };
	double s1 = 1;
	double s2 = s1;
	double f1, f2;
	double sh, fh;
	double r;
	double chi_vec[2];
	//	double chi_B = 1*pow(rho_0 / (rho + rho_0), 1 / (1 - n_chi));
	double chi_00, mlt;
	mlt = 0.8;
	chi_00 = 1 + 0.8*pot_NFW(r_min); // 0.8 safe factor

	//	double a;
	r = 0;
	chi_vec[0] = chi_00*s1;
	chi_vec[1] = -(chi_vec[0]-1)/(2*R); //Newtonian


	while (r < 0.1 / m_inf){
	Runge_Kutta_adap_step(r, chi_vec, h, err, chi_vec_eq, 2);
	}

	f1 = (chi_vec[0] - 1) / r + chi_vec[1];
	f2 = f1;
	do{ // second guess of opposite sign
	r = 0;
	s1 = s2;
	f1 = f2;
	s2 *= mlt;
	chi_vec[0] = chi_00*s2;
	chi_vec[1] = -(chi_vec[0] - 1) / (2 * R); //Newtonian

	while (r < 0.1 / m_inf){
	Runge_Kutta_adap_step(r, chi_vec, h, err, chi_vec_eq, 2);
	}

	f2 = (chi_vec[0] - 1) / r + chi_vec[1];
	} while (f1*f2 > 0);

	bool y_0_acc = true;
	bool y_max_acc = true;

	while ((y_0_acc) && (y_max_acc)){
	sh = s2;
	s2 = (s2 + s1) / 2;

	r = 0;
	chi_vec[0] = chi_00*s2;
	chi_vec[1] = -(chi_vec[0] - 1) / (2 * R); //Newtonian
	while (r < 0.1 / m_inf){
	Runge_Kutta_adap_step(r, chi_vec, h, err, chi_vec_eq, 2);
	}

	fh = f2;
	f2 = (chi_vec[0] - 1) / r + chi_vec[1];

	if (f1*f2 > 0){
	s1 = sh;
	f1 = fh;
	}

	y_0_acc = (abs(s1 - s2) / s2) > 1e-6;
	if (chi_max == 0) y_max_acc = abs(f2) > 1e-5;
	else y_max_acc = abs(f2 / chi_max) > 1e-4;
	}
	return chi_00*s2;
	*/
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
	
	BOOST_LOG_TRIVIAL(debug) << "Writing potential into file " << file_name;
	File << "# r/r_s	-Phi_N	(phi_inf-phi)/(2*beta*M_PL)" << endl;
    File << std::scientific;
	for (size_t j = 0; j < size; j++){
		File << chi[0][j] / R << "	" << -pot_new(chi[0][j]) << "	" << 1 - chi[1][j] << endl;
	}
	File.close();

	// forces
	file_name = param.out_opt.out_dir + "forces.dat";
	File.open(file_name);
	
	BOOST_LOG_TRIVIAL(debug) << "Writing forces into file " << file_name;
	File << "# r/r_s	-F_N/M_PL	-F_\phi/M_PL" << endl;
    File << std::scientific;
	for (size_t j = 1; j < size; j++){
		File << chi[0][j] / R << "	" << chi[2][j] / (2 * beta*beta*force_new(chi[0][j])) << endl;
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
	File << "# r/r_s	-\Phi_N(r)	\delta\phi(r)/(2\beta\Mpl)" << endl;
    File << std::scientific;
	for (size_t j = 0; j < size; j++){
		File << chi[0][j] / R << "	" << -pot_NFW(chi[0][j]) << "	" << 1 - chi[1][j] << endl;
	}
	File.close();

	// forces
	file_name = param.out_opt.out_dir + "forces.dat";
	File.open(file_name);
	
	BOOST_LOG_TRIVIAL(debug) << "Writing forces into file " << file_name;
	File << "# r/r_s	-F_N	-F_\phi" << endl;
    File << std::scientific;
	for (size_t j = 1; j < size; j++){
		File << << chi[0][j] / R << "	" << chi[2][j] / (2*beta*beta*force_NFW(chi[0][j])) << endl;
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
	if (r<R) return -2 * M_PI*rho_c*(R*R - r*r / 3.0);
	else return -4 * M_PI*rho_c*R*R*R / (3 * r);
}

double force_star(double r){
	if (r < R) return -4 * M_PI*rho_c*r / 3.0;
	else return -4 * M_PI*rho_c*R*R*R / (3 * r*r);
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
	switch (mod){
	case 0:return rho_star(r);
	case 1:return rho_NFW(r);
	}
}

double pot_NFW(double r){
	if (r == 0) return -Ms / R;
	double x = r / R;
	return -Ms / R*log(1 + x) / x;
}

double force_NFW(double r){
	double x = r / R;
	if (x < 1e-10) return -Ms / (R*R)*(1 / 2.0 - 2 * x / 3 + 3 * x*x / 4); // Taylor expansion to prevent great numerical errors
	return Ms / (R*R) * (x - (1 + x)*log(1 + x)) / (x*x*(1 + x));
}

double pot_new(double r){
	switch (mod){
	case 0:return pot_star(r);
	case 1:return pot_NFW(r);
	}
}

double force_new(double r){
	switch (mod){
	case 0:return force_star(r);
	case 1:return force_NFW(r);
	}
}

double r_eq_star(){
	// linear case
	if (pot_star(0) + Ys > 0){
		double r = Ys / abs(pot_star(0)) - 1;
		return -r * R;
	}
	// non-linear case, inside the star
	else if (pot_star(R) + Ys > 0){
		return sqrt(3*(R*R - Ys / (2 * M_PI*rho_c)));
	}
	// non-linear case, outside the star
	else
	{
		return 4 * M_PI*rho_c / 3 * pow(R, 3) / Ys;
	}

}

double r_eq_NFW(){
	// linear case
	if (pot_NFW(0) + Ys > 0){
		double r = Ys / abs(pot_NFW(0)) - 1;
		return r / R;
	}
	else
	{
		double r1 = R * 2 * (1 - R*Ys / (Ms));
		double r2;

		get_x1_x2(r1, r2, pot_NFW, -Ys, 1.5);
		root_finder(r1, r2, pot_NFW, -Ys, Ys);

		return r2;
	}
}

double get_r_eq(){
	switch (mod){
	case 0:return r_eq_star();
	case 1:return r_eq_NFW();
	}
}

double get_rho_c(double _r_eq){
	if (_r_eq >= 0) return -rho_c*Ys / pot_new(_r_eq);
	double x = -_r_eq / R;
	return get_rho_c(0) / (1 + x);
}

double get_Ys(){
	if (R_eq < 0) return (1 - R_eq / R)*abs(pot_new(abs(0)));
	else return abs(pot_new(R_eq));
}