// NFW_Chameleon.cpp : Defines the entry point for the console application.
//

#include <math.h>
#include "limits.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <functional>
using namespace std;

#include "constants.h"

#include "integrator.hpp"
#include "io.hpp"

#include "planar_slab.h"
#include "rho_pot_force.h"
#include "stars_NFW.h"

int main(int argc, char* argv[])
{
	if (status_t::status_OK == handle_cmd_line(argc, argv))
	{
		return 0;
	}
	else
	{
		return 0;
	}


	cout.precision(15);
	////////// PARAMETERS ///////////
	R = 1.57E27; // characteristic scale - radius of star / scale radius for NFW
//	R = 1.57E17;
	r_eq =5*R; // valid for r_eq < 0 (linear regime), rho_c decreasing as 1/(1-r_eq/R)
	mod = 2; // 0...star, 1,...NFW, 2...planar slab
	int reg = 1; // 0...potential, 1...forces
	int reg_halo = 0; // 0..c,M200=(rho_c, rho_0, r_s), 1..rho_c(c, rho_0) r_s (M200)
	int reg_rho = 0; // 0..rho_c(R_eq, Ys), 1..R_eq(rho_c, Ys), 2..Ys(R_eq, rho_c)
	double err = 1e-8; // absolute / relative tolerances
	double step = R/10;

	/////////////////////////////////
	Ms = 4 * PI*rho_c*pow(R, 3);	 
	n = 1 / 2.0;
	Ys = 1e-15;
	// Ys = 4.58335e-13;
	double M200_sun;  // mass of the halo in solar units / 1e12
	double M200; // mass of the halo in eV
	
	switch (reg_halo){
	case 1:{
			   rho_0 = 4e-11;
			   c = 15; //16.4397
			   rho_c = 200 * rho_0*c*(1 + c)*(1 + c);
			   M200_sun = 10.0;
			   M200 = M200_sun * Mp/4.34e-9 * 1.9891e42;
			   R = pow(M200/(4*PI*rho_c)/(log(1+c)-c/(c+1)) , 1 / 3.0);
			   break;
	}
	default:{
				rho_c = 3e-005;
			//	rho_c = 5e7;
				rho_0 = rho_c*1e-6;
			}
	}
	Ms = 4 * PI*rho_c*pow(R, 3);
	switch (reg_rho){
	case 0: {
				for (int i = 0; i < 3; i++){ // iteration to get the right r_eq
					rho_c = get_rho_c(r_eq);
					rho_0 = rho_c*1e-6;
					Ms = 4 * PI*rho_c*pow(R, 3);
					switch (mod){
					case 0: { // star
								R200 = 3 * R;
								break;
					}
					case 1: { // NFW halo
								R200 = get_R200();
								break;
					}
					case 2:{ // Planar slab
							   R200 = 1 / K;
							   break;
					}
					}
					c = R200 / R;
				}
				break;
	}
	case 1:{
			   r_eq = get_r_eq();
			  
			   break;
	}
	case 2:{
			   Ys = get_Ys();
	}
	}
	switch (mod){
	case 0: { // star
				R200 = 10 * R;
				break;
	}
	case 1: { // NFW halo
				R200 = get_R200();
				break;
	}
	// case 2:{ // Planar slab
			   R200 = 1 / K;
			   step = 1e-2 / K;
			   break;
	// }
	}
	c = R200 / R;
	
	M200 = Ms*(log(1+c)-c/(1+c));
	M200_sun = M200/Mp*4.34e-9/1.9891e42;
	chi_0 = 2 * B*Mp*Ys;
	chi_B = chi_bulk(rho_c + rho_0);
	m_inf = sqrt(V_eff_2nd_derr(chi_0));

	double r_0=R*pow(rho_c/rho_0,1/3.0);
	double r_max;

	r_max = 10 / m_inf;
	r_max = 1e1 / K;
	// double M = 4 * PI / 3 * rho_c*pow(R, 3);
	//////////////////////////////////////

	////////////////////////////////////////////////////
	////////ALLOCATION, REALLOCATION ///////////////////
	double **chi = (double**)malloc(2 * sizeof(void*));

	h_N = (int)(R200/step+log(r_max/step));
	h_re = (int)(log(r_max/step));
	for (int i = 0; i < 2; i++){
		chi[i] = (double*)malloc(h_N * sizeof(double));
	}
	////////////////////////////////////////////////////
	int N_i;


	switch (mod){
	case 0: { // star
				slv_Chameleon_star(r_max, err, chi, step, N_i, reg, &N_i);
				break;
	}
	case 1: { // NFW halo
				slv_Chameleon_NFW(r_max, err, chi, step, N_i, reg, &N_i);
				break;
	}
	// case 2: { // planar slab
	// 			slv_planar_slab(r_max, err, chi, step, N_i, &N_i);
	// 			break;
	// }
	}

	/*
	for (double step_j = R / 100; step_j < 10/m_inf; step_j *= 1.3){
		cout << step_j / R << "	" << (1-rho_0/rho_NFW(step_j))*(1 / step_j + 2 / (R*(1 + step_j / R)))/Mp << "	" << sqrt(V_eff_2nd_derr(chi_bulk_r(step_j)))/Mp << endl;
	}
	*/
	/*
	for (int j = 1; j<40 * R; j++){
		cout << j / R << "	" << B / Mp*(chi_0 - chi_bulk_laplace(j)) << "	" << -pot_NFW(j) << endl;
		}
	*/
	/*
	for (int j = 1; j<40 * R; j++){
	//	cout << j  << "	" << Mp / B*chi_bulk_laplace(j) / rho_NFW(j) << "	" << V_eff_2nd_derr(chi_bulk_r(j)) << endl;
		cout << j << "	" << chi_bulk_r(j)  << "	" << chi_bulk_der(j)<< endl;
	}
	*/


	cin.get();
	return 0;
}

