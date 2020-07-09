// main.cpp : Defines the entry point for the console application.
//

#include "integrator.hpp"
#include "io.hpp"
#include "stars_NFW.hpp"

int main(int argc, char* argv[])
{
	status_t status = handle_cmd_line(argc, argv);

	if (status == status_t::status_OK)
	{
		param.init();
		param.print_info();

	}

	return int(status);



	// switch (mod){
	// case 0: { // star
	// 			R200 = 10 * R;
	// 			break;
	// }
	// case 1: { // NFW halo
	// 			R200 = get_R200();
	// 			break;
	// }
	// }
	// c = R200 / R;
	
	// M200 = Ms*(log(1+c)-c/(1+c));
	// M200_sun = M200/Mp*4.34e-9/1.9891e42;

	// double r_0=R*pow(rho_c/rho_0,1/3.0);
	// double r_max;

	// // double M = 4 * PI / 3 * rho_c*pow(R, 3);
	// //////////////////////////////////////

	// ////////////////////////////////////////////////////
	// ////////ALLOCATION, REALLOCATION ///////////////////
	// double **chi = (double**)malloc(2 * sizeof(void*));

	// for (int i = 0; i < 2; i++){
	// 	chi[i] = (double*)malloc(h_N * sizeof(double));
	// }
	// ////////////////////////////////////////////////////
	// int N_i;


	// switch (mod){
	// case 0: { // star
	// 			slv_Chameleon_star(r_max, err, chi, step, N_i, reg, &N_i);
	// 			break;
	// }
	// case 1: { // NFW halo
	// 			slv_Chameleon_NFW(r_max, err, chi, step, N_i, reg, &N_i);
	// 			break;
	// }
	// // case 2: { // planar slab
	// // 			slv_planar_slab(r_max, err, chi, step, N_i, &N_i);
	// // 			break;
	// // }
	// }

	// /*
	// for (double step_j = R / 100; step_j < 10/m_inf; step_j *= 1.3){
	// 	cout << step_j / R << "	" << (1-rho_0/rho_NFW(step_j))*(1 / step_j + 2 / (R*(1 + step_j / R)))/Mp << "	" << sqrt(V_eff_2nd_derr(chi_bulk_r(step_j)))/Mp << endl;
	// }
	// */
	// /*
	// for (int j = 1; j<40 * R; j++){
	// 	cout << j / R << "	" << B / Mp*(chi_0 - chi_bulk_laplace(j)) << "	" << -pot_NFW(j) << endl;
	// 	}
	// */
	// /*
	// for (int j = 1; j<40 * R; j++){
	// //	cout << j  << "	" << Mp / B*chi_bulk_laplace(j) / rho_NFW(j) << "	" << V_eff_2nd_derr(chi_bulk_r(j)) << endl;
	// 	cout << j << "	" << chi_bulk_r(j)  << "	" << chi_bulk_der(j)<< endl;
	// }
	// */


	return 0;
}

