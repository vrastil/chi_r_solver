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
	// }

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


	return 0;
}

