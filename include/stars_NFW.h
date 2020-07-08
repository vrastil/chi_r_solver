//////// STARS AND NFW HALO //////////////////
#pragma once
double chi_bulk(double rho){
	return chi_0*pow(rho_0 / rho, 1 / (1 - n));
}

double chi_bulk_der(double r){
	double x = r / R;
	switch (mod){
	case 0:return 0;
	case 1:return 1 / R*chi_0*pow(rho_0 / rho_c, 1 / (1 - n)) / (1 - n)*pow(x*(1 + x)*(1 + x), n / (1 - n))*((1 + x)*(1 + x) + 2 * x*(1 + x));
	}

}

double chi_bulk_laplace(double r){
	double x = r / R;

	double A = chi_0 / (R*R) / pow(1 - n, 2);
	double B = pow(rho_0 / rho_c, 1 / (1 - n));
	double C = pow(x, (2 * n - 1) / (1 - n));
	double D = pow(1 + x, 2 * n / (1 - n));
	double E = -n*n*(pow(1 + x, 2)) + 2 * n*(1 + x + x*x) + 2 * x*(4 * x + 3);
	return A*B*C*D*E;
}

double chi_bulk_r(double r){
	return chi_0*pow(rho_0 / rho_r(r), 1 / (1 - n));
}


double V_eff_der(double r, double chi){
	double v = B / Mp*(rho_r(r) - rho_0*pow(chi_0 / abs(chi), 1 - n));
	if ((v != v) || (isinf(v))){
		return 0;
	}
	return v;
}

double V_eff_2nd_derr(double chi){
	return B / Mp*rho_0*(1 - n)*pow(chi / chi_0, n)*chi_0 / (chi*chi);
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
	double m2 = V_eff_2nd_derr(chi[0]);
	if (m2 != m2) return 0;
	*/
	return chi[1];
}

double chi_1_der_NFW_lin(double r, double *chi){
	double m2 = V_eff_2nd_derr(chi_bulk_r(r));
	if (m2 != m2) m2 = 0;
	if (r == 0) return (m2*chi[0] - chi_bulk_laplace(r)) / 3;
	return m2*chi[0] - chi_bulk_laplace(r) - 2 / r*chi[1];
}

bool is_integrate(double r, double *y, double r_max){
	if ((r_max != 0) && (r > r_max)) return false;
	if (r < R200) return true;
	return abs(chi_0 - y[0]) / chi_0 > 1e-2;
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
			   m = sqrt(V_eff_2nd_derr(chi_bulk_r(t_0)));
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
			   y_0[1] = -(s - chi_0) / (2 * R)*(c + 1) / c; // Newtonian relationship between potential and derivative for NFW halo
	}
	}
}

double fce_max(double r, double *y){
	double dchi = y[0] - chi_0;
	//	double mlt_err=1;
	//	if (dchi > 0) mlt_err *= 1000;
	if (dchi < 0) return (y[1] * r + dchi*(m_inf*r + 1));
	else return 10 * (y[1] * r + dchi*(m_inf*r + 1));
}

double shoot_meth_star(double err, double eps){
	double m = sqrt(V_eff_2nd_derr(chi_B));
	double h = 1;
	double r2a = R - Mp*(chi_0 - chi_B) / (B*rho_c*R*R)*R; // analytival value of r2
	double s1;
	if (r_eq<r2a) s1 = r_eq;
	else s1 = r2a;
	if (s1 > R) return R;
	//	s1 = R;
	auto fce_min_star = bind(fce_min_r2, placeholders::_1, placeholders::_2, placeholders::_3, eps, m);
	auto integrate_star = bind(is_integrate, placeholders::_1, placeholders::_2, 1 / m_inf);
	t_function f_diff[] = { chi_0_der, chi_1_der };

	return shoot_meth(s1, 0.7, integrate_star, fce_min_star, fce_max, err, f_diff, 2, chi_0*eps);

	/*
	double M = 4 * PI / 3 * rho_c*pow(R, 3);
	double a;


	while (r < 1 / m_inf){
	Runge_Kutta_adap_step(r, y, h, err, f, 2);
	if (y[0] <= 0){
	break;
	}
	}



	//		f1 = y[0] - y_max;
	f1 = (m_inf*r + 1) / r*(y[0] - chi_0) + y[1];
	f2 = f1;
	do{ // second guess of opposite sign
	guess1 = guess2;
	f1 = f2;
	guess2 /= 2;
	if (guess2 < 1e-10){
	guess2 = 0;
	r = guess2;
	y[0] = chi_B*(1 + eps);
	//y[1] = chameleon_der(r, m, eps);
	y[1] = 0;
	}
	else{
	r = guess2;
	y[0] = chi_B*(1 + eps);
	//y[1] = chameleon_der(r, m, eps);
	y[1] = chi_B*eps*(m*r - 1) / r;
	}

	while (r < 1 / m_inf){
	Runge_Kutta_adap_step(r, y, h, err, f, 2);
	if (y[0] <= 0){
	break;
	}
	}

	// f2 = y[0] - y_max;
	f2 = (m_inf*r + 1) / r*(y[0] - chi_0) + y[1];
	} while (f1*f2 > 0);

	while ((y_0_acc) && (y_max_acc)){
	guess_h = guess2;
	guess2 = (guess2 + guess1) / 2;

	r = guess2;
	y[0] = chi_B*(1 + eps);
	//y[1] = chameleon_der(r, m, eps);
	y[1] = chi_B*eps*(m*r - 1) / r;

	while (r < 1 / m_inf){
	Runge_Kutta_adap_step(r, y, h, err, f, 2);
	if (y[0] <= 0){
	break;
	}
	}

	fh = f2;
	// f2 = y[0] - y_max;
	f2 = (m_inf*r + 1) / r*(y[0] - chi_0) + y[1];
	if (f1*f2 > 0){
	guess1 = guess_h;
	f1 = fh;
	}

	y_0_acc = (abs(guess1 - guess2) / guess2) > 1e-6;
	if (y_max == 0) y_max_acc = abs(f2) > 1e-5;
	else y_max_acc = abs(f2 / y_max) > 1e-4;
	}

	return guess2;

	*/
}

double shoot_meth_NFW(double err, double eps){
	// double m = sqrt(V_eff_2nd_derr(chi_B));
	double h = 1;
	double s1 = r_eq;
	//	s1 = R200;
	auto fce_min_star = bind(fce_min_r2, placeholders::_1, placeholders::_2, placeholders::_3, eps, 0);
	auto integrate_star = bind(is_integrate, placeholders::_1, placeholders::_2, 1 / m_inf);
	t_function f_diff[] = { chi_0_der, chi_1_der };

	return shoot_meth(s1, 0.5, integrate_star, fce_min_star, fce_max, err, f_diff, 2, chi_0*eps);
}

void slv_Chameleon(double r_max, double chi_pot_0, double chi_der_0, double err, double **chi, double step, int &N_i, int reg){

	t_function chi_vec_eq[] = { chi_0_der, chi_1_der };
	double h = step;
	double r = 0;
	double chi_vec[2] = { chi_pot_0, chi_der_0 };
	int i = 0;

	chi[0][i] = r;
	chi[1][i] = chi_vec[reg];
	while (r < r_max){
		Runge_Kutta_adap_step(r, chi_vec, h, err, chi_vec_eq, 2);
		if ((r - chi[0][i]) >= step){
			i++;
			chi[0][i] = r;
			chi[1][i] = chi_vec[reg];
			if (reg == 1) chi[1][i] *= -B / Mp;
		}
	}
	N_i = i;
}

double get_chi_0_star(double err, double eps){
	double h = 1;
	double s1 = chi_0 + 1 * 2 * B*Mp*pot_star(0); // guess from the analytical solution
	//	auto fce_min_star = bind(fce_min_chi_0, placeholders::_1, placeholders::_2, placeholders::_3);
	auto integrate_star = bind(is_integrate, placeholders::_1, placeholders::_2, 1 / m_inf);
	t_function f_diff[] = { chi_0_der, chi_1_der };

	return shoot_meth(s1, 0.95, integrate_star, fce_min_chi_0, fce_max, err, f_diff, 2, chi_0*eps);

	/*
	double f1, f2;
	double sh, fh;
	double r;
	double chi_vec[2];
	//	double chi_B = chi_0*pow(rho_0 / (rho + rho_0), 1 / (1 - n));
	double chi_00, mlt;
	if ((2*pot_new(0) + Ys) > 0){
	chi_00 = chi_0;
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

	f1 = (m_inf*r + 1) / r*(chi_vec[0] - chi_0) + chi_vec[1];
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


	f2 = (m_inf*r + 1) / r*(chi_vec[0] - chi_0) + chi_vec[1];
	//	f2 = chi_vec[0] - chi_max;
	} while (f1*f2 > 0);

	return shoot_meth(0, r_max, chi_der_min, chi_max, s1*chi_00, s2*chi_00, err, chi_vec_eq, 1);
	*/
}

double get_chi_0_NFW(double err, double eps){
	double h = 1;
	double s1 = chi_0 + 2 * B*Mp*pot_NFW(0); // guess from the analytical solution
	//	auto fce_min_NFW = bind(fce_min_chi_0, placeholders::_1, placeholders::_2, placeholders::_3);
	auto integrate_star = bind(is_integrate, placeholders::_1, placeholders::_2, 1 / (1 * m_inf));
	t_function f_diff[] = { chi_0_der, chi_1_der };

	return shoot_meth(s1, 0.95, integrate_star, fce_min_chi_0, fce_max, err, f_diff, 2, eps*chi_0*m_inf*1e3);
	/*
	double h = step;
	t_function chi_vec_eq[] = { chi_0_der, chi_1_der };
	double s1 = 1;
	double s2 = s1;
	double f1, f2;
	double sh, fh;
	double r;
	double chi_vec[2];
	//	double chi_B = chi_0*pow(rho_0 / (rho + rho_0), 1 / (1 - n));
	double chi_00, mlt;
	mlt = 0.8;
	chi_00 = chi_0 + 0.8*2 * B*Mp*pot_NFW(r_min); // 0.8 safe factor

	//	double a;
	r = 0;
	chi_vec[0] = chi_00*s1;
	chi_vec[1] = -(chi_vec[0]-chi_0)/(2*R); //Newtonian


	while (r < 0.1 / m_inf){
	Runge_Kutta_adap_step(r, chi_vec, h, err, chi_vec_eq, 2);
	}

	f1 = (chi_vec[0] - chi_0) / r + chi_vec[1];
	f2 = f1;
	do{ // second guess of opposite sign
	r = 0;
	s1 = s2;
	f1 = f2;
	s2 *= mlt;
	chi_vec[0] = chi_00*s2;
	chi_vec[1] = -(chi_vec[0] - chi_0) / (2 * R); //Newtonian

	while (r < 0.1 / m_inf){
	Runge_Kutta_adap_step(r, chi_vec, h, err, chi_vec_eq, 2);
	}

	f2 = (chi_vec[0] - chi_0) / r + chi_vec[1];
	} while (f1*f2 > 0);

	bool y_0_acc = true;
	bool y_max_acc = true;

	while ((y_0_acc) && (y_max_acc)){
	sh = s2;
	s2 = (s2 + s1) / 2;

	r = 0;
	chi_vec[0] = chi_00*s2;
	chi_vec[1] = -(chi_vec[0] - chi_0) / (2 * R); //Newtonian
	while (r < 0.1 / m_inf){
	Runge_Kutta_adap_step(r, chi_vec, h, err, chi_vec_eq, 2);
	}

	fh = f2;
	f2 = (chi_vec[0] - chi_0) / r + chi_vec[1];

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

void slv_Chameleon_star(double r_max, double err, double **chi, double &step, int &N_i, int reg){
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
		double m = sqrt(V_eff_2nd_derr(chi_B));
		double r2;
		if (r_eq>0){
			r2 = shoot_meth_star(err, eps);
			chi[0][i] = 0;
			chi[1][i] = (1 - reg)*chi_B*(1 + eps*m*r2 / sinh(m*r2));
		}
		else{
			r2 = 0;
			chi[0][i] = 0;
			chi[1][i] = (1-reg)*chi_B*(1 + eps);
		}

		double shr_d_shr2, chr_d_shr2; // sinh(r) / sinhr(r2), cosh(r) / sinhr(r2)

		for (r = step; r < r2; r += step){
			i++;
			chi[0][i] = r;

			if (m*r>30){
				shr_d_shr2 = exp(m*(r - r2)); // approximation 2*sinh(x)=e^x
				chr_d_shr2 = shr_d_shr2;
			}
			else{
				shr_d_shr2 = sinh(m*r) / sinh(m*r2);
				chr_d_shr2 = cosh(m*r) / sinh(m*r2);
			}

			if (reg == 0) chi[1][i] = chi_B*(1 + eps*r2 / r * shr_d_shr2);
			else chi[1][i] = -B / Mp*chi_B*eps * r2 / (r*r)*(m*r*chr_d_shr2 - shr_d_shr2);
		}
		if ((R - r2) / R < 1e-6){ // no-shell solution
			r = R;
			chi_vec[0] = chi_B + (chi_0 - chi_B) / (m*R)*tanh(m*R);
			chi_vec[1] = (chi_0 - chi_B)*(m_inf*R + 1) / (R*R)*(R - tanh(m*R) / (m*R));
		}
		else{
			r = r2;
			chi_vec[0] = chi_B*(1 + eps);
			chi_vec[1] = chi_B*eps / r2*(m*r2 / tanh(m*r2) - 1);
		}
	}
	else i--;
	//	auto fce_min_star = bind(fce_min, placeholders::_1, placeholders::_2, placeholders::_3, 0, -1);
	auto integrate_star = bind(is_integrate, placeholders::_1, placeholders::_2, 1 / m_inf);
	//	double s=0;
	double mlt = 1;
	if (reg == 1) mlt *= -B / Mp;

	integrate_cout(r, chi_vec, integrate_star, err, chi_vec_eq, 2, chi, step, mlt, reg, i);

	double a = chi_vec[1] * r*r*exp(m_inf*r) / (m_inf*r + 1);
	chi_vec_eq[0] = bind(chi_0_der_lin, placeholders::_1, placeholders::_2, a);
	while (r < r_max){
		Runge_Kutta_adap_step(r, chi_vec, h, err, chi_vec_eq, 1);
		chi_vec[1] = a*exp(-m_inf*r) / (r*r)*(m_inf*r + 1);
		if ((r - chi[0][i]) >= step){
			step *= 2;
			i++;
			chi[0][i] = r;
			chi[1][i] = chi_vec[reg] * mlt;
		}
	}
	N_i = i;
}

void slv_Chameleon_star(double r_max, double err, double **chi, double &step, int &N_i, int reg, int *cout_max){
	slv_Chameleon_star(r_max, err, chi, step, N_i, reg);
	switch (reg){
	case 0:{
			   cout << "# r/r_s	-Phi_N	(phi_inf-phi)/(2*B*Mp)" << endl;
			   for (int j = 0; j<*cout_max; j++){
				   cout << chi[0][j] / R << "	" << -pot_new(chi[0][j]) << "	" << 1 / (2 * B * Mp)*(chi_0 - chi[1][j]) << endl;
			   }
			   break;
	}
	case 1:{
			   cout << "# r/r_s	-F_N/Mp	-F_\phi/Mp" << endl;
			   for (int j = 1; j < *cout_max; j++){
				   //	   cout << chi[0][j]/R << "	" << -force_new(chi[0][j]) / Mp << "	" << -1 / (2 * B*B) * chi[1][j] / Mp << endl;
				   //	   cout << chi[0][j] / R << "	" << -force_new(chi[0][j]) / Mp << "	" << -chi[1][j] / Mp << endl;
				   cout << chi[0][j] / R << "	" << chi[1][j] / (2 * B*B*force_new(chi[0][j])) << endl;
			   }
			   break;
	}
	}
}

void slv_Chameleon_NFW(double r_max, double err, double **chi, double &step, int &N_i, int reg){
	double h = step;
	double eps = 1e-2;
	t_function chi_vec_eq[] = { chi_0_der, chi_1_der };
	double s;
	if ((pot_NFW(0) + Ys) <= 0) s = 0;
	else s = get_chi_0_NFW(err, eps);
	// chi_0 + 2 * B*Mp*pot_NFW(0);// 
	double r;
	int i = 0;
	double chi_vec[2];
	fce_min_chi_0(s, r, chi_vec);

	if (s == 0){ // screening regime
		double m;
		double r2 = shoot_meth_NFW(err, eps);
		chi[0][i] = 0;
		chi[1][i] = 0;
		double shr_d_shr2, chr_d_shr2;
		for (r = step; r < r2; r += step){
			i++;
			if_low_memory(i, h_N, h_re, chi, 2);
			m = sqrt(V_eff_2nd_derr(chi_bulk_r(r)));
			chi[0][i] = r;

			if (m*r>30){ // approximation for m*r>>1
				shr_d_shr2 = exp(m*(r - r2)); // approximation 2*sinh(x)=e^x
				chr_d_shr2 = shr_d_shr2;
			}
			else{
				shr_d_shr2 = sinh(m*r) / sinh(m*r2);
				chr_d_shr2 = cosh(m*r) / sinh(m*r2);
			}

			if (reg == 0) chi[1][i] = chi_bulk_r(r)*(1 + eps*r2 / r * shr_d_shr2);
			else chi[1][i] = -B / Mp*(chi_bulk_der(r)*(1 + eps*r2 / r * shr_d_shr2) + chi_bulk_r(r)*eps * r2 / (r*r)*(m*r*chr_d_shr2 - shr_d_shr2));
		}
		r = r2;
		m = sqrt(V_eff_2nd_derr(chi_bulk_r(r)));
		chi_vec[0] = chi_bulk_r(r)*(1 + eps);
		chi_vec[1] = chi_bulk_der(r)*(1 + eps) + chi_bulk_r(r)*eps / r*(m*r / tanh(m*r) - 1);
	}
	else i--;

	//	auto fce_min_star = bind(fce_min, placeholders::_1, placeholders::_2, placeholders::_3, 0, -1);
	auto integrate_star = bind(is_integrate, placeholders::_1, placeholders::_2, 1 / (1 * m_inf));
	//	double s=0;
	double mlt = 1;
	if (reg == 1) mlt *= -B / Mp;

	integrate_cout(r, chi_vec, integrate_star, err, chi_vec_eq, 2, chi, step, mlt, reg, i);

	double a = chi_vec[1] * r*r*exp(m_inf*r) / (m_inf*r + 1);
	// double b = -(chi_vec[0] - chi_0)*r*exp(m_inf*r);
	//	a = b;
	//	a = (a+10*b)/11;

	//if (reg == 0) a = -(chi_vec[0] - chi_0)*r*exp(m_inf*r);
	//else  a = chi_vec[1] * r*r*exp(m_inf*r) / (m_inf*r + 1);
	/*
	while (r < r_max){
	r += step;
	step *= 2;
	i++;
	if_low_memory(i, h_N, h_re, chi, 2);
	chi[0][i] = r;
	if (reg == 0) chi[1][i] = chi_0 - a*exp(-m_inf*r) / r;
	else chi[1][i] = -mlt*a*exp(-m_inf*r) / (r*r)*(m_inf*r + 1);
	}
	*/

	chi_vec_eq[0] = bind(chi_0_der_lin, placeholders::_1, placeholders::_2, a);
	while (r < r_max){
		Runge_Kutta_adap_step(r, chi_vec, h, err, chi_vec_eq, 1);
		//chi_vec[0] = chi_0 - a*exp(-m_inf*r) / r;
		chi_vec[1] = a*exp(-m_inf*r) / (r*r)*(m_inf*r + 1);
		if ((r - chi[0][i]) >= step){
			step *= 2;
			i++;
			if_low_memory(i, h_N, h_re, chi, 2);
			chi[0][i] = r;
			chi[1][i] = chi_vec[reg] * mlt;
		}
	}

	N_i = i;
}

void slv_Chameleon_NFW(double r_max, double err, double **chi, double &step, int &N_i, int reg, int *cout_max){
	slv_Chameleon_NFW(r_max, err, chi, step, N_i, reg);
	switch (reg){
	case 0:{
			   cout << "# r/r_s	-\Phi_N(r)	\delta\phi(r)/(2\beta\Mpl)" << endl;
			   for (int j = 0; j<*cout_max; j++){
				   cout << chi[0][j] / R << "	" << -pot_NFW(chi[0][j]) << "	" << 1 / (2 * B*Mp)*(chi_0 - chi[1][j]) << endl;
				   //	   cout << chi[0][j] << "	" << chi[1][j] << "	" << chi_0 << endl;
				   //	   cout << chi[0][j] << "	" << chi_0 - chi[1][j] << "	" << chi_0 - chi_bulk_r(chi[0][j]) << endl;
			   }
			   break;
	}
	case 1:{
			   double r, Mr;
			   cout << "# r/r_s	-F_N	-F_\phi" << endl;
			   for (int j = 0; j < *cout_max; j++){
				   // cout << chi[0][j] / R << "	" << -force_NFW(chi[0][j])/Mp << "	" << -chi[1][j]/Mp << endl;
				   cout << chi[0][j] / R << "	" << chi[1][j] / (2*B*B*force_NFW(chi[0][j])) << endl;
				   //  cout << chi[0][j] / R*10 << "	" << sqrt(chi[0][j] * abs(force_NFW(chi[0][j])))*3e5 << "	" << sqrt(chi[0][j] * abs(force_NFW(chi[0][j]) + chi[1][j]))*3e5 << endl; // [km/s] vs [kpc]
				   //   cout << chi[0][j] << "	" << -chi[1][j] << "	" << B / Mp*chi_bulk_der(chi[0][j]) << endl;
				 //  r = chi[0][j];
				//   Mr = abs(force_NFW(r))*r*r / G;
				//   cout << r / 1.57E26 << "	" << Mr / Mp*4.34e-9 / 1.9891e42 << "	" << (Mr + r*r / G*abs(chi[1][j])) / Mp*4.34e-9 / 1.9891e42 << endl;
			   }
			   break;
	}
	}
}