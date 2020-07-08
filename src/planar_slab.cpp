//////// PLANAR SLAB ///////////////

double pot_plan(double z){
	return phi_s*pow(1 + K*z, -2.0 / (n_pl - 2));
}

double V_pl_der(double phi){
	if (phi > 0) return n_pl*gamma_chi*pow(M_L, 4 - n_pl)*pow(phi, n_pl - 1);
	else return	-n_pl*gamma_chi*pow(M_L, 4 - n_pl)*pow(-phi, n_pl - 1);
}

bool is_integrate_pl(double r, double *y, double r_max){
	return (r < r_max);
}

void fce_min_phi_0(double s, double &t_0, double *y_0){
	t_0 = 0;
	y_0[0] = phi_s;
	y_0[1] = s;
}

double fce_max_pl(double r, double *y){
//	return y[0] - pot_plan(r);
//	return y[1];
 return	y[1] + 2*y[0] / ((n_pl - 2)*r);
}

double phi_0_der(double r, double *phi){
	return phi[1];
}

double phi_1_der(double r, double *phi){
	return m_pl*m_pl*phi[0];
}

double phi_pl_1_der(double r, double *phi){
	double tmp = V_pl_der(phi[0]);
	if (tmp != tmp){
		return 0;
	}
	return tmp;
}

double get_phi_der(double r_max, double phi_max, double phi_min, double err){
	double h = 1;
	t_function phi_vec_eq[] = { phi_0_der, phi_1_der };
	double s1 = 1;
	double s2 = s1;
	double f1, f2;
	double sh, fh;
	double r;
	double phi_vec[2];

	r = 0;
	phi_vec[0] = phi_min;
	phi_vec[1] = s1;
	while (r < r_max){ // first guess
		Runge_Kutta_adap_step(r, phi_vec, h, err, phi_vec_eq, 2);
		if (abs(phi_vec[0])>1e3) break;
	}
	f1 = phi_vec[0] - phi_max;

	do{ // second guess of opposite sign
		r = 0;
		s1 = s2;
		s2 -= 0.3;
		phi_vec[0] = phi_min;
		phi_vec[1] = s2;
		while (r < r_max){
			Runge_Kutta_adap_step(r, phi_vec, h, err, phi_vec_eq, 2);
			if (abs(phi_vec[0])>1e3) break;
		}
		f2 = phi_vec[0] - phi_max;
	} while (f1*f2 > 0);


	while (abs(s1 - s2) > 1e-10){
		sh = s2;
		//	s2 = (s1*f2 - s2*f1) / (f2-f1);
		//	s1 = sh;
		s2 = (s2 + s1) / 2;

		r = 0;
		phi_vec[0] = phi_min;
		phi_vec[1] = s2;
		while (r < r_max){
			Runge_Kutta_adap_step(r, phi_vec, h, err, phi_vec_eq, 2);
			if (abs(phi_vec[0])>1e3) break;
		}
		fh = f2;
		f2 = phi_vec[0] - phi_max;
		if (f1*f2 > 0){
			s1 = sh;
			f1 = fh;
		}
	}
	return s2;
}

double get_phi_pl_der(double r_max, double err){
	double h = 1;
	double s1;
	t_function phi_vec_eq[] = { phi_0_der, phi_pl_1_der };
	double r_mlt = 1;
	double f1;

	if (n_pl > 2) s1 = -K*phi_s * 2 / (n_pl - 2); // exact initial condition
	else {
		s1 = K*phi_s; // first guess is phi(0)/z_c
		do{
			auto integrate_pl = bind(is_integrate_pl, placeholders::_1, placeholders::_2, r_max*r_mlt);
			s1 = shoot_meth(s1, 0.9, integrate_pl, fce_min_phi_0, fce_max_pl, err, phi_vec_eq, 2, phi_s / r_max);
			f1 = integrate_fmax(s1, integrate_pl, fce_min_phi_0, fce_max_pl, err, phi_vec_eq, 2);
			r_mlt *= 0.9;
		} while ((abs(f1 / (phi_s / r_max)) > 1) && (r_mlt > 0.1)); // if the integration fails
	}	
	return s1;

	/*
	double s2 = 0;
	// s1 = -2.0 / (n_pl - 2)*sqrt(1.0 / 2 * (n_pl - 2)*(n_pl - 2)*gamma_chi*pow(M_L, 4 - n_pl))*pow(phi_s, n_pl / 2);
	double s1 = s2;
	double f1, f2;
	double sh, fh;
	double r;
	double phi_vec[2];

	r = 0;
	phi_vec[0] = phi_min;
	phi_vec[1] = s1*phi_min;
	while (r < r_max){ // first guess
		Runge_Kutta_adap_step(r, phi_vec, h, err, phi_vec_eq, 2);
		if (abs(phi_vec[0])>1e3*phi_min) break;
	}
	f1 = phi_vec[0] - phi_max;

	do{ // second guess of opposite sign
		r = 0;
		s1 = s2;
		s2 -= 1 * sgn(n_pl - 2);
		phi_vec[0] = phi_min;
		phi_vec[1] = phi_min*s2;
		while (r < r_max){
			Runge_Kutta_adap_step(r, phi_vec, h, err, phi_vec_eq, 2);
			if (abs(phi_vec[0])>1e3*phi_min) break;
		}
		f2 = phi_vec[0] - phi_max;
	} while (f1*f2 > 0);

	while ((abs((s1 - s2) / s2) > 1e-5) && (abs(f2) > 1e-5)){
		sh = s2;
		//	s2 = (s1*f2 - s2*f1) / (f2-f1);
		//	s1 = sh;
		s2 = (s2 + s1) / 2;

		r = 0;
		phi_vec[0] = phi_min;
		phi_vec[1] = phi_min*s2;
		while (r < r_max){
			Runge_Kutta_adap_step(r, phi_vec, h, err, phi_vec_eq, 2);
			if (abs(phi_vec[0])>1e3*phi_min) break;
		}
		fh = f2;
		f2 = phi_vec[0] - phi_max;
		if (f1*f2 > 0){
			s1 = sh;
			f1 = fh;
		}
	}
	return s2*phi_min;
	*/	
}

void slv_planar_slab(double r_max, double err, double **phi, double step, int &N_i){
	t_function phi_vec_eq[] = { phi_0_der, phi_pl_1_der };
	double h = 1;
	int i = -1;

	double phi_der_0 = get_phi_pl_der(r_max, err);
	double phi_vec[2] = { phi_s, phi_der_0 };

	double r = 0;
	auto integrate_pl = bind(is_integrate_pl, placeholders::_1, placeholders::_2, r_max);
	integrate_cout(r, phi_vec, integrate_pl, err, phi_vec_eq, 2, phi, step, 1, 0, i);
	N_i = i;
}

void slv_planar_slab(double r_max, double err, double **phi, double &step, int &N_i, int *cout_max){
	slv_planar_slab(r_max, err, phi, step, N_i);
	double err_pl=0;
	cout << "# z/z_c	\phi_A	\phi_N" << endl;
	for (int j = 0; j<*cout_max; j++){
		cout << phi[0][j] * K << "	" << abs(pot_plan(phi[0][j])) / Mp*1e6 << "	" << abs(phi[1][j]) / Mp *1e6 << endl;
		err_pl += pow((phi[1][j] - pot_plan(phi[0][j])) / pot_plan(phi[0][j]), 2);
	}
	if (phi[0][*cout_max-1] + step < r_max){
		for (double r_next = phi[0][*cout_max - 1] + step; r_next < r_max;r_next+=step){
			cout << r_next * K << "	" << abs(pot_plan(r_next)) / Mp*1e6 << endl;
		}
	}
	cout.precision(2);
	cout << "# Err = " << 100*sqrt(err_pl/ *cout_max) << "%";
}