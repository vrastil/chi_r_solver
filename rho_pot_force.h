// rho_pot_force.h
//
#pragma once
double get_R200(){
	double A = rho_c / (200 * rho_0);
	double B = sqrt(27 * A*A + 4 * A);
	double C = 3 * sqrt(3)*B + 27 * A + 2;
	double D = pow(C / 2, 1 / 3.0);
	return R*(D + 1 / D - 2) / 3;
}

double pot_star(double r){
	if (r<R) return -2 * PI*G*rho_c*(R*R - r*r / 3.0);
	else return -4 * PI*G*rho_c*R*R*R / (3 * r);
}

double force_star(double r){
	if (r < R) return -4 * PI*G*rho_c*r / 3.0;
	else return -4 * PI*G*rho_c*R*R*R / (3 * r*r);
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
//	double phi_0 = Ms*G / R / (c + 1);
//	if (r == 0) return -Ms*G / R *(1 - 1 / (c + 1));
	if (r == 0) return -Ms*G / R;
	double x = r / R;
	return -Ms*G / R*log(1 + x) / x;
//	if (x < c) return -Ms*G / R*(log(1 + x) / x - 1 / (1 + c));
	// else return -Ms*G / (R*x)*(log(1 + c) - c / (c + 1));
}

double force_NFW(double r){
	double x = r / R;
//	if (x == 0) return -Ms*G / (2 * R*R);
	if (x < 1e-10) return -Ms*G / (R*R)*(1 / 2.0 - 2 * x / 3 + 3 * x*x / 4); // Taylor expansion to prevent great numerical errors
	return Ms*G / (R*R) * (x - (1 + x)*log(1 + x)) / (x*x*(1 + x));
//	if (x<c) return Ms*G / (R*R) * (x - (1 + x)*log(1 + x)) / (x*x*(1 + x));
//	else return -Ms*G / (R*R)*(log(1 + c) - c / (c + 1)) / (x*x);
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
	if (pot_star(0) + Ys > 0){
		double r = Ys / abs(pot_star(0)) - 1;
		return -r * R;
	}
	else if (pot_star(R) + Ys > 0){
		return sqrt(3*(R*R - Ys / (2 * PI*G*rho_c)));
	}
	else{
		return 4 * PI*G*rho_c / 3 * pow(R, 3) / Ys;
	}

}

double r_eq_NFW(){
	if (pot_NFW(0) + Ys > 0){
		double r = Ys / abs(pot_NFW(0)) - 1;
		return r / R;
	}
	double r1 = R * 2 * (1 - R*Ys / (Ms*G));
	double r2;

	get_x1_x2(r1, r2, pot_NFW, -Ys, 1.5);
	root_finder(r1, r2, pot_NFW, -Ys, Ys);

	return r2;
}

double get_r_eq(){
	switch (mod){
	case 0:return r_eq_star();
	case 1:return r_eq_NFW();
	}
}

double get_rho_c(double r_eq){
	if (r_eq >= 0) return -rho_c*Ys / pot_new(r_eq);
	double x = -r_eq / R;
	return get_rho_c(0) / (1 + x);
}

double get_Ys(){
	if (r_eq < 0) return (1 - r_eq / R)*abs(pot_new(abs(0)));
	else return abs(pot_new(r_eq));
}