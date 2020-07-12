// runge_integ_root.h
//

#include <math.h>
#include "integrator.hpp"

static void Runge_Kutta_adap_step(double &t_0, double *y_0, double &h_in, double h_min, double h_max, double atol,
                           double rtol, t_function *f, int dim, int &N_iter);
static void Runge_Kutta_adap_step(double &t_0, double *y_0, double &h_in, double h_min, double h_max, double err, t_function *f, int dim);
static void Runge_Kutta_adap_step(double &t_0, double *y_0, double &h_in, double h_min, double h_max, double err, t_function *f, int dim, int &N_iter);


static void integrate(double s, double &t, double *y, std::function<bool(double, double*)>t_max, std::function<void(double, double &, double*)>fce_min,
               double err, t_function *f_diff, int dim);
static void root_finder(double &x1, double &x2, std::function<double(double)> fce_x, double scale);
static void root_finder(double &x1, double &x2, std::function<double(double)> fce_x, double f_eq, double x_rtol, double x_atol, double f_atol);
static void get_x1_x2(double &x1, double &x2, std::function<double(double)> fce, double mlt);

static void shift(double *y, double *y_0, double *shift, double h, int dim){
	for (int i = 0; i < dim; i++){
		y[i] = y_0[i] + h*shift[i];
	}
}

void fce(double *y, t_function *f, double t_0, double *x_0, int dim){
	for (int i = 0; i < dim; i++){
		y[i] = f[i](t_0, x_0);
	}
}

void Runge_Kutta_step(double &t_0, double *y_0, double h, t_function *f, int dim){

	double *k1 = new double[dim];
	double *k2 = new double[dim];
	double *k3 = new double[dim];
	double *k4 = new double[dim];
	double *y_h = new double[dim];

	fce(k1, f, t_0, y_0, dim);

	shift(y_h, y_0, k1, h / 2, dim);
	fce(k2, f, t_0 + h / 2, y_h, dim);

	shift(y_h, y_0, k2, h / 2, dim);
	fce(k3, f, t_0 + h / 2, y_h, dim);

	shift(y_h, y_0, k3, h, dim);
	fce(k4, f, t_0 + h, y_h, dim);

	t_0 += h;
	shift(y_0, y_0, k1, h / 6, dim);
	shift(y_0, y_0, k2, h / 3, dim);
	shift(y_0, y_0, k3, h / 3, dim);
	shift(y_0, y_0, k4, h / 6, dim);

	delete[] k1;
	delete[] k2;
	delete[] k3;
	delete[] k4;
	delete[] y_h;
}

void Runge_Kutta_adap_step(double &t_0, double *y_0, double &h_in, double h_min, double h_max, double atol, double rtol, t_function *f, int dim, int &N_iter){
	if (N_iter>30) Runge_Kutta_step(t_0, y_0, h_in, f, dim); // prevents from an infinite loop, h_out = h_in *1e-21
	else {
		if ((h_in > h_max) && h_max != 0) h_in = h_max;
		if (h_in < h_min) h_in = h_min;

		double err_h = 0;
		double scale;
		double dh;
		double * yh = new double[dim];
		double * y1 = new double[dim];
		for (int i = 0; i < dim; i++){
			yh[i] = y_0[i];
			y1[i] = y_0[i];
		}

		Runge_Kutta_step(t_0, y1, h_in, f, dim); // one step of h
		t_0 -= h_in;
		Runge_Kutta_step(t_0, y_0, h_in / 2, f, dim); // two steps of h/
		Runge_Kutta_step(t_0, y_0, h_in / 2, f, dim);
		
		for (int i = 0; i < dim; i++){
			scale = atol + rtol*abs(y_0[i]);
			err_h += pow((y_0[i] - y1[i]) / scale, 2);
		}
		err_h = sqrt(err_h / dim);
		if (err_h == 0) dh = 10;
		else dh = 0.9*pow(1 / err_h, 1.0 / 5);
		if (dh > 10) dh = 10; // maximal increase of the stepsize by a factor 10 or
		if (dh < 1 / 5.0) dh = 1 / 5.0; // decrease by a factor 5
		if (err_h < 1){
			if (N_iter == 0) h_in *= dh; // no increase if previous step fails
		}
		else{
			t_0 -= h_in; // go back
			for (int i = 0; i < dim; i++){
				y_0[i] = yh[i];
			}
			delete[] y1;
			delete[] yh;
			h_in *= dh; // adjust h
			N_iter++;
			Runge_Kutta_adap_step(t_0, y_0, h_in, h_min, h_max, atol, rtol, f, dim, N_iter); // reiterate
		}
	}
}

void Runge_Kutta_adap_step(double &t_0, double *y_0, double &h_in, double h_min, double h_max, double err, t_function *f, int dim, int &N_iter){
	Runge_Kutta_adap_step(t_0, y_0, h_in, h_min, h_max, err, err, f, dim, N_iter);
}
void Runge_Kutta_adap_step(double &t_0, double *y_0, double &h_in, double err, t_function *f, int dim){
	int N_iter = 0;
	Runge_Kutta_adap_step(t_0, y_0, h_in, 0, 0, err, f, dim, N_iter);
}

void Runge_Kutta_adap_step(double &t_0, double *y_0, double &h_in, double h_min, double h_max, double err, t_function *f, int dim){
	int N_iter = 0;
	Runge_Kutta_adap_step(t_0, y_0, h_in, h_min, h_max, err, f, dim, N_iter);
}

void integrate(double s, double &t, double *y, std::function<bool(double, double*)>t_max, std::function<void(double, double &, double*)>fce_min, double err, t_function *f_diff, int dim){
	double h = 1;
	fce_min(s, t, y);
	while (t_max(t, y)){
		Runge_Kutta_adap_step(t, y, h, err, f_diff, dim);
	}
}

void integrate_cout(double &t, double *y, std::function<bool(double, double*)>t_max, double err, t_function *f_diff, int dim, chi_t& chi, double step, double mlt, int &i){
	double h = step;
	double h_max = step;
	i++;
	chi.push_back(t, y[0], y[1]);
	while (t_max(t, y)){
		Runge_Kutta_adap_step(t, y, h, 0,h_max, err, f_diff, dim);
		if (t > param.spatial.R) {
			h_max = param.spatial.R;
		}
		if (t > param.spatial.R) { 
			h_max = 0;
		}
		if ((t - chi[0][i]) >= step){
			if(t>param.spatial.R) step *= 1.5;
			i++;
			chi.push_back(t, y[0], y[1]);
		}
	}
}

double integrate_fmax(double s, std::function<bool(double, double*)>t_max, std::function<void(double, double &, double*)>fce_min, t_function fce_max, double err, t_function *f_diff, int dim){
	double *y = (double*)malloc(dim*sizeof(double));
	double t = 0;
	integrate(s, t, y, t_max, fce_min, err, f_diff, dim);
	return fce_max(t, y);
}

void root_finder(double &x1, double &x2, std::function<double(double)> fce_x, double f_eq, double x_rtol, double x_atol, double f_atol){
	// return x such as fce(x) - f_eq = 0
	// fce(x1)*fce(x2) < 0;
	double f1, f2, fh, xh;
	bool f_acc = true;
	bool x_acc = true;
	double scale;
	f1 = fce_x(x1) - f_eq;
	f2 = fce_x(x2) - f_eq;
	while ((f_acc) && (x_acc)){
		fh = f2;
		xh = x2;
		x2 = (x2 + x1) / 2;
		f2 = fce_x(x2) - f_eq;
		if (f1*f2 > 0){
			x1 = xh;
			f1 = fh;
		}
		scale = x_atol + x2*x_rtol;
		x_acc = (abs(x1 - x2) / scale) > 1;
		scale = f_atol;
		f_acc = abs(f2 / scale) > 1;
	}
	if (abs(f1) > abs(f2)){
		xh = x2;
		x2 = x1; // x1 <=> |f1| < |f2|
		x1 = xh;
	}	
}
void root_finder(double &x1, double &x2, std::function<double(double)> fce_x, double f_eq, double scale){
	root_finder(x1, x2, fce_x, f_eq, 1e-8, 0, scale*1e-8);
}

void root_finder(double &x1, double &x2, std::function<double(double)> fce_x, double scale){
	root_finder(x1, x2, fce_x, 0, scale);
}

void get_x1_x2(double &x1, double &x2, std::function<double(double)> fce, double f_eq, double mlt){
	// input x1 and mlt
	// return x1 and x2 with (fce(x_1) - f_eq)*(fce(x_2) - f_eq) < 0
	double xh = x1;
	int N_h = 0;
	double f1, f2;

	x2 = x1;
	f2 = fce(x2) - f_eq;
	do {
		x1 = x2;
		x2 *= mlt;
		f1 = f2;
		f2 = fce(x2) - f_eq;
		N_h++;
	} while ((f1*f2 > 0) && (N_h < 100));
	if (N_h < 100) return;

	x2 = xh;
	f2 = fce(x2) - f_eq;
	mlt = 1 / mlt;
	N_h = 0;
	do {
		x1 = x2;
		x2 *= mlt;
		f1 = f2;
		f2 = fce(x2) - f_eq;
		N_h++;
	} while ((f1*f2 > 0) && (N_h < 100));
}

void get_x1_x2(double &x1, double &x2, std::function<double(double)> fce, double mlt){
	get_x1_x2(x1, x2, fce, 0, mlt);
}

double shoot_meth(double s1, double mlt, std::function<bool(double, double*)>t_max, std::function<void(double, double &, double*)>fce_min,
	t_function fce_max, double err, t_function *f_diff, int dim, double scale){
	// s1, mlt -- first guess and a multiplier to obtain other guesses
	// initial conditions is specified in fce_min with parameter s
	// integrate till t_max is true
	// boundary conditions at t_max are specified in fce_max, return 0 if achieved
	double s2 = 0;

	auto fce_s = bind(integrate_fmax, std::placeholders::_1, t_max, fce_min, fce_max, err, f_diff, dim);
	get_x1_x2(s1, s2, fce_s, mlt);
	root_finder(s1, s2, fce_s, scale);
	return s1;
}

double shoot_meth(double s1, double mlt, std::function<bool(double, double*)>t_max, std::function<void(double, double &, double*)>fce_min, t_function fce_max, double err, t_function *f_diff, int dim){
	return shoot_meth(s1, mlt, t_max, fce_min, fce_max, err, f_diff, dim, 0);
}