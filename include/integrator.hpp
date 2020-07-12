// runge_integ_root.h
//
#pragma once

#include <functional>
#include "io.hpp"

typedef std::function<double(double, double*)> t_function;

void Runge_Kutta_adap_step(double &t_0, double *y_0, double &h_in, double err, t_function *f, int dim);
void integrate_cout(double &t, double *y, std::function<bool(double, double*)>t_max, double err, t_function *f_diff, int dim, chi_t& chi,
               double step, double mlt, int &i);
double shoot_meth(double s1, double mlt, std::function<bool(double, double*)>t_max, std::function<void(double, double &, double*)>fce_min,
	t_function fce_max, double err, t_function *f_diff, int dim, double scale);
void root_finder(double &x1, double &x2, std::function<double(double)> fce_x, double f_eq, double scale);
void get_x1_x2(double &x1, double &x2, std::function<double(double)> fce, double f_eq, double mlt);