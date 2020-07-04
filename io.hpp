#pragma once

class Parameters
{
    int mod; // 0 = star, 1 = NFW
    double n;
    double chi_0;
    double Ys;
    double R; // scale radius
    double R200; // virial radius
    double c; // concentration parameter
    double rho_c;
    double Ms;
    double r_eq;
    double rho_0;
    double chi_B;
    double m_inf;

    int h_N;
    int h_re;

};

void handle_cmd_line(int ac, const char* const av[], Parameters& param);