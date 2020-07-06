#pragma once

class Parameters
{
public:
    void print_info() const;

    struct {
        int mod; // 0 = star, 1 = NFW
    } generic;

    struct {
        double n;
    } chi_opt;

    struct {
        std::string out_dir;
    } out_opt;

    // double chi_0;
    // double Ys;
    // double R; // scale radius
    // double R200; // virial radius
    // double c; // concentration parameter
    // double rho_c;
    // double Ms;
    // double r_eq;
    // double rho_0;
    // double chi_B;
    // double m_inf;

    // int h_N;
    // int h_re;

};

extern Parameters param;

#define p_mod param.generic.mod
#define p_out_dir param.out_opt.out_dir


enum class status_t {
    status_OK = 0,
    cmd_opt_HELP,
    cmd_opt_BAD_CFG,
};

status_t handle_cmd_line(int ac, const char* const av[]);