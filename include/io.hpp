#pragma once

#include <cstddef>
#include <cmath>
#include <string>
#include <array>
#include <vector>

typedef std::size_t mod_t;

#define MOD_STAR 0
#define MOD_NFW  1

class Parameters
{
public:
    void print_info() const;
    void init();

    /********************
     * GENERIC / INTEGRATION
     ********************/
    struct {
        mod_t mod;
        double err;
        double step;
        double r_max;
        size_t h_N;
        size_t h_re;
        // chi[0] -- radius
        // chi[1] -- potential
        // chi[2] -- derivative (force WITHOUT -beta/M_PL factor)
        std::array<std::vector<double>, 3> chi;
    } integration;

    /********************
     * SPATIAL PROPERTIES
     ********************/
    struct {
        double c;
        double R;
        double R_eq;
        double R200;
        double rho_c;
        double rho_0;
        double M200; // mass of the halo in eV
        double Ms; // 4 PI rho R3 in [kpc]
        double M200_sun; // mass of the halo in solar units / 1e12
    } spatial;

    /********************
     * CHAMELEON
     ********************/
    struct {
        double n;
        double beta;
        double Ys;

        double chi_B;
        double m_inf;
    } chi_opt;

    /********************
     * INPUT / OUTPUT
     ********************/
    struct {
        std::string out_dir;
        bool print_par;
    } out_opt;
};

extern Parameters param;

enum class status_t {
    status_OK = 0,
    cmd_opt_HELP,
    cmd_opt_BAD_CFG,
};

status_t handle_cmd_line(int ac, const char* const av[]);