#pragma once

#include <cmath>
#include <cstdint>
#include <string>

typedef uint8_t mod_t;


// const double G = 6.67384E-11; // gravitation constatnt
static const double G = 6.71E-57; // gravitation constatnt, eV
static const double M_PL = sqrt(1 / (8 * M_PI*G));// Planck mass

#define M_SUN_TO_EV(m) (m*M_PL/4.34e-9 * 1.9891e42)


class Parameters
{
public:
    void print_info() const;
    void init();

    /********************
     * GENERIC
     ********************/
    struct Generic_t {
        mod_t mod;
        mod_t reg_halo;
        mod_t reg_rho;

        double err;
        double step;
        double r_max;
        size_t h_N;
        size_t h_re;
    } generic;
    /********************/

    /********************
     * SPATIAL PROPERTIES
     ********************/
    struct Spatial_t {
        double c;
        double R;
        double R_eq;
        double R200;
        double rho_c;
        double rho_0;
        double Ms;
        double M200; // mass of the halo in eV
        double M200_sun; // mass of the halo in solar units / 1e12
    } spatial;
    /********************/

    /********************
     * CHAMELEON
     ********************/
    struct Chi_opt_t {
        double n;
        double beta;
        double Ys;

        double chi_0;
        double chi_B;
        double m_inf;
    } chi_opt;
    /********************/


    /********************
     * INPUT / OUTPUT
     ********************/
    struct {
        std::string out_dir;
    } out_opt;
    /********************/

    // int h_re;

};

extern Parameters param;

enum class status_t {
    status_OK = 0,
    cmd_opt_HELP,
    cmd_opt_BAD_CFG,
};

status_t handle_cmd_line(int ac, const char* const av[]);