#pragma once

#include <cstddef>
#include <cmath>
#include <string>
#include <array>
#include <vector>
#include <fstream>

typedef std::size_t mod_t;

#define MOD_STAR 0
#define MOD_NFW  1

class chi_t
{
private:
    typedef typename std::array<std::vector<double>, 3> array_type;
    typedef typename array_type::iterator iterator;
    typedef typename array_type::const_iterator const_iterator;
    array_type m_vec;

public:
    iterator begin() noexcept { return m_vec.begin(); }
    const_iterator cbegin() const noexcept { return m_vec.cbegin(); }
    iterator end() noexcept { return m_vec.end(); }
    const_iterator cend() const noexcept { return m_vec.cend(); }
    void reserve(size_t n) { for (auto vec : m_vec) vec.reserve(n); }
    void push_back(const double& r, const double& pot, const double& force)
    {
        m_vec[0].push_back(r);
        m_vec[1].push_back(pot);
        m_vec[2].push_back(force);
    }
    size_t size() { return m_vec[0].size(); }
    std::vector<double>& operator[] (size_t n) { return m_vec[n]; }
    const std::vector<double>& operator[] (size_t n) const { return m_vec[n]; }
};

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
        // chi[0] -- radius
        // chi[1] -- potential
        // chi[2] -- derivative (force WITHOUT -beta/M_PL factor)
        chi_t chi;
    } integration;

    /********************
     * SPATIAL PROPERTIES
     ********************/
    struct {
        double c;
        double R;
        double R_eq;
        double R200;
        double Omega_m;
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


/**
 * class  Ofstream handles opening and closing files, has 16MB buffer for output 
 */

class Ofstream : public std::ofstream
{
public:
    Ofstream(std::string file_name);
    char* buf;
    ~Ofstream();
};

/**
 * class  Ifstream handles opening and closing files
 */

class Ifstream : public std::ifstream
{
public:
    Ifstream(std::string file_name);
    ~Ifstream();
};

std::string currentDateTime();
std::string std_out_dir();

void create_dir(const std::string& out_dir);
void remove_dir(const std::string &out_dir);
void remove_all_files(const std::string &out_dir);