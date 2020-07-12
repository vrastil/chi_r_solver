/**
 * @brief command line arguments manipulation
 * 
 * @file core_cmd.cpp
 * @author Michal Vrastil
 * @date 2018-07-11
 */
#define BOOST_LOG_DYN_LINK 1

#include <boost/log/trivial.hpp>
#include <boost/any.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include <cstddef>
#include <iomanip>
#include <iostream>
#include <fstream>

#include "io.hpp"
#include "stars_NFW.hpp"

namespace po = boost::program_options;

#define MAX_MEMORY_SIZE (std::size_t)(10000)

/****************************//**
 * PUBLIC FUNCTIONS DEFINITIONS *
 ********************************/

Parameters param;

status_t handle_cmd_line(int ac, const char* const av[]){
    std::string config_file;
    status_t status = status_t::status_OK;

    // GENERIC
    po::options_description generic("Generic options");
    generic.add_options()
        ("help,h", "produce this help message")
        ("mod", po::value<mod_t>(&param.integration.mod)->default_value(0), "Mode; 0: star, 1: NFW")
        ("err", po::value<double>(&param.integration.err)->default_value(1E-8), "absolute / relative tolerances")
        ;


    // SPATIAL PROPERTIES
    po::options_description spatial("Spatial properties");
    spatial.add_options()
        ("c", po::value<double>(&param.spatial.c)->default_value(5), "concentration of the halo")
        ("R", po::value<double>(&param.spatial.R)->default_value(25), "characteristic scale - radius of star in solar units / scale radius in kpc for NFW")
        ("M200_sun", po::value<double>(&param.spatial.M200_sun)->default_value(10), "mass of the star/halo in solar units, for halo with additional factor of 1E12")
        ("rho_0", po::value<double>(&param.spatial.rho_0)->default_value(4E-11), "background density")
        ;
        
    // CHAMELEON
    po::options_description chi_opt("Chameleon parameters options");
    chi_opt.add_options()
        ("n", po::value<double>(&param.chi_opt.n)->default_value(0.5), "chameleon power-law potential exponent,0 < n < 1")
        ("beta", po::value<double>(&param.chi_opt.beta)->default_value(1/sqrt(6.)), "chameleon coupling constant")
        ("Ys", po::value<double>(&param.chi_opt.Ys)->default_value(1E-5), "screening potential")
        ("R_eq", po::value<double>(&param.spatial.R_eq)->default_value(1), "distance where screening potential equals gravitational")
        ;
    
    // INPUT / OUTPUT
    po::options_description config_output("Output options");
    config_output.add_options()
        ("out_dir,o", po::value<std::string>(&param.out_opt.out_dir)->default_value("output/"), "output folder name")
        ("print_par", po::value<bool>(&param.out_opt.print_par)->default_value(true), " to print parameters at the start")
        ("config,c", po::value<std::string>(&config_file)->default_value("input.cfg"), "configuration file name (optional)")
        ;

    // ALL IN ONE
    po::options_description cmdline_options("\n" "PROJECT" " v" "PROJECT_VERSION");//< store all normal parameters
    cmdline_options.add(generic).add(spatial).add(chi_opt).add(config_output);

    // PARSE PARAMETERS
    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, cmdline_options), vm);
    po::notify(vm);

    // HANDLE SPECIAL CASES
    if (vm.count("help")) {
        std::cout << std::setprecision(3) << cmdline_options;
        status = status_t::cmd_opt_HELP;
    }
    else if (vm.count("config")) {
        std::ifstream ifs(config_file.c_str());
        if (ifs){
        BOOST_LOG_TRIVIAL(debug) << "Using configuration options defined in file " << config_file << " (given command line options have higher priority).";
            po::store(po::parse_config_file(ifs, cmdline_options), vm);
            notify(vm);
        } else{
            BOOST_LOG_TRIVIAL(warning) << "Cannot open config file '" << config_file << "'. Using default and command lines values.";
        }
    }

    return status;
}

void Parameters::print_info() const
{
    BOOST_LOG_TRIVIAL(info) <<
        "Parameters given through command line options, through configuration file, or derived:\n"
        "\nGeneric options:\n"
        "\tmod = " << param.integration.mod << "\n"

        "\nIntegration options:\n"
        "\terr = " << param.integration.err << "\n"
        "\tstep = " << param.integration.step << "\n"
        "\tr_max = " << param.integration.r_max << "\n"
        "\th_N = " << param.integration.h_N << "\n"
        "\th_re = " << param.integration.h_re << "\n"

        "\nSpatial parameters:\n"
        "\tc = " << param.spatial.c << "\n"
        "\tR = " << param.spatial.R << "\n"
        "\tR200 = " << param.spatial.R200 << "\n"
        "\trho_c = " << param.spatial.rho_c << "\n"
        "\trho_0 = " << param.spatial.rho_0 << "\n"
        "\tM200 = " << param.spatial.M200 << "\n"
        "\tMs = " << param.spatial.Ms << "\n"
        "\tM200_sun = " << param.spatial.M200_sun << "\n"

        "\nChameleon parameters:\n"
        "\tn = " << param.chi_opt.n << "\n"
        "\tbeta = " << param.chi_opt.beta << "\n"
        "\tYs = " << param.chi_opt.Ys << "\n"
        "\tR_eq = " << param.spatial.R_eq << "\n"
        "\tchi_B = " << param.chi_opt.chi_B << "\n"
        "\tm_inf = " << param.chi_opt.m_inf << "\n"

        "\nOutput options:\n"
        "\tout dir" << param.out_opt.out_dir;
}


void Parameters::init()
{
    // SPATIAL PROPERTIES
    spatial_init();

    // CHI PROPERTIES
    chi_init();

    // INTEGRATION
    integration.step = spatial.R / 10;
	integration.r_max = 10 / chi_opt.m_inf;
	integration.h_N = std::min((size_t)(spatial.R200/integration.step+log(integration.r_max/integration.step)), MAX_MEMORY_SIZE);
	integration.h_re = (size_t)(log(integration.r_max/integration.step));

    // ALLOCATION
    for (auto vec : integration.chi)
    {
        vec.reserve(integration.h_N);
    }

}