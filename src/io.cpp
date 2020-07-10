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

#include <iomanip>
#include <iostream>
#include <fstream>

#include "io.hpp"
#include "stars_NFW.hpp"

namespace po = boost::program_options;

#define MOD_STAR 0
#define MOD_NFW  1

#define MOD_M200_DERIVED 0
#define MOD_RHOC_DERIVED 1

#define MOD_CHI_RHOC_DERIVED 0
#define MOD_CHI_REQ_DERIVED  1
#define MOD_CHI_YS_DERIVED   2

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
        ("mod", po::value<mod_t>(&param.generic.mod)->default_value(0), "Mode; 0: star, 1: NFW")
        ("reg_halo", po::value<mod_t>(&param.generic.reg_halo)->default_value(0), "Mode; 0: c, M200=(rho_c, rho_0, R_s), 1: rho_c(c, rho_0), R_s (M200)")
        ("reg_rho", po::value<mod_t>(&param.generic.reg_rho)->default_value(0), "Mode; 0: rho_c(R_eq, Ys); 1: R_eq(rho_c, Ys); 2: Ys(R_eq, rho_c)")
        ("err", po::value<double>(&param.integration.err)->default_value(1E-8), "absolute / relative tolerances")
        ;


    // SPATIAL PROPERTIES
    po::options_description spatial("Spatial properties");
    spatial.add_options()
        ("c", po::value<double>(&param.spatial.c)->default_value(5), "concentration of the halo")
        ("R", po::value<double>(&param.spatial.R)->default_value(25), "characteristic scale - radius of star / scale radius for NFW")
        ("R200", po::value<double>(&param.spatial.R200)->default_value(200), "virial radius of the halo, radius where rho = 200 * rho_background")
        ("M200_sun", po::value<double>(&param.spatial.M200_sun)->default_value(10), "mass of the halo in solar units / 1e12")
        ("rho_c", po::value<double>(&param.spatial.rho_c)->default_value(3E-5), "bulk density of the star / NFW")
        ("rho_0", po::value<double>(&param.spatial.rho_0)->default_value(4E-11), "background densit")
        ;
        
    // CHAMELEON
    po::options_description chi_opt("Chameleon parameters options");
    chi_opt.add_options()
        ("n", po::value<double>(&param.chi_opt.n)->default_value(0.5), "chameleon power-law potential exponent,0 < n < 1")
        ("beta", po::value<double>(&param.chi_opt.beta)->default_value(1/sqrt(6.)), "chameleon coupling constant")
        ("Ys", po::value<double>(&param.chi_opt.Ys)->default_value(1E-5), "screening potential?")
        ;
    
    // INPUT / OUTPUT
    po::options_description config_output("Output options");
    config_output.add_options()
        ("out_dir,o", po::value<std::string>(&param.out_opt.out_dir)->default_value("output/"), "output folder name")
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
        "\tmod = " << param.generic.mod << "\n"
        "\treg_halo = " << param.generic.reg_halo << "\n"
        "\treg_rho = " << param.generic.reg_rho << "\n"

        "\nIntegration options:\n"
        "\terr = " << param.integration.err << "\n"
        "\tstep = " << param.integration.step << "\n"
        "\tr_max = " << param.integration.r_max << "\n"
        "\th_N = " << param.integration.h_N << "\n"
        "\th_re = " << param.integration.h_re << "\n"

        "\nSpatial parameters:\n"
        "\tc = " << param.spatial.c << "\n"
        "\tR = " << param.spatial.R << "\n"
        "\tR_eq = " << param.spatial.R_eq << "\n"
        "\tR200 = " << param.spatial.R200 << "\n"
        "\trho_c = " << param.spatial.rho_c << "\n"
        "\trho_0 = " << param.spatial.rho_0 << "\n"
        "\tMs = " << param.spatial.Ms << "\n"
        "\tM200 = " << param.spatial.M200 << "\n"
        "\tM200_sun = " << param.spatial.M200_sun << "\n"

        "\nChameleon parameters:\n"
        "\tn = " << param.chi_opt.n << "\n"
        "\tbeta = " << param.chi_opt.beta << "\n"
        "\tYs = " << param.chi_opt.Ys << "\n"
        "\tchi_0 = " << param.chi_opt.chi_0 << "\n"
        "\tchi_B = " << param.chi_opt.chi_B << "\n"
        "\tm_inf = " << param.chi_opt.m_inf << "\n"

        "\nOutput options:\n"
        "\tout dir" << param.out_opt.out_dir;
}


void Parameters::init()
{
    integration.step = spatial.R / 10;
    spatial.R_eq = 5*spatial.R;

    // HALO REGIME
    switch (generic.reg_halo)
    {
	    case MOD_RHOC_DERIVED:
        {
			   spatial.rho_c = 200 * spatial.rho_0*spatial.c*(1 + spatial.c)*(1 + spatial.c);
			   spatial.M200 = M_SUN_TO_EV(spatial.M200_sun);
			   spatial.R = pow(spatial.M200/(4*M_PI*spatial.rho_c)/(log(1+spatial.c)-spatial.c/(spatial.c+1)) , 1 / 3.0);
			   break;
	    }
        case MOD_M200_DERIVED:
        default:
        {
            break;
        }
	}

    // CHI REGIME
    switch (generic.reg_rho)
    {
        case MOD_CHI_RHOC_DERIVED:
        {
                    for (int i = 0; i < 3; i++){ // iteration to get the right r_eq
                        spatial.rho_c = get_rho_c(spatial.R_eq);
                        spatial.rho_0 = spatial.rho_c*1e-6;
                        spatial.Ms = 4 * M_PI*spatial.rho_c*pow(spatial.R, 3);
                        switch (generic.mod)
                        {
                            case MOD_STAR:
                            {
                                        spatial.R200 = 3 * spatial.R;
                                        break;
                            }
                            case MOD_NFW:
                            {
                                        spatial.R200 = get_R200();
                                        break;
                            }
                        }
                        spatial.c = spatial.R200 / spatial.R;
                    }
                    break;
        }
        case MOD_CHI_YS_DERIVED:
        {
                chi_opt.Ys = get_Ys();
                break;
        }
        case MOD_CHI_REQ_DERIVED:
        default:
        {
                spatial.R_eq = get_r_eq();
                break;
        }
	}


    spatial.Ms =  M_PI*spatial.rho_c*pow(spatial.R, 3);
    chi_opt.chi_0 = 2 * chi_opt.beta*M_PL*chi_opt.Ys;
	chi_opt.chi_B = chi_bulk(spatial.rho_c + spatial.rho_0);
    chi_opt.m_inf = sqrt(V_eff_2nd_derr(chi_opt.chi_0));

    
	integration.r_max = 10 / chi_opt.m_inf;
	integration.h_N = (int)(spatial.R200/integration.step+log(integration.r_max/integration.step));
	integration.h_re = (int)(log(integration.r_max/integration.step));

    // R200
    switch (generic.mod){
	case MOD_STAR:
    {
				spatial.R200 = 10 * spatial.R;
				break;
	}
	case MOD_NFW:
    {
				spatial.R200 = get_R200();
				break;
	}
	}

	spatial.c = spatial.R200 / spatial.R;
	spatial.M200 = spatial.Ms*(log(1+spatial.c)-spatial.c/(1+spatial.c));
	spatial.M200_sun = M_SUN_TO_EV(spatial.M200);

    // ALLOCATION
    for (auto vec : integration.chi)
    {
        vec.reserve(integration.h_N);
    }

}