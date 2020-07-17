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
#include "units.hpp"
#include "stars_NFW.hpp"

namespace po = boost::program_options;

constexpr size_t MAX_MEMORY_SIZE = 1024UL * 1024UL * 200UL; // 200 MB

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
        ("Omega_m", po::value<double>(&param.spatial.Omega_m)->default_value(4E4), "background density")
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
        "\tmod = " << param.integration.mod <<

        "\nIntegration options:\n"
        "\terr = " << param.integration.err <<
        "\tstep / R_s = " << param.integration.step / get_radius_scale() <<
        "\tr_max / R_s = " << param.integration.r_max / get_radius_scale() <<

        "\nSpatial parameters:\n"
        "\tc = " << param.spatial.c <<
        "\tR / R_s = " << param.spatial.R / get_radius_scale() <<
        "\tR200 / R_s = " << param.spatial.R200 / get_radius_scale() <<
        "\trho_c = " << param.spatial.rho_c <<
        "\trho_0 = " << param.spatial.rho_0 <<
        "\tM200 / M_s = " << param.spatial.M200 / get_mass_scale() <<
        "\tMs / M_s = " << param.spatial.Ms / get_mass_scale() <<
        "\tM200_sun = " << param.spatial.M200_sun <<

        "\nChameleon parameters:\n"
        "\tn = " << param.chi_opt.n <<
        "\tbeta = " << param.chi_opt.beta <<
        "\tYs = " << param.chi_opt.Ys <<
        "\tR_eq = " << param.spatial.R_eq <<
        "\tchi_B = " << param.chi_opt.chi_B <<
        "\tm_inf = " << param.chi_opt.m_inf <<

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

    // ALLOCATION
    size_t capacity = std::min((size_t)(spatial.R200/integration.step+log(integration.r_max/integration.step)), MAX_MEMORY_SIZE);
    BOOST_LOG_TRIVIAL(info) << "Capacity: " << capacity;
    integration.chi.reserve(capacity);

    // DIRECTORY STRUCT
    out_opt.out_dir = std_out_dir();
    create_dir(out_opt.out_dir);
}

namespace fs = boost::filesystem;

#define BUFF_SIZE 1024 * 1024 * 16 // 16 MB buffer

Ofstream::Ofstream(std::string file_name) : std::ofstream(file_name), buf(new char[BUFF_SIZE])
{
    if (!this->is_open()) throw std::runtime_error("Error while opening '" + file_name + "'");
    this->rdbuf()->pubsetbuf(buf, sizeof(buf));
}

Ofstream::~Ofstream()
{
    delete[] buf;
    if (this->is_open()) this->close();
}

Ifstream::Ifstream(std::string file_name) : std::ifstream(file_name)
{
    if (!this->is_open()) throw std::runtime_error("Error while opening '" + file_name + "'");
}

Ifstream::~Ifstream()
{
    if (this->is_open()) this->close();
}

// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
std::string currentDateTime()
{
	const time_t     now = time(0);
	struct tm  tstruct;
	char       buf[80];
	gmtime_r(&now, &tstruct);
	strftime(buf, sizeof(buf), "%y%m%d_%H%M%S", &tstruct);
	
	std::string returnval(buf);
    return returnval;
}

std::string std_out_dir()
{
    /// directory name = YYMMDD_HHMMSS_m_p_M_b
    std::string out_dir = param.out_opt.out_dir;
    // switch (param.integration.mod)
    // {
    //     case MOD_STAR: out_dir += "STAR/"; break;
    //     case MOD_NFW: out_dir += "NFW/"; break;
    // }
    // out_dir += currentDateTime();

    // const std::string out_dir = param.out_opt.out_dir + std::to_string(param.integration.mod) +"m_" +
    //                   std::to_string(param.spatial.R) + "R_" + std::to_string(param.spatial.c) +"c_" + 
    //                   std::to_string(param.chi_opt.Ys) + "Y";

    // TODO: for now do not check existing directories, overwrite
    return out_dir;

    /// check if directory exists
    if (!fs::exists(fs::path(out_dir.c_str()))) return out_dir + "/";

    /// try appending numbers starting at 2
    else
    {
        for(size_t i = 2; ; ++i)
        {
            const std::string out_dir_new = out_dir  + "_" + std::to_string(i);
            if (!fs::exists(fs::path(out_dir_new.c_str()))) return out_dir_new + "/";
        }
    }
}

void create_dir(const std::string &out_dir)
{
	const fs::path dir(out_dir.c_str());
	if(fs::create_directories(dir))
    {
        BOOST_LOG_TRIVIAL(debug) << "Directory created: "<< out_dir;
    }
}

void remove_dir(const std::string &out_dir)
{
    const fs::path dir(out_dir.c_str());
    if (fs::remove_all(dir))
    {
        BOOST_LOG_TRIVIAL(debug) << "Directory removed: "<< out_dir;
    }
}

void remove_all_files(const std::string &out_dir)
{
    const fs::path dir(out_dir.c_str());
    size_t i = 0;

    for(auto & p : fs::directory_iterator(dir))
    {
        if (fs::is_regular_file(p))
        {
            fs::remove(p);
            ++i;
        }
    }
    BOOST_LOG_TRIVIAL(debug) << "Removed " << i << " file(s) in directory: "<< out_dir;
}