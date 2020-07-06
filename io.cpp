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

namespace po = boost::program_options;

/****************************//**
 * PUBLIC FUNCTIONS DEFINITIONS *
 ********************************/

 Parameters param;

status_t handle_cmd_line(int ac, const char* const av[]){
    std::string config_file;
    status_t status = status_t::status_OK;

    // options ONLY on command line
    po::options_description generic("Generic options");
    generic.add_options()
        ("help,h", "produce this help message")			
        ("config,c", po::value<std::string>(&config_file)->default_value("input.cfg"), "configuration file name (optional)")
        ;
        
    // options both on command line	and in configuration file
    po::options_description chi_opt("Chameleon parameters options");
    chi_opt.add_options()
        ("mod", po::value<int>(&param.generic.mod)->default_value(0), "Mode; 0 = star, 1 = NFW")
        ("n", po::value<double>(&param.chi_opt.n)->default_value(0.5), "chameleon power-law potential exponent,0 < n < 1")
        ;
    
    po::options_description config_output("Output options");
    config_output.add_options()
        ("out_dir,o", po::value<std::string>(&param.out_opt.out_dir)->default_value("output/"), "output folder name")
        ;


    po::options_description cmdline_options("\n" "PROJECT" " v" "PROJECT_VERSION");//< store all normal parameters
    cmdline_options.add(generic).add(chi_opt).add(config_output);

    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, cmdline_options), vm);
    po::notify(vm);

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
        "Parameters given through command line options or through configuration file:\n\n"
        "Generic options:\n"
        "\tmod = " << p_mod << "\n"

        "Chameleon parameters options:\n"
        "\tn = " << param.chi_opt.n << "\n"

        "Output options:\n"
        "\tout dir" << p_out_dir;
}