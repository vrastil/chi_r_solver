/**
 * @brief command line arguments manipulation
 * 
 * @file core_cmd.cpp
 * @author Michal Vrastil
 * @date 2018-07-11
 */

#include <boost/any.hpp>
#include <boost/program_options.hpp>
#include <iomanip>
#include <fstream>

#include "io.hpp"

namespace po = boost::program_options;

/****************************//**
 * PUBLIC FUNCTIONS DEFINITIONS *
 ********************************/

void handle_cmd_line(int ac, const char* const av[], Parameters& param){
    std::string config_file;

    // options ONLY on command line
    po::options_description generic("Generic options");
    generic.add_options()
        ("help,h", "produce this help message")			
        ("config,c", po::value<std::string>(&config_file)->default_value("INPUT.cfg"), "configuration file name (optional)")
        ;
        
    // options both on command line	and in configuration file
    po::options_description chi_opt("Chameleon parameters options");
    config_mesh.add_options()
        ("mesh_num,m", po::value<size_t>(&sim.box_opt.mesh_num)->default_value(128), "number of mesh cells per dimension (potential)")
        ("mesh_num_pwr,M", po::value<size_t>(&sim.box_opt.mesh_num_pwr)->default_value(256), "number of mesh cells per dimension (power spectrum)")
        ("par_num,p", po::value<size_t>(&sim.box_opt.par_num_1d)->default_value(128), "number of particles per dimension")
        ("box_size,L", po::value<FTYPE_t>(&sim.box_opt.box_size)->default_value(512, "512"), "box size in units of Mpc/h")
        ;
    
    po::options_description config_output("Output options");
    config_output.add_options()
        ("print_every", po::value<size_t>(&sim.out_opt.print_every)->default_value(1, "1"), "save particle positions and power spectrum "
                                                                                        "every n-th step, set 0 for no printing")
        ("pwr_bins", po::value<size_t>(&sim.out_opt.bins_per_decade)->default_value(30), "number of bins per decade in power spectrum")
        ("corr_pt", po::value<size_t>(&sim.out_opt.points_per_10_Mpc)->default_value(10), "number of points per 10 Mpc in correlation function")
        ("out_dir,o", po::value<std::string>(&sim.out_opt.out_dir)->default_value("output/"), "output folder name")
        ("print_par_pos", po::value<bool>(&sim.out_opt.print_par_pos)->default_value(false), "print particles positions")
        ("print_dens", po::value<bool>(&sim.out_opt.print_dens)->default_value(false), "print density map and histogram")
        ("print_pwr", po::value<bool>(&sim.out_opt.print_pwr)->default_value(false), "print power spectrum")
        ("print_extrap_pwr", po::value<bool>(&sim.out_opt.print_extrap_pwr)->default_value(false), "print extrapolated power spectrum")
        ("print_corr", po::value<bool>(&sim.out_opt.print_corr)->default_value(false), "print correlation function")
        ("print_vel_pwr", po::value<bool>(&sim.out_opt.print_vel_pwr)->default_value(false), "print velocity power spectrum")
        ;


    po::options_description cmdline_options("\n" "PROJECT" " v" "PROJECT_VERSION");//< store all normal parameters
    cmdline_options.add(generic).add(chi_opt).add(config_output);

    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, cmdline_options), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << std::setprecision(3) << cmdline_options;
        no_run = true;
    }
    else if (vm.count("config")) {
        std::ifstream ifs(config_file.c_str());
        if (ifs){
        BOOST_LOG_TRIVIAL(debug) << "Using configuration options defined in file " << config_file << " (given command line options have higher priority).";
            po::store(po::parse_config_file(ifs, config_test), vm);
            // po::store(po::parse_config_file(ifs, cmdline_options), vm);
            notify(vm);
        } else{
            BOOST_LOG_TRIVIAL(warning) << "Cannot open config file '" << config_file << "'. Using default and command lines values.";
        }
    }
}