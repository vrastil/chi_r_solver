// main.cpp : Defines the entry point for the console application.
//

#include "integrator.hpp"
#include "io.hpp"
#include "stars_NFW.hpp"

int main(int argc, char* argv[])
{
	status_t status = handle_cmd_line(argc, argv);

	if (status == status_t::status_OK)
	{
		param.init();
		if (param.out_opt.print_par) param.print_info();
		solve();
	}

	return int(status);
}

