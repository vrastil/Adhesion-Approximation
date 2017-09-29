
#include "stdafx.h"
#include "core.h"
#include <boost/program_options.hpp>

namespace po = boost::program_options;

using namespace std;

int handle_cmd_line(int ac, char* av[], Sim_Param* sim){
	try {
		string config_file;
		int int_pwr_type;
		// options ONLY on command line
		po::options_description generic("Generic options");
		generic.add_options()
			("help,h", "produce this help message")			
			("config,c", po::value<string>(&config_file)->default_value("INPUT.cfg"), "configuration file name (optional)")
		//	("version,v", "print current version of the program")
			;
			
		// options both on command line	and in configuration file
		po::options_description config_mesh("Mesh configuration");
		config_mesh.add_options()
			("mesh_num,m", po::value<unsigned>(&sim->mesh_num)->default_value(32), "number of mesh cells in the box per dimension (potential)")
            ("mesh_num_pwr,M", po::value<unsigned>(&sim->mesh_num_pwr)->default_value(1024), "number of mesh cells in the box per dimension (power spectrum)")
			("par_num,p", po::value<unsigned>(&sim->Ng)->default_value(2), "ratio of mesh cells and number of particles per dimension")
			("box_size,L", po::value<unsigned>(&sim->box_size)->default_value(512), "box size in units of Mpc/h")
			;
			
		po::options_description config_integ("Integration options");
		config_integ.add_options()
			("redshift,z", po::value<double>(&sim->z_in)->default_value(200.), "redshift at the start of the simulation")
			("redshift_0,Z", po::value<double>(&sim->z_out)->default_value(10.), "redshift at the end of the simulation")
            ("time_step,b", po::value<double>(&sim->db)->default_value(0.1, "0.1"), "dimensionless time-step (scale factor)")
            ("print_every", po::value<unsigned>(&sim->print_every)->default_value(1, "1"), "save particle positions and power spectrum every n-th step")
			;
		
		po::options_description config_app("Approximations");
		config_app.add_options()
			("comp_ZA", po::value<bool>(&sim->comp_ZA)->default_value(false), "compute Zeldovich approximation")
			("comp_FF", po::value<bool>(&sim->comp_FF)->default_value(false), "compute Frozen-flow approximation")
			("comp_FP", po::value<bool>(&sim->comp_FP)->default_value(false), "compute Frozen-potential approximation")
			("comp_AA", po::value<bool>(&sim->comp_AA)->default_value(false), "compute Adhesion approximation")
			("comp_FP_pp", po::value<bool>(&sim->comp_FP_pp)->default_value(false),
								"compute Frozen-potential approximation (particle-particle interaction)")
			;
			
		po::options_description config_power("Power spectrum options");
		config_power.add_options()
			("pwr_type,P", po::value<int>(&int_pwr_type)->default_value(0), "power spectrum type")
			("index,n", po::value<double>(&sim->power.ns)->default_value(1.), "spectral index of the scale-free power spectrum")
			("sigma8,s", po::value<double>(&sim->power.s8)->default_value(1.), "normalization of the power spectrum at R = 8 Mpc/h")
			("smoothing_k,k", po::value<double>(&sim->power.k2_G)->default_value(0.),
								"smoothing wavenumber of TZA in units of h/Mpc, set 0 for ZA")
			;
		
		po::options_description config_run("Run options");
		config_run.add_options()
			("out_dir,o", po::value<string>(&sim->out_dir)->default_value("output/"), "output folder name")
			("num_thread,t", po::value<unsigned>(&sim->nt)->default_value(0), "number of threads the program will use, set 0 for max. available")
			;
		
		po::options_description config_other("Other options");
		config_other.add_options()
			("viscosity,v", po::value<double>(&sim->nu)->default_value(1.), "'viscozity' for adhesion approximation in units of pixel^2")
            ("cut_radius,r", po::value<double>(&sim->rs)->default_value(2.7, "2.7"), "short-range force cutoff radius in units of mesh cells")
			;
			
		po::options_description cmdline_options("\nCOSMOLOGICAL APPROXIMATION");
        cmdline_options.add(generic).add(config_mesh).add(config_power).add(config_integ).add(config_app).add(config_run).add(config_other);
		
		po::variables_map vm;
		store(po::command_line_parser(ac, av).options(cmdline_options).run(), vm);		
		po::notify(vm);

		sim->power.pwr_type = static_cast<e_power_spec>(int_pwr_type);

		if (vm.count("help")) {
            cout << setprecision(3) << cmdline_options << "\n";
            return 2;
        }
		
		if (vm.count("config")) {
            ifstream ifs(config_file.c_str());
			if (ifs){
			cout << "\nUsing configuration options defined in file " << config_file << " (given command line options have higher priority).\n";
				store(po::parse_config_file(ifs, cmdline_options), vm);
				notify(vm);
			} else{
				cout << "can not open config file: " << config_file << "\n";
				return 0;
			}
        }
		
	}
	catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        return 1;
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
    }
	return 0;
}