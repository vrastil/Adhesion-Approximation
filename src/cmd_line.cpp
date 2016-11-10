#include "stdafx.h"
#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <iterator>
using namespace std;

extern int par_num, mesh_num, box_size, nt;
extern double ns, k2_G, s8, z;

extern string out_dir;
string config_file;

int handle_cmd_line(int ac, char* av[]){
	try {
		// options ONLY on command line
		po::options_description generic("Generic options");
		generic.add_options()
			("help,h", "produce this help message")			
			("config,c", po::value<string>(&config_file)->default_value("INPUT.cfg"), "configuration file name (optional")
		//	("version,v", "print current version of the program")
			;
			
		// options both on command line	and in configuration file
		po::options_description config("Configuration");
		config.add_options()
			("mesh_num,m", po::value<int>(&mesh_num)->default_value(32), "number of mesh cells in the box per dimension")
			("par_num,p", po::value<int>(&par_num)->default_value(32), "number of particles in the box per dimension")
			("box_size,b", po::value<int>(&box_size)->default_value(512), "box size in units of Mpc/h")
			("redshift,z", po::value<double>(&z)->default_value(200.), "redshift at the start of the simulation")
			("index,n", po::value<double>(&ns)->default_value(1.), "spectral index of the scale-free power spectrum")
			("sigma8,s", po::value<double>(&s8)->default_value(1.), "normalization of the power spectrum at R = 8 Mpc/h")
			("smoothing_k,k", po::value<double>(&k2_G)->default_value(0.), "smoothing wavenumber of TZA in units of h/Mpc, set 0 for ZA")
			("out_dir,o", po::value<string>(&out_dir)->default_value("output/"), "output folder name")
			("num_thread,t", po::value<int>(&nt)->default_value(2), "number of threads the program will use")
			;
			
		po::options_description cmdline_options("\nADHESION APPROXIMATION");
        cmdline_options.add(generic).add(config);
		
		po::variables_map vm;
		store(po::command_line_parser(ac, av).options(cmdline_options).run(), vm);		
		po::notify(vm);
		

		if (vm.count("help")) {
            cout << cmdline_options << "\n";
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