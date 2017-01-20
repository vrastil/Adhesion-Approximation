
#include "stdafx.h"
#include <boost/filesystem.hpp>
#include "core.h"

namespace fs = boost::filesystem;
using namespace std;

void work_dir_over(string out_dir){
	string wrk_dir;
	fs::path dir(out_dir.c_str());
	if(fs::create_directory(dir)){
        cout << "Directory Created: "<< out_dir << endl;
    }
	
	wrk_dir = out_dir + "par_cut/";
	dir = wrk_dir.c_str();
	if(fs::create_directory(dir)){
        cout << "Directory Created: "<< wrk_dir << endl;
    }
	
	wrk_dir = out_dir + "pwr_diff/";
	dir = wrk_dir.c_str();
	if(fs::create_directory(dir)){
        cout << "Directory Created: "<< wrk_dir << endl;
    }
	
	wrk_dir = out_dir + "pwr_spec/";
	dir = wrk_dir.c_str();
	if(fs::create_directory(dir)){
        cout << "Directory Created: "<< wrk_dir << endl;
    }
	
	wrk_dir = out_dir + "rho_map/";
	dir = wrk_dir.c_str();
	if(fs::create_directory(dir)){
        cout << "Directory Created: "<< wrk_dir << endl;
    }
	
	wrk_dir = out_dir + "supp/";
	dir = wrk_dir.c_str();
	if(fs::create_directory(dir)){
        cout << "Directory Created: "<< wrk_dir << endl;
    }
	
	wrk_dir = out_dir + "rho_bin/";
	dir = wrk_dir.c_str();
	if(fs::create_directory(dir)){
        cout << "Directory Created: "<< wrk_dir << endl;
    }
}

void print_pow_spec(const vector<fftw_complex> &pwr_spec_binned, string out_dir, string suffix)
{
	out_dir += "pwr_spec/";
	string file_name = out_dir + "pwr_spec" + suffix + ".dat";
	FILE* pFile;
	pFile = fopen(file_name.c_str(), "w");
	if (pFile == NULL)
	{
		printf("Error while opening %s\n", file_name.c_str());
		perror("Error");
		return;
	}
	
	cout << "Writing power spectrum into file " << file_name << endl;
	fprintf (pFile, "# This file contains power spectrum P(k) in units [(Mpc/h)^3] depending on wavenumber k in units [h/Mpc].\n");
	fprintf (pFile, "# k [h/Mpc]\tP(k) [(Mpc/h)^3]\n");
	
	for (unsigned j = 0; j < pwr_spec_binned.size(); j++){
		if (pwr_spec_binned[j][1]) fprintf (pFile, "%f\t%f\n",  pwr_spec_binned[j][0], pwr_spec_binned[j][1]);
	}

	fclose (pFile);
}

void print_par_pos_cut_small(Particle_x* particles, const Sim_Param &sim, string out_dir, string suffix)
{
	out_dir += "par_cut/";
	string file_name = out_dir + "par_cut" + suffix + ".dat";
	FILE* pFile;
	pFile = fopen(file_name.c_str(), "w");
	if (pFile == NULL)
	{
		printf("Error while opening %s\n", file_name.c_str());
		perror("Error");
		return;
	}
	
	cout << "Writing small cut through the box of particles into file " << out_dir + "par_cut" + suffix + ".dat\n";
	fprintf (pFile, "# This file contains positions of particles in units [Mpc/h].\n");
	double x, y, z, dx;
	for(int i=0; i < sim.par_num; i++)
	{
		x = particles[i].position.x;
		y = particles[i].position.y;
		z = particles[i].position.z;			
		dx = abs(y - sim.mesh_num/2.);
		if ((dx < 0.5) && (x < sim.mesh_num/4.) && (z < sim.mesh_num/4.))
		{
			// cut (L/4 x L/4 x 0.5)
			fprintf (pFile, "%f\t%f\t%f\n", x/sim.mesh_num*sim.box_size , z/sim.mesh_num*sim.box_size, y/sim.mesh_num*sim.box_size);
		}
	}
	fclose (pFile);
}

void print_track_par(const Tracking& track, const Sim_Param &sim, string out_dir, string suffix){
	out_dir += "par_cut/";
	FILE* ofile = fopen((out_dir + "track_par_pos" + suffix + ".dat").c_str(), "w");
		
	cout << "Writing positons of " << track_num*track_num << " tracked particles into file " << out_dir + "track_par_pos" + suffix + ".dat\n";
	fprintf (ofile, "# This file contains positions of particles in units [Mpc/h].\n");
	fprintf (ofile, "# x [Mpc/h]\tz [Mpc/h]\n");
	for (int i=0; i<track_num*track_num; i++){
		for (int j=0; j<step;j++){
			fprintf (ofile, "%f\t%f\n", track_pos[i][j][0]/mesh_num*L, track_pos[i][j][1]/mesh_num*L);
			fprintf (pFile, "%f\t%f\t%f\n", x/sim.mesh_num*sim.box_size , z/sim.mesh_num*sim.box_size, y/sim.mesh_num*sim.box_size);
		}
		fprintf (ofile, "\n\n");
	}
	fclose (ofile);
}