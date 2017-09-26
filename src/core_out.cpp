
#include "stdafx.h"
#include <boost/filesystem.hpp>
#include "core.h"

namespace fs = boost::filesystem;
using namespace std;

// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
string currentDateTime()
{
	time_t     now = time(0);
	struct tm  tstruct;
	char       buf[80];
	tstruct = *gmtime(&now);
	strftime(buf, sizeof(buf), "%y%m%d_%H%M%S", &tstruct);
	
	string returnval(buf);
    return returnval;
}

string std_out_dir(string pre_subdir, const Sim_Param &sim)
{
    return sim.out_dir + pre_subdir + currentDateTime() + "_" + to_string(sim.mesh_num) +"m_" +
           to_string(sim.Ng) + "p_" + to_string(sim.mesh_num_pwr) +"M_" + to_string(sim.box_size) + "b/";
}

void create_dir(string out_dir)
{
	fs::path dir(out_dir.c_str());
	if(fs::create_directories(dir)){
        cout << "Directory Created: "<< out_dir << endl;
    }
}
void work_dir_over(string out_dir)
{
	create_dir(out_dir + "par_cut/");
	create_dir(out_dir + "pwr_diff/");
	create_dir(out_dir + "pwr_spec/");
	create_dir(out_dir + "rho_map/");
	create_dir(out_dir + "rho_bin/");
}

void print_pow_spec(const vector<double_2> &pwr_spec_binned, string out_dir, string suffix)
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

void print_pow_spec_diff(const vector<double_2> &pwr_spec_binned, const vector<double_2> &pwr_spec_binned_0,
	double b, string out_dir, string suffix)
{
	out_dir += "pwr_diff/";
	FILE* pFile;
	pFile = fopen((out_dir + "pwr_spec_diff" + suffix + ".dat").c_str(), "w");
	cout << "Writing power spectrum difference into file " << out_dir + "pwr_spec_diff" + suffix + ".dat\n";
	fprintf (pFile, "# This file contains relative difference between power spectrum P(k) and lineary extrapolated power spectrum depending on wavenumber k in units [h/Mpc].\n");
	fprintf (pFile, "# k [h/Mpc]\t(P(k, z)-P_lin(k, z))/P_lin(k, z)\n");
	
	double P_k, P_ZA;
	
	for (unsigned j = 0; j < pwr_spec_binned.size(); j++){
		if (pwr_spec_binned[j][0] == pwr_spec_binned_0[j][0]){
			P_k = pwr_spec_binned[j][1];
			P_ZA = pwr_spec_binned_0[j][1] * pow(b, 2.);
			if((P_ZA) && (P_k)) fprintf (pFile, "%f\t%f\n", pwr_spec_binned[j][0], (P_k-P_ZA)/P_ZA);
		} else printf ("WARNING! Binned power spectra don`t match each other! k = %f, while k_lin = %f\n", pwr_spec_binned[j][0], pwr_spec_binned_0[j][0]);
	}

	fclose (pFile);
}

void print_track_par(const Tracking& track, const Sim_Param &sim, string out_dir, string suffix)
{
	out_dir += "par_cut/";
	FILE* pFile = fopen((out_dir + "track_par_pos" + suffix + ".dat").c_str(), "w");
    double x,y,z;
    double x_0 = sim.x_0();
	cout << "Writing positons of " << track.num_track_par << " tracked particles into file " << out_dir + "track_par_pos" + suffix + ".dat\n";
	fprintf (pFile, "# This file contains positions of particles in units [Mpc/h].\n");
	fprintf (pFile, "# x [Mpc/h]\tz [Mpc/h]\n");
	for (int i=0; i<track.num_track_par; i++){
		for (int j=0; j<track.num_step();j++){
			x = track.par_pos[j][i].position[0];
			y = track.par_pos[j][i].position[1];
			z = track.par_pos[j][i].position[2];
			fprintf (pFile, "%f\t%f\t%f\n", x*x_0 , z*x_0, y*x_0);
		}
		fprintf (pFile, "\n\n");
	}
	fclose (pFile);
}

void print_rho_map(const Mesh& delta, const Sim_Param &sim, string out_dir, string suffix)
{
    out_dir += "rho_map/";
    double x_0 = sim.x_0_pwr();
	FILE* pFile;
	pFile = fopen((out_dir + "rho_map" + suffix + ".dat").c_str(), "w");
	cout << "Writing density map into file " << out_dir + "rho_map" + suffix + ".dat\n";
	fprintf (pFile, "# This file contains density map delta(x).\n");
	fprintf (pFile, "# x [Mpc/h]\tz [Mpc/h]\tdelta\n");
	for (int i = 0; i < sim.mesh_num_pwr; i++){
		for (int j = 0; j < sim.mesh_num_pwr; j++){
			fprintf (pFile, "%f\t%f\t%f\n", i*x_0, j*x_0, delta(i, sim.mesh_num_pwr/2, j));
		}
		fprintf (pFile, "\n");
	}

	fclose (pFile);
}

void print_projected_rho(const Mesh& delta, const Sim_Param &sim, string out_dir, string suffix)
{
    out_dir += "rho_map/";
    double x_0 = sim.x_0_pwr();
	FILE* pFile;
	pFile = fopen((out_dir + "rho_map_projected" + suffix + ".dat").c_str(), "w");
	cout << "Writing density map into file " << out_dir + "rho_map" + suffix + ".dat\n";
	fprintf (pFile, "# This file contains density map delta(x).\n");
	fprintf (pFile, "# x [Mpc/h]\tz [Mpc/h]\tdelta\n");
	double rho, rho_tmp;
	for (int i = 0; i < sim.mesh_num_pwr; i++){
		for (int j = 0; j < sim.mesh_num_pwr; j++){
			rho = 0;
			for (int k = 0; k < sim.mesh_num_pwr; k++){
				rho_tmp = delta(i, k, j);
			//	if (rho_tmp != -1) printf("Density in (%i, %i, %i) = %f\n", i, j, k, rho_tmp);
				rho+=rho_tmp + 1;
			}
			fprintf (pFile, "%f\t%f\t%f\n", i*x_0, j*x_0, rho -1);
		}
		fprintf (pFile, "\n");
	}

	fclose (pFile);
}

void print_dens_bin(const vector<int> &dens_binned, int mesh_num, string out_dir, string suffix){
	out_dir += "rho_bin/";
	FILE* pFile;
	pFile = fopen((out_dir + "rho_bin" + suffix + ".dat").c_str(), "w");
	cout << "Writing binned density into file " << out_dir + "rho_bin" + suffix + ".dat\n";
	fprintf (pFile, "# This file contains binned density field.\n");
	fprintf (pFile, "# dens\tbin_num\n");
	
	double dens;
	for (unsigned j = 0; j < dens_binned.size(); j++)
	{
		dens = j*0.2-0.9;
		fprintf (pFile, "%f\t%i\n", dens, dens_binned[j]);
        //fprintf (pFile, "%f\t%f\n", dens, dens_binned[j] / pow(mesh_num, 3));        
	}

	fclose (pFile);
}