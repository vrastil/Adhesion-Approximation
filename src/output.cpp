#include "stdafx.h"
#include <boost/filesystem.hpp>
#include <fftw3.h>

namespace fs = boost::filesystem;
using namespace std;
typedef double(*t_power)(double, double*);
const double PI = acos(-1.);

/*string work_dir(string out_dir){
	fs::path dir(out_dir.c_str());
	if(fs::create_directory(dir)){
        cout << "Directory Created: "<< out_dir << endl;
    }
	int i = 0;
	string try_dir;
	
	do{
		i++;
		try_dir = out_dir + "run_" + to_string(i) + "/";
	} while(fs::is_directory(try_dir));
	
	dir = try_dir.c_str();
	fs::create_directory(dir);
	cout << "Directory Created: "<< try_dir << endl;
	return try_dir;
}*/

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

void print_par_pos(int par_num, int mesh_num, int L, double** displ_vec, string out_dir){
	out_dir += "par_cut/";
	FILE* ofile = fopen((out_dir + "par_pos.dat").c_str(), "w");
	cout << "Writing particle positons into file " << out_dir + "par_pos.dat\n";
	fprintf (ofile, "# This file contains positions of particles in units [Mpc/h].\n");
	fprintf (ofile, "# There is in total %i particles in a box of size %i Mpc/h.\n", (int)pow(par_num, 3), L);
	fprintf (ofile, "# x [Mpc/h]\ty [Mpc/h]\tz [Mpc/h]\n");
	for(int i=0; i< pow(par_num, 2); i++){
		for (int j=0; j < par_num; j++){
			fprintf (ofile, "%f\t%f\t%f\n",
				displ_vec[0][i*(mesh_num + 2)+j]/mesh_num*L,
				displ_vec[1][i*(mesh_num + 2)+j]/mesh_num*L,
				displ_vec[2][i*(mesh_num + 2)+j]/mesh_num*L);
		}
	}
	fclose (ofile);
}

void print_par_pos_cut(int par_num, int mesh_num, int L, double** displ_vec, string out_dir, string suffix){
	out_dir += "par_cut/";
	FILE* ofile = fopen((out_dir + "par_cut" + suffix + ".dat").c_str(), "w");
	cout << "Writing cut through the box of particles into file " << out_dir + "par_cut" + suffix + ".dat\n";
	fprintf (ofile, "# This file contains positions of particles in units [Mpc/h].\n");
	double x, y, z, dx;
	for(int i=0; i< pow(par_num, 2); i++){
		for (int j=0; j < par_num; j++){
	
			x = displ_vec[0][i*(par_num + 2)+j];
			y = displ_vec[1][i*(par_num + 2)+j];
			z = displ_vec[2][i*(par_num + 2)+j];
			dx = abs(y - mesh_num/2.);
			if (dx < 0.5){
				// cut (L/4 x L/4 x 0.5)
				fprintf (ofile, "%f\t%f\n", x/mesh_num*L , z/mesh_num*L);
			}
		}
	}
	fclose (ofile);
}

void print_par_pos_cut_small(int par_num, int mesh_num, int L, double** displ_vec, string out_dir, string suffix){
	out_dir += "par_cut/";
	FILE* ofile = fopen((out_dir + "par_cut" + suffix + ".dat").c_str(), "w");
	cout << "Writing small cut through the box of particles into file " << out_dir + "par_cut" + suffix + ".dat\n";
	fprintf (ofile, "# This file contains positions of particles in units [Mpc/h].\n");
	double x, y, z, dx;
	for(int i=0; i< pow(par_num, 2); i++){
		for (int j=0; j < par_num; j++){
	
			x = displ_vec[0][i*(par_num + 2)+j];
			y = displ_vec[1][i*(par_num + 2)+j];
			z = displ_vec[2][i*(par_num + 2)+j];			
			dx = abs(y - mesh_num/2.);
			if ((dx < 0.5) && (x < mesh_num/4.) && (z < mesh_num/4.)){
				// cut (L/4 x L/4 x 0.5)
				fprintf (ofile, "%f\t%f\n", x/mesh_num*L , z/mesh_num*L);
			}
		}
	}
	fclose (ofile);
}

void gen_pow_spec_binned(int mesh_num, fftw_complex* pwr_spec, fftw_complex* pwr_spec_binned, double k_min, double k_max, int bin_num){

	double log_bin = pow(k_max / k_min, 1./bin_num);
	double k;
	int bin;
	printf("Computing binned power spectrum P(k)...\n");
	for (int j = 0; j < bin_num; j++){
		pwr_spec_binned[j][0] = 0;
		pwr_spec_binned[j][1] = 0;
	}	
	for (int i = 0; i < pow(mesh_num, 2)*(mesh_num/2 + 1); i++){
		k = pwr_spec[i][0];
		if ((k <=k_max) && (k>=k_min)){
			bin = (int)(log(k/k_min)/log(log_bin));
			pwr_spec_binned[bin][1] += pwr_spec[i][1]; // P(k)
			pwr_spec_binned[bin][0]++;
		}
	}
		
	k = k_min*sqrt(log_bin);
	for (int j = 0; j < bin_num; j++){
		if (pwr_spec_binned[j][0]) pwr_spec_binned[j][1] /= pwr_spec_binned[j][0];
		pwr_spec_binned[j][0] = k;
		k *= log_bin;
	}
}

void print_pow_spec(int mesh_num, fftw_complex* pwr_spec_binned, string out_dir, string suffix, double k_min, double k_max, int bin_num){
	out_dir += "pwr_spec/";	
	FILE* pFile;
	pFile = fopen((out_dir + "pwr_spec" + suffix + ".dat").c_str(), "w");
	cout << "Writing power spectrum into file " << out_dir + "pwr_spec" + suffix + ".dat\n";
	fprintf (pFile, "# This file contains power spectrum P(k) in units [(Mpc/h)^3] depending on wavenumber k in units [h/Mpc].\n");
	fprintf (pFile, "# k [h/Mpc]\tP(k) [(Mpc/h)^3]\n");
	
	for (int j = 0; j < bin_num; j++){
		if (pwr_spec_binned[j][1]) fprintf (pFile, "%f\t%f\n",  pwr_spec_binned[j][0], pwr_spec_binned[j][1]);
	}

	fclose (pFile);
}

void print_pow_spec_diff(fftw_complex* pwr_spec_binned, fftw_complex* pwr_spec_binned_0, string out_dir, string suffix, int bin_num, double b){
	out_dir += "pwr_diff/";
	FILE* pFile;
	pFile = fopen((out_dir + "pwr_spec" + suffix + ".dat").c_str(), "w");
	cout << "Writing power spectrum into file " << out_dir + "pwr_spec" + suffix + ".dat\n";
	fprintf (pFile, "# This file contains relative difference between power spectrum P(k) and lineary extrapolated power spectrum depending on wavenumber k in units [h/Mpc].\n");
	fprintf (pFile, "# k [h/Mpc]\t(P(k, z)-P_ZA(k, z))/P_ZA(k, z)\n");
	
	double P_k, P_ZA;
	
	for (int j = 0; j < bin_num; j++){
		if (pwr_spec_binned[j][0] == pwr_spec_binned_0[j][0]){
			P_k = pwr_spec_binned[j][1];
			P_ZA = pwr_spec_binned_0[j][1] * pow(b, 2.);
			if((P_ZA) && (P_k)) fprintf (pFile, "%f\t%f\n", pwr_spec_binned[j][0], (P_k-P_ZA)/P_ZA);
		} else printf ("WARNING! Binned power spectra don`t match each other! k = %f, while k_ZA = %f\n", pwr_spec_binned[j][0], pwr_spec_binned_0[j][0]);
	}

	fclose (pFile);
}

typedef double t_supp[2];

void upd_suppresion(t_supp* supp, fftw_complex* pwr_spec_binned, fftw_complex* pwr_spec_binned_0, double b, int step)
{
	double P_k, P_ZA, supp_tmp = 0;
	int i = 0;
	for (int j = 0; j < 100; j++){
		P_k = pwr_spec_binned[j][1];
		P_ZA = pwr_spec_binned_0[j][1] * pow(b, 2.);
		if((P_ZA) && (P_k))
		{
			supp_tmp += (P_k-P_ZA)/P_ZA;
			i++;
		}
		if (i == 10) break;
	}
	supp[step][0] = b;
	supp[step][1] = supp_tmp / i;
}

void print_suppression(t_supp* supp, string out_dir, int num_step, string suffix){
	FILE* pFile;
	out_dir += "supp/";
	pFile = fopen((out_dir + "supp" + suffix + ".dat").c_str(), "w");
	cout << "Writing power spectrum suppresion into file " << out_dir + "suppresion" + ".dat\n";
	fprintf (pFile, "# This file contains power spectrum suppresion, i.e. relative difference between power spectrum P(k) and lineary extrapolated power spectrum depending on time.\n");
	
	for (int j = 0; j < num_step; j++){
		fprintf (pFile, "%f\t%f\t%f\n", supp[j][0], supp[j][1], 1/supp[j][0]-1);
	}

	fclose (pFile);
}

void get_track_par_id(int par_num, int* id, int track_num){
	int x, y, z, k;
	double s;
	y = par_num / 2; // middle of the cube
	s = par_num / (4.*(track_num+1.)); // quarter of the cube
	k = 0;
	for (int i=1; i<=track_num;i++){
		z = (int)(s*i);
		for (int j=1; j<=track_num;j++){
			x = (int)(s*j);
			id[k] = z*par_num*(par_num + 2)+y*(par_num + 2)+x;
			k++;
		}
	}
}

void update_track_par(double*** track_pos, double** displ_vec, int step, int* id, int track_num){
	for (int i=0; i<track_num*track_num; i++){			
		track_pos[i][step][0] = displ_vec[0][id[i]]; // z
		track_pos[i][step][1] = displ_vec[2][id[i]]; // x
	}
}

void print_track_par(double*** track_pos, int step, int track_num, int mesh_num, int L, string out_dir, string suffix){
	out_dir += "par_cut/";
	FILE* ofile = fopen((out_dir + "track_par_pos" + suffix + ".dat").c_str(), "w");
		
	cout << "Writing positons of " << track_num*track_num << " tracked particles into file " << out_dir + "track_par_pos" + suffix + ".dat\n";
	fprintf (ofile, "# This file contains positions of particles in units [Mpc/h].\n");
	fprintf (ofile, "# x [Mpc/h]\tz [Mpc/h]\n");
	for (int i=0; i<track_num*track_num; i++){
		for (int j=0; j<step;j++){
			fprintf (ofile, "%f\t%f\n", track_pos[i][j][0]/mesh_num*L, track_pos[i][j][1]/mesh_num*L);
		}
		fprintf (ofile, "\n\n");
	}
	fclose (ofile);
}

void print_rho_map(int mesh_num, int L,  double* delta, string out_dir, string suffix){
	out_dir += "rho_map/";
	FILE* pFile;
	pFile = fopen((out_dir + "rho_map" + suffix + ".dat").c_str(), "w");
	cout << "Writing density map into file " << out_dir + "rho_map" + suffix + ".dat\n";
	fprintf (pFile, "# This file contains density map delta(x).\n");
	fprintf (pFile, "# x [Mpc/h]\tz [Mpc/h]\tdelta\n");
	int pos;
	for (int i = 0; i < mesh_num; i++){
		pos = i*mesh_num*(mesh_num+2)+mesh_num/2*(mesh_num+2); // (z, y, ?) position
		for (int j = 0; j < mesh_num; j++){
			fprintf (pFile, "%i\t%i\t%f\n", j*L/mesh_num, i*L/mesh_num, delta[pos+j]);
		}
		fprintf (pFile, "\n");
	}

	fclose (pFile);
}


void gen_dens_binned(int mesh_num, int Ng, double* rho, int* dens_binned, int bin_num)
{
	int bin;
	double rho_avg;
	printf("Computing binned density field...\n");
	for (int j = 0; j < bin_num; j++){
		dens_binned[j] = 0;
	}
	
	for (int i = 0; i < mesh_num; i+=Ng)
	{
		for (int j = 0; j < mesh_num; j+=Ng)
		{
			for (int k = 0; k < mesh_num; k+=Ng)
			{
				// Need to go through all mesh cells [i, i+Ng-1]*[j, j+Ng-1], [k, k+Ng, -1]
				rho_avg = 0;
				for (int ii = i; ii < i+Ng; ii++)
				{
					for (int jj = j; jj  < j+Ng; jj++)
					{
						for (int kk = k; kk < k+Ng; kk++)
						{
							rho_avg+=rho[ii*mesh_num*(mesh_num+2) + jj*(mesh_num+2) + kk];
						}
					}
				}
				rho_avg /= pow(Ng, 3);
				bin = (int)((rho_avg+1)/0.2);
				if (bin >= bin_num) bin = bin_num - 1;
				dens_binned[bin]++;
			}
		}
	}
}

void print_dens_bin(int* dens_binned, string out_dir, string suffix, int bin_num, int mesh_num){
	out_dir += "rho_bin/";
	FILE* pFile;
	pFile = fopen((out_dir + "rho_bin" + suffix + ".dat").c_str(), "w");
	cout << "Writing power spectrum into file " << out_dir + "rho_bin" + suffix + ".dat\n";
	fprintf (pFile, "# This file contains binned density field.\n");
	fprintf (pFile, "# dens\tbin_num\n");
	
	double dens;
	for (int j = 0; j < bin_num; j++)
	{
		dens = j*0.2-1;
		fprintf (pFile, "%f\t%f\n", dens, dens_binned[j] / pow(mesh_num, 3));

	}

	fclose (pFile);
}