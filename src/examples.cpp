
#include "stdafx.h"
#include "core.h"
#include "core_out.h"
#include "core_app.h"
#include "core_mesh.h"


using namespace std;
const double PI = acos(-1.);

void print_par_pos_cut(Particle_v* particles, const Sim_Param &sim, string out_dir, string suffix)
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
	
	cout << "Writing particles into file " << out_dir + "par_cut" + suffix + ".dat\n";
	fprintf (pFile, "# This file contains positions of particles in units [Mpc/h].\n");
	double x, y, z, dx;
	for(int i=0; i < sim.par_num; i++)
	{
		x = particles[i].position.x;
		y = particles[i].position.y;
		z = particles[i].position.z;
		fprintf (pFile, "%f\t%f\t%f\n", x*sim.x_0() , z*sim.x_0(), y*sim.x_0());
	}
	fclose (pFile);
}

void print_force(Particle_v* particles, const std::vector<Mesh>& app_field, const Sim_Param &sim, string out_dir, string suffix,
				const LinkedList& linked_list)
{
	out_dir += "force/";
	string file_name = out_dir + "force" + suffix + ".dat";
	FILE* pFile;
	pFile = fopen(file_name.c_str(), "w");
	if (pFile == NULL)
	{
		printf("Error while opening %s\n", file_name.c_str());
		perror("Error");
		return;
	}
	
	cout << "Writing force into file " << out_dir + "par_cut" + suffix + ".dat\n";
	fprintf (pFile, "# This file contains positions of particles in units [Mpc/h].\n");
	
	Vec_3D<double> f_tmp, f_total;
	Vec_3D<double> position (sim.mesh_num/2., sim.mesh_num/2., sim.mesh_num/2.);
	Vec_3D<double> position_0(sim.mesh_num/2., sim.mesh_num/2., sim.mesh_num/2.);
	double dr = 0.09;
	double r, rper;
	double fl, fs = 0, ft;
	int i_max = 90;
	position.x+=dr;
	fprintf(pFile, "# r\tfs\tfl\tfs+fl\tft\n");
	for (int i = 1; i < i_max; i++)
	{
		position.x=position_0.x+pow(1+dr, i) - 1 + 0.5;
		r = position.x -position_0.x;
		rper = get_distance(position_0, position, sim.mesh_num);
	//	printf("\tCalculating force at distance r = %f / %f\n", r, rper);
	//	z = (Vec_3D<int>)(position / sim.Hc); // chain position of particle
		f_tmp.assign(0., 0., 0.);
		f_total.assign(0., 0., 0.);
		
		assign_from(app_field, position, &f_tmp, sim.order);
		fl = f_tmp.norm();
		printf("F_long = (%f, %f, %f)\t", f_tmp.x, f_tmp.y, f_tmp.z);
		f_total+=f_tmp;
		f_tmp.assign(0., 0., 0.);
		
		force_short(sim, linked_list, particles, position, &f_tmp);
		printf("F_short = (%f, %f, %f)\n", f_tmp.x, f_tmp.y, f_tmp.z);
		fs = f_tmp.norm();
		f_total+=f_tmp;
		ft = f_total.norm();
		if (r > sim.mesh_num/2.) { printf("distance greater than N/2, dst = %f\n", r); break;}
//		printf("# Position = (%f, %f, %f)\n", position.x, position.y, position.z);
//		printf("# HOC (position) %i\n", APP.linked_list.HOC(z));
		fprintf(pFile, "%f\t%f\t%f\t%f\t%f\t%f\n",rper , fs, fl, fs+fl, ft, r);
	}
	fclose (pFile);
}

int example(Sim_Param &sim)
{
//	sim.par_num = 1;
//	sim.Ng = sim.mesh_num;
//	sim.order = 2;
	
	sim.print_info(); // print simulation parameters

	cout << "\n"
	"*******\n"
	"EXAMPLE\n"
	"*******\n";
	
	string out_dir_app = sim.out_dir + "example/";
	work_dir_over(out_dir_app);
	create_dir(out_dir_app + "force/");
	
	/** ALLOCATION OF MEMORY + FFTW PREPARATION **/
	App_Var_FP_mod APP(sim, "_exmp_");
	printf("Initialization completed...\n");

	const Vec_3D<double> pos_0 (sim.mesh_num/2., sim.mesh_num/2., sim.mesh_num/2.);
	Vec_3D<double> pos;
	const Vec_3D<double> vel (0., 0., 0.);
	
	for(int i = 0; i < sim.par_num; i++)
	{
		for(int j = 0; j < 3; j++) pos[j] = pos_0[j]+0.01*rand() / RAND_MAX -0.005;
		APP.particles[i] = Particle_v(pos, vel);
	}
	
	get_rho_from_par(APP.particles, &APP.power_aux, sim);
	print_projected_rho(APP.power_aux, sim, out_dir_app, "_0");
	printf("Computing density in k-sapce...\n");
	fftw_execute_dft_r2c(APP.p_F, APP.power_aux);
	
	/* Computing initial potential in k-space */
	gen_pot_k(APP.power_aux, &APP.power_aux);
	
	/* Computing displacement in k-space with S2 shaped particles */
	gen_displ_k_S2(&APP.app_field, APP.power_aux, sim.a);
	
	/* Computing force in q-space */
	printf("Computing force for S2 shaped particles in q-space...\n");
	fftw_execute_dft_c2r_triple(APP.p_B, APP.app_field);
	
	printf("Creating linked list...\n");
	APP.linked_list.get_linked_list(APP.particles);
	
	print_par_pos_cut(APP.particles, sim, out_dir_app, APP.z_suffix());
//	print_rho_map(APP.power_aux, sim, out_dir_app, APP.z_suffix());
	print_force(APP.particles, APP.app_field,sim, out_dir_app, APP.z_suffix(), APP.linked_list);
	

	
	printf("Example ended successfully.\n");
	return 0;
}