#include "app_var.hpp"
#include "core_out.h"
#include "core_app.h"
#include "core_mesh.h"

using namespace std;

static const char *humanSize(uint64_t bytes){
	char const *suffix[] = {"B", "KB", "MB", "GB", "TB"};
	char length = sizeof(suffix) / sizeof(suffix[0]);

	int i = 0;
	double dblBytes = bytes;

	if (bytes > 1024) {
		for (i = 0; (bytes / 1024) > 0 && i<length-1; i++, bytes /= 1024)
			dblBytes = bytes / 1024.0;
	}

	static char output[200];
	sprintf(output, "%.02lf %s", dblBytes, suffix[i]);
	return output;
}

/**
 * @class:	Tracking
 * @brief:	class storing info about tracked particles
 */


Tracking::Tracking(int sqr_num_track_par, int par_num_per_dim):
sqr_num_track_par(sqr_num_track_par), num_track_par(sqr_num_track_par*sqr_num_track_par)
{
    printf("Initializing IDs of tracked particles...\n");
    par_ids.reserve(num_track_par);
    int x, y, z;
    FTYPE_t s;
    y = par_num_per_dim / 2; // middle of the cube
    s = par_num_per_dim / FTYPE_t(4*(sqr_num_track_par+1)); // quarter of the cube
    for (int i=1; i<=sqr_num_track_par;i++)
    {
        z = (int)(s*i);
        for (int j=1; j<=sqr_num_track_par;j++)
        {
            x = (int)(s*j);
            par_ids.push_back(x*par_num_per_dim*par_num_per_dim+y*par_num_per_dim+z);
        }
    }
}

template <class T> void Tracking::update_track_par(const std::vector<T>& particles)
{
    std::vector<Particle_x<FTYPE_t>> par_pos_step;
    par_pos_step.reserve(num_track_par);
    for (int i=0; i<num_track_par; i++){
        par_pos_step.emplace_back(particles[par_ids[i]].position);
    }
    par_pos.push_back(par_pos_step);
}

void Tracking::print_track_par(const Sim_Param &sim, string out_dir, string suffix)
{
    out_dir += "par_cut/";
    string file_name = out_dir + "track_par_pos" + suffix + ".dat";
    Ofstream File(file_name);

    FTYPE_t x,y,z;
    const FTYPE_t x_0 = sim.x_0();
    cout << "Writing positons of " << num_track_par << " tracked particles into file " << file_name << "\n";
    File << "# This file contains positions of particles in units [Mpc/h].\n"
            "# x [Mpc/h]\tz [Mpc/h]\n";
    for (int i=0; i<num_track_par; i++){
        for (unsigned j=0; j<num_step();j++){
            x = par_pos[j][i].position[0];
            y = par_pos[j][i].position[1];
            z = par_pos[j][i].position[2];
            File << x*x_0 << "\t" << z*x_0 << "\t" << y*x_0 << "\n";
        }
        File << "\n\n";
    }
}

/**
 * @class:	App_Var<T>
 * @brief:	class containing variables for approximations
 */
 
template <class T> 
App_Var<T>::App_Var(const Sim_Param &sim, string app_str):
	sim(sim), step(0), print_every(sim.out_opt.print_every),
    b(sim.integ_opt.b_in), b_out(sim.integ_opt.b_out), db(sim.integ_opt.db),
    app_str(app_str), z_suffix_const("_" + app_str + "_"), out_dir_app(std_out_dir(app_str + "_run/", sim)),
	track(4, sim.box_opt.par_num_1d),
    dens_binned(500), is_init_pwr_spec_0(false), is_init_vel_pwr_spec_0(false)
{    
    // EFFICIENTLY ALLOCATE VECTOR OF MESHES
    app_field.reserve(3);
    power_aux.reserve(3);
    for(size_t i = 0; i < 3; i++){
        app_field.emplace_back(sim.box_opt.mesh_num);
        power_aux.emplace_back(sim.box_opt.mesh_num_pwr);
    }
    memory_alloc = sizeof(FTYPE_t)*(app_field[0].length*app_field.size()+power_aux[0].length*power_aux.size());

    // RESERVE MEMORY FOR BINNED POWER SPECTRA
    unsigned bin_num = (unsigned)ceil(log10(sim.box_opt.mesh_num_pwr)*sim.out_opt.bins_per_decade);
    pwr_spec_binned.reserve(bin_num);
    pwr_spec_binned_0.reserve(bin_num);
    vel_pwr_spec_binned_0.reserve(bin_num);

    // RESERVE MEMORY FOR BINNED CORRELATION FUNCTION
    bin_num =  (unsigned)ceil((sim.other_par.x_corr.upper - sim.other_par.x_corr.lower)/ 10. * sim.out_opt.points_per_10_Mpc);
    corr_func_binned.reserve(bin_num);
    
    // CREAT SUBDIR STRUCTURE
    if (print_every) work_dir_over(out_dir_app);
    else create_dir(out_dir_app);

    // PARTICLES ALLOCATION
    particles.reserve(sim.box_opt.par_num);
    memory_alloc += sizeof(T)*sim.box_opt.par_num;

	// FFTW PREPARATION
	if (!FFTW_PLAN_OMP_INIT()){
		throw runtime_error("Errors during multi-thread initialization");
	}
	FFTW_PLAN_OMP(sim.run_opt.nt);
	p_F = FFTW_PLAN_R2C(sim.box_opt.mesh_num, sim.box_opt.mesh_num, sim.box_opt.mesh_num, app_field[0].real(),
        app_field[0].complex(), FFTW_ESTIMATE);
	p_B = FFTW_PLAN_C2R(sim.box_opt.mesh_num, sim.box_opt.mesh_num, sim.box_opt.mesh_num, app_field[0].complex(),
        app_field[0].real(), FFTW_ESTIMATE);
    p_F_pwr = FFTW_PLAN_R2C(sim.box_opt.mesh_num_pwr, sim.box_opt.mesh_num_pwr, sim.box_opt.mesh_num_pwr, power_aux[0].real(),
		power_aux[0].complex(), FFTW_ESTIMATE);
	p_B_pwr = FFTW_PLAN_C2R(sim.box_opt.mesh_num_pwr, sim.box_opt.mesh_num_pwr, sim.box_opt.mesh_num_pwr, power_aux[0].complex(),
		power_aux[0].real(), FFTW_ESTIMATE);
}

template <class T> 
App_Var<T>::~App_Var()
{	// FFTW CLEANUP
	FFTW_DEST_PLAN(p_F);
    FFTW_DEST_PLAN(p_B);
    FFTW_DEST_PLAN(p_F_pwr);
	FFTW_DEST_PLAN(p_B_pwr);
	FFTW_PLAN_OMP_CLEAN();
}

template <class T> 
string App_Var<T>::z_suffix()
{
	z_suffix_num.str("");
	z_suffix_num << fixed << setprecision(2) << z();
	return z_suffix_const + "z" + z_suffix_num.str();
}

template <class T> 
void App_Var<T>::upd_time()
{
	step++;
	if ((b_out - b) < db) db = b_out - b;
	b += db;
}

template <class T> 
void App_Var<T>::print_position()
{/* Printing positions */
    print_par_pos_cut_small(particles, sim, out_dir_app, z_suffix());
    track.print_track_par(sim, out_dir_app, z_suffix());
}

template <class T> 
void App_Var<T>::print_density()
{/* Printing density */
    gen_dens_binned(power_aux[0], dens_binned, sim);    
    print_rho_map(power_aux[0], sim, out_dir_app, z_suffix());
    print_dens_bin(dens_binned, out_dir_app, z_suffix());
}

template <class T> 
void App_Var<T>::get_binned_power_spec()
{/* Compute power spectrum and bin it */
    fftw_execute_dft_r2c(p_F_pwr, power_aux[0]);
    pwr_spec_k(power_aux[0], power_aux[0]);
    gen_pow_spec_binned(sim, power_aux[0], pwr_spec_binned);
}

template <class T> 
void App_Var<T>::print_power_spec()
{/* Printing power spectrum */
    print_pow_spec(pwr_spec_binned, out_dir_app, "_par" + z_suffix());
    if (!is_init_pwr_spec_0){
        pwr_spec_binned_0 = pwr_spec_binned;
        D_init = growth_factor(b, sim.cosmo);
        is_init_pwr_spec_0 = true;
    }
    FTYPE_t D_now = growth_factor(b, sim.cosmo);
    print_pow_spec_diff(pwr_spec_binned, pwr_spec_binned_0, D_now / D_init, out_dir_app, "_par" + z_suffix());
    print_pow_spec_diff(pwr_spec_binned, pwr_spec_input, D_now, out_dir_app, "_input" + z_suffix());
    print_pow_spec_diff(pwr_spec_binned, pwr_spec_binned_0, pwr_spec_input, D_now, D_init,
                        out_dir_app, "_hybrid" + z_suffix());
}

template <class T> 
void App_Var<T>::print_extrap_pwr(const Extrap_Pk<FTYPE_t, 2>& P_k)
{/* Print extrapolated power spectrum */
    gen_pow_spec_binned_from_extrap(sim, P_k, pwr_spec_binned);
    print_pow_spec(pwr_spec_binned, out_dir_app, "_extrap" + z_suffix());
}

template <class T> 
void App_Var<T>::print_corr(const Extrap_Pk<FTYPE_t, 2>& P_k)
{/* Printing correlation function */
    gen_corr_func_binned_gsl_qawf(sim, P_k, corr_func_binned);
    print_corr_func(corr_func_binned, out_dir_app, "_gsl_qawf_par" + z_suffix());
    gen_corr_func_binned_gsl_qawf_lin(sim, b, corr_func_binned);
    print_corr_func(corr_func_binned, out_dir_app, "_gsl_qawf_lin" + z_suffix());
}

template <class T> 
void App_Var<T>::print_vel_pwr()
{/* Print velocity power spectrum */
    fftw_execute_dft_r2c_triple(p_F_pwr, power_aux);
    vel_pwr_spec_k(power_aux, power_aux[0]);
    gen_pow_spec_binned(sim, power_aux[0], pwr_spec_binned);
    print_vel_pow_spec(pwr_spec_binned, out_dir_app, z_suffix());
    if (!is_init_vel_pwr_spec_0){
        vel_pwr_spec_binned_0 = pwr_spec_binned;
        is_init_vel_pwr_spec_0 = true;
        dDda_init = growth_change(b, sim.cosmo);
    }
    print_vel_pow_spec_diff(pwr_spec_binned, vel_pwr_spec_binned_0, growth_change(b, sim.cosmo) / dDda_init, out_dir_app, z_suffix());
}

template <class T> 
void App_Var<T>::print_output()
{
    /* Printing positions */
    if (sim.out_opt.print_par_pos) print_position();

    /* Get discrete density from particles */
    if (sim.out_opt.get_rho) get_rho_from_par(particles, power_aux[0], sim);
    
    /* Printing density */
    if (sim.out_opt.print_dens) print_density();

    /* Compute power spectrum and bin it */
    if (sim.out_opt.get_pwr) get_binned_power_spec();

    /* Printing power spectrum */
    if (sim.out_opt.print_pwr) print_power_spec();

    /* Extrapolate power spectrum beyond range of simulation box */
    if (sim.out_opt.get_pk_extrap){
        const Extrap_Pk<FTYPE_t, 2> P_k(pwr_spec_binned, sim);
        /* Print extrapolated power spectrum */
        if (sim.out_opt.print_extrap_pwr) print_extrap_pwr(P_k);
        /* Printing correlation function */
        if (sim.out_opt.print_corr) print_corr(P_k);
    }

    /* Velocity power spectrum */
    if (sim.out_opt.print_vel_pwr && get_vel_from_par(particles, power_aux, sim)) print_vel_pwr();
}

template <class T> 
void App_Var<T>::print_mem() const
{
    printf("Allocated %s of memory.\n", humanSize(memory_alloc));
}

template <class T> 
void App_Var<T>::print_info() const
{
    sim.print_info(out_dir_app, app_str);
}

/**
 * @class:	App_Var_AA
 * @brief:	class containing variables for adhesion approximation
 */
 
 App_Var_AA::App_Var_AA(const Sim_Param &sim, string app_str):
    App_Var<Particle_v<FTYPE_t>>(sim, app_str), expotential (sim.box_opt.mesh_num)
{
    memory_alloc += sizeof(FTYPE_t)*expotential.length;
}

/**
 * @class:	App_Var_FP_mod
 * @brief:	class containing variables for modified Frozen-potential approximation
 */
 
 App_Var_FP_mod::App_Var_FP_mod(const Sim_Param &sim, string app_str):
    App_Var<Particle_v<FTYPE_t>>(sim, app_str), linked_list(sim.box_opt.par_num, sim.app_opt.M, sim.app_opt.Hc)
{
    memory_alloc += sizeof(int)*(linked_list.HOC.length+linked_list.par_num);

    // precompute short range force
    size_t res = size_t(sim.app_opt.rs/0.05)+1; // force resolution 5% of mesh cell
    const FTYPE_t r0 = sim.app_opt.rs / (res-1);
    Data_Vec<FTYPE_t, 2> data(res);
    FTYPE_t r;
    const FTYPE_t e2 = pow2(sim.box_opt.Ng*0.1); // softening of 10% of average interparticle length

    #pragma omp parallel for private(r)
    for(unsigned i = 0; i < res; i++)
    {
        r = i*r0;
        data[0][i] = pow2(r); // store square of r
        data[1][i] = (force_tot(r, e2) - force_ref(r, sim.app_opt.a))/(4*PI);
    }
    fs_interp.init(data);
}

/**
 * @class LinkedList
 * @brief class handling linked lists
 */


LinkedList::LinkedList(unsigned par_num, int m, FTYPE_t hc):
	par_num(par_num), Hc(hc), LL(par_num), HOC(m, m, m) {}
	
void LinkedList::get_linked_list(const std::vector<Particle_v<FTYPE_t>>& particles)
{
	HOC.assign(-1);
	for (unsigned i = 0; i < par_num; i++)
	{
		LL[i] = HOC(particles[i].position/Hc);
		HOC(particles[i].position/Hc) = i;
	}
}

template void Tracking::update_track_par(const std::vector<Particle_x<FTYPE_t>>& particles);
template void Tracking::update_track_par(const std::vector<Particle_v<FTYPE_t>>& particles);
template class App_Var<Particle_x<FTYPE_t>>;
template class App_Var<Particle_v<FTYPE_t>>;

#ifdef TEST

#include <catch.hpp>
#include "test.hpp"

TEST_CASE( "UNIT TEST: tracking class {Tracking}", "[core]" )
{
    print_unit_msg("tracking class {Tracking}");

    // 2, 12, 2
    Tracking track(2, 24);

    CHECK( track.num_track_par == 4 );
    CHECK( track.par_ids.size() == 4 );
    CHECK( track.par_ids[0] == 1442 );
}

#endif