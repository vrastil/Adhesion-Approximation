/**
 * @file:	core.cpp
 * @brief:	class definitions
 */
 
#include "core.h"
#include "core_cmd.h"
#include "core_out.h"
#include "core_app.h"
#include "core_mesh.h"
#include "core_power.h"

using namespace std;
using json = nlohmann::json;

const char *humanSize(uint64_t bytes){
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
 * @class:	Cosmo_Param
 * @brief:	class storing parameters for power spectrum
 */

Cosmo_Param::Cosmo_Param():
    // cosmo == NULL as indicator of uninitialization
    // config first initialize to default (in case new configuration options are added)
    config(default_config), cosmo(NULL)
    {
    #ifdef TEST
    cout << ">>> Debug: Creating Cosmo_Param via call to Cosmo_Param()\n";
    cout << "\tconfig.transfer_function_method = " << config.transfer_function_method << "\n";
    #endif
    }

void Cosmo_Param::init()
{
    k2_G *= k2_G;
    h = H0/100;
    #ifdef TEST
    cout << ">>> Debug: Creating Cosmo_Param via call to init()\n";
    cout << "\tconfig.transfer_function_method = " << config.transfer_function_method << "\n";
    #endif

    // create flat LCDM cosmology
    int status = 0;
    cosmo = ccl_cosmology_create_with_lcdm_params(Omega_c(), Omega_b, 0, h, sigma8, ns, config, &status);
    if (status) throw std::runtime_error(cosmo->status_message);

    #ifdef TEST
    cout << ">>> Debug: Created ccl_cosmology*: " << cosmo << "\n";
    #endif
    
    // PRECOMPUTED VALUES
    D_norm = norm_growth_factor(*this); //< use only when outside CCL range

    // normalize power spectrum
    norm_pwr(*this);
}


Cosmo_Param::~Cosmo_Param()
{
    #ifdef TEST
    cout << ">>> Debug: Cosmo_Param destructor: " << this << "\n";
    #endif
    if(cosmo){
        #ifdef TEST
        cout << ">>> Debug: Free space after ccl_cosmology: " << cosmo << "\n";
        #endif
        ccl_cosmology_free(cosmo);
    }
}

Cosmo_Param::operator void*() const
{
    return const_cast<Cosmo_Param*>(this);
}


// convert to pyccl transfer_function_types keys
const map<string, transfer_function_t> transfer_function_method = {
    {"emulator", ccl_emulator},
    {"eisenstein_hu", ccl_eisenstein_hu},
    {"bbks", ccl_bbks},
    {"boltzmann_class", ccl_boltzmann_class},
    {"boltzmann_camb", ccl_boltzmann_camb}
};
// convert to pyccl matter_power_spectrum_types keys
const map<string, matter_power_spectrum_t> matter_power_spectrum_method = {
    {"linear", ccl_linear},
    {"halofit", ccl_halofit},
    {"halo_model", ccl_halo_model}
};
// convert to pyccl mass_function_types keys
const map<string, mass_function_t> mass_function_method = {
    {"tinker", ccl_tinker},
    {"tinker10", ccl_tinker10},
    {"watson", ccl_watson},
    {"angulo", ccl_angulo}
};
// convert to pyccl baryons_power_spectrum keys
const map<string, baryons_power_spectrum_t> baryons_power_spectrum_method = {
    {"nobaryons", ccl_nobaryons},
    {"bcm", ccl_bcm}
};

/**
 * return first occurence of 'value' in std::map
 */
template<typename T, typename U>
T find_value(const std::map<T, U>& map, const U& value)
{
    for(auto x : map) if (x.second == value) return x.first;
    throw std::out_of_range("Value not found");
}

void to_json(json& j, const Cosmo_Param& cosmo)
{
    j = json{
        {"A", cosmo.A},
        {"index", cosmo.ns},
        {"sigma8", cosmo.sigma8},
        {"smoothing_k", cosmo.k2_G},
        {"Omega_c", cosmo.Omega_c()},
        {"Omega_b", cosmo.Omega_b},
        {"Omega_m", cosmo.Omega_m},
        {"h", cosmo.h}
    };
    
    j["transfer_function_method"] = find_value(transfer_function_method, cosmo.config.transfer_function_method);
    j["matter_power_spectrum_method"] = find_value(matter_power_spectrum_method, cosmo.config.matter_power_spectrum_method);
    j["mass_function_method"] = find_value(mass_function_method, cosmo.config.mass_function_method);
    j["baryons_power_spectrum_method"] = find_value(baryons_power_spectrum_method, cosmo.config.baryons_power_spectrum_method);
}

void from_json(const json& j, Cosmo_Param& cosmo)
{
    #ifdef TEST
    cout << ">>> Debug: Loading 'Cosmo_Param& cosmo' from json file\n";
    #endif
    cosmo.A = j.at("A").get<FTYPE>();
    cosmo.ns = j.at("index").get<FTYPE>();
    cosmo.sigma8 = j.at("sigma8").get<FTYPE>();
    cosmo.k2_G = j.at("smoothing_k").get<FTYPE>();
    cosmo.Omega_b = j.at("Omega_b").get<FTYPE>();
    cosmo.Omega_m = j.at("Omega_m").get<FTYPE>();
    cosmo.h = j.at("h").get<FTYPE>();
    cosmo.H0 = cosmo.h * 100;
    
    string tmp;
    try{
        tmp = j.at("transfer_function_method").get<string>();
        cosmo.config.transfer_function_method = transfer_function_method.at(tmp);
    }catch(const out_of_range& oor){
        cosmo.config.transfer_function_method = ccl_boltzmann_class;
    }
    try{
        tmp = j.at("matter_power_spectrum_method").get<string>();
        cosmo.config.matter_power_spectrum_method = matter_power_spectrum_method.at(tmp);
    }catch(const out_of_range& oor){
        cosmo.config.matter_power_spectrum_method = ccl_halofit;
    }
    try{
        tmp = j.at("mass_function_method").get<string>();
        cosmo.config.mass_function_method = mass_function_method.at(tmp);
    }catch(const out_of_range& oor){
        cosmo.config.mass_function_method = ccl_tinker10;
    }
    try{
        tmp = j.at("baryons_power_spectrum_method").get<string>();
        cosmo.config.baryons_power_spectrum_method = baryons_power_spectrum_method.at(tmp);
    }catch(const out_of_range& oor){
        cosmo.config.baryons_power_spectrum_method = ccl_nobaryons;
    }
    cosmo.init();
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
	FTYPE s;
	y = par_num_per_dim / 2; // middle of the cube
	s = par_num_per_dim / FTYPE(4*(sqr_num_track_par+1)); // quarter of the cube
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

 template <class T>
 void Tracking::update_track_par(T* particles)
 {
     std::vector<Particle_x> par_pos_step;
     par_pos_step.reserve(num_track_par);
     for (int i=0; i<num_track_par; i++){
         par_pos_step.emplace_back(particles[par_ids[i]].position);
     }
     par_pos.push_back(par_pos_step);
 }

void Run_Opt::init()
{
    if(nt == 0) nt = omp_get_max_threads();
    else omp_set_num_threads(nt);
    if (seed == 0){
        srand(time(NULL));
        seed = (static_cast<long>(rand()) << (sizeof(int) * 8)) | rand();
    } else mlt_runs = 1;
    phase = true;
}

void Box_Opt::init()
{
    Ng = mesh_num / par_num_1d;
    Ng_pwr = mesh_num_pwr/par_num_1d;
    par_num = par_num_1d*par_num_1d*par_num_1d;
}

void Integ_Opt::init()
{
    b_in = 1/(z_in + 1);
	b_out = 1/(z_out + 1);
}

void Out_Opt::init()
{
    get_pk_extrap = print_corr|| print_extrap_pwr;
    get_pwr = get_pk_extrap || print_pwr;
    get_rho = get_pwr || print_dens;
}

void App_Opt::init(const Box_Opt& box_opt)
{
    a = rs / FTYPE(0.735);
    M = (int)(box_opt.mesh_num / rs);
    Hc = FTYPE(box_opt.mesh_num) / M;
    nu_dim = nu;
    nu /= pow_(box_opt.box_size/box_opt.mesh_num, 2); // converting to dimensionless units
}

void Other_par::init(const Box_Opt& box_opt)
{
    FTYPE tmp = PI/box_opt.box_size;

    nyquist["analysis"] = tmp*box_opt.mesh_num_pwr;
    nyquist["potential"] = tmp*box_opt.mesh_num;
    nyquist["particle"] = tmp*box_opt.par_num_1d;
    k_print.lower = 2*tmp;
    k_print.upper = 2*tmp*box_opt.mesh_num_pwr;
    x_corr.lower = 0.1;
    x_corr.upper = 200;
}

/**
 * @class:	Sim_Param
 * @brief:	class storing simulation parameters
 */

Sim_Param::Sim_Param(int ac, char* av[])
{
	handle_cmd_line(ac, av, *this);//< throw if anything happend
    run_opt.init();
    box_opt.init();
    integ_opt.init();
    out_opt.init();
    app_opt.init(box_opt);
    cosmo.init();
    other_par.init(box_opt);
}

Sim_Param::Sim_Param(string file_name)
{
    try{
        Ifstream i(file_name);
        json j;
        i >> j;
        try{ run_opt = j.at("run_opt"); } // sim_param.json has run_opt
        catch(const out_of_range& oor){ // stack_info.json does not have run_opt
            run_opt.nt = 0; // max
            run_opt.seed = 0; // random
            run_opt.init();
        }
        box_opt = j.at("box_opt");
        integ_opt = j.at("integ_opt");
        try{ out_opt = j.at("out_opt"); } // new format of json files
        catch(const out_of_range& oor){ // old format does not store Out_Opt
            try {out_opt.out_dir = j.at("out_dir"); } // stack_info.json doesn`t store out_dir
            catch(const out_of_range& oor){out_opt.out_dir = "~/home/FASTSIM/output/"; } // do not need it, set to some default
            out_opt.bins_per_decade = 20;
            out_opt.points_per_10_Mpc = 10;
        }
        app_opt = j.at("app_opt");
        from_json(j.at("cosmo"), cosmo); //< call explicitly, SWIG has some issues with json.hpp
        app_opt.init(box_opt);
        other_par.init(box_opt);
    }
    catch(const out_of_range& oor){
        string err = string(oor.what()) + " in file '" + file_name + "'";
        throw out_of_range(err);
    }
}

void to_json(json& j, const Box_Opt& box_opt)
{
    j = json{
        {"mesh_num", box_opt.mesh_num},
        {"mesh_num_pwr", box_opt.mesh_num_pwr},
        {"Ng", box_opt.Ng},
        {"par_num", box_opt.par_num_1d},
        {"box_size", box_opt.box_size}
    };
}

void from_json(const json& j, Box_Opt& box_opt)
{
    box_opt.mesh_num = j.at("mesh_num").get<unsigned>();
    box_opt.mesh_num_pwr = j.at("mesh_num_pwr").get<unsigned>();
    box_opt.par_num_1d = j.at("par_num").get<unsigned>();
    box_opt.box_size = j.at("box_size").get<FTYPE>();

    box_opt.init();
}

void to_json(json& j, const Integ_Opt& integ_opt)
{
    j = json{
        {"redshift", integ_opt.z_in},
        {"redshift_0", integ_opt.z_out},
        {"time_step", integ_opt.db}
    };
}

void from_json(const json& j, Integ_Opt& integ_opt)
{
    integ_opt.z_in = j.at("redshift").get<FTYPE>();
    integ_opt.z_out = j.at("redshift_0").get<FTYPE>();
    integ_opt.db = j.at("time_step").get<FTYPE>();

    integ_opt.init();
}

void to_json(json& j, const App_Opt& app_opt)
{
    j = json{
        {"viscosity", app_opt.nu_dim},
        {"cut_radius", app_opt.rs}
    };
}

void from_json(const json& j, App_Opt& app_op)
{
    app_op.nu_dim = j.at("viscosity").get<FTYPE>();
    app_op.nu = app_op.nu_dim;
    app_op.rs = j.at("cut_radius").get<FTYPE>();
}

void to_json(json& j, const Run_Opt& run_opt)
{
    j = json{
        {"num_thread", run_opt.nt},
        {"seed", run_opt.seed}
    };
}

void from_json(const json& j, Run_Opt& run_opt)
{
    run_opt.nt = j.at("num_thread").get<unsigned>();
    run_opt.seed = j.at("seed").get<unsigned long>();
    run_opt.init();
}

void to_json(json& j, const Out_Opt& out_opt)
{
    j = json{
        {"bins_per_decade", out_opt.bins_per_decade},
        {"points_per_10_Mpc", out_opt.points_per_10_Mpc},
        {"out_dir", out_opt.out_dir}
    };
}

void from_json(const json& j, Out_Opt& out_opt)
{
    out_opt.bins_per_decade = j.at("bins_per_decade").get<unsigned>();
    out_opt.points_per_10_Mpc = j.at("points_per_10_Mpc").get<unsigned>();
    out_opt.out_dir = j.at("out_dir").get<string>();
}

void Sim_Param::print_info(string out, string app) const
{
    if (out == "")
    {
        printf("\n*********************\n");
        printf("SIMULATION PARAMETERS\n");
        printf("*********************\n");
        printf("Ng:\t\t%i\n", box_opt.Ng);
        printf("Num_par:\t%i^3\n", box_opt.par_num_1d);
        printf("Num_mesh:\t%i^3\n", box_opt.mesh_num);
        printf("Num_mesh_pwr:\t%i^3\n", box_opt.mesh_num_pwr);
        printf("Box size:\t%.0f Mpc/h\n", box_opt.box_size);
        printf("Redshift:\t%G--->%G\n", integ_opt.z_in, integ_opt.z_out);
        printf("Pk:\t\t[sigma_8 = %G, As = %G, ns = %G, k_smooth = %G]\n", 
            cosmo.sigma8, cosmo.A, cosmo.ns, sqrt(cosmo.k2_G));
        cout <<"\t\t[transfer_function_method = " << find_value(transfer_function_method, cosmo.config.transfer_function_method) << "]\n";
        cout <<"\t\t[matter_power_spectrum_method = " << find_value(matter_power_spectrum_method, cosmo.config.matter_power_spectrum_method) << "]\n";
        cout <<"\t\t[mass_function_method = " << find_value(mass_function_method, cosmo.config.mass_function_method) << "]\n";
        cout << "\t\t[baryons_power_spectrum_method = " << find_value(baryons_power_spectrum_method, cosmo.config.baryons_power_spectrum_method) << "]\n";
        printf("AA:\t\t[nu = %G (Mpc/h)^2]\n", app_opt.nu_dim);
        printf("LL:\t\t[rs = %G, a = %G, M = %i, Hc = %G]\n", app_opt.rs, app_opt.a, app_opt.M, app_opt.Hc);
        printf("num_thread:\t%i\n", run_opt.nt);
        printf( "Output:\t\t'%s'\n", out_opt.out_dir.c_str());
    }
    else
    {
        string file_name = out + "sim_param.json";
        ofstream o(file_name);

        json j = {
            {"box_opt", box_opt},
            {"integ_opt", integ_opt},
            {"cosmo", cosmo},
            {"app_opt", app_opt},
            {"run_opt", run_opt},
            {"out_opt", out_opt},
            {"k_nyquist", other_par.nyquist},
            {"results", {}},
            {"app", app}
        };
        o << setw(2) << j << endl;
        o.close();
    }
}


void Sim_Param::print_info() const
{
	Sim_Param::print_info("", "");
}

bool Run_Opt::simulate()
{
    if (!pair || !phase)
    {
        mlt_runs--;
        seed = (static_cast<long>(rand()) << (sizeof(int) * 8)) | rand();
    }
    if (pair) phase = !phase;
    return mlt_runs;
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
    memory_alloc = sizeof(FTYPE)*(app_field[0].length*app_field.size()+power_aux[0].length*power_aux.size());

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

    // PARTICLES INITIALIZATION
    particles = new T[sim.box_opt.par_num];
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
{
    delete[] particles;

	// FFTW CLEANUP
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
void App_Var<T>::print_output()
{
    /* Printing positions */
    if (sim.out_opt.print_par_pos){
        print_par_pos_cut_small(particles, sim, out_dir_app, z_suffix());
        print_track_par(track, sim, out_dir_app, z_suffix());
    }

    /* Get discrete density from particles */
    if (sim.out_opt.get_rho)
        get_rho_from_par(particles, &power_aux[0], sim);
    
    /* Printing density */
    if (sim.out_opt.print_dens){
        gen_dens_binned(power_aux[0], dens_binned, sim);    
        print_rho_map(power_aux[0], sim, out_dir_app, z_suffix());
        print_dens_bin(dens_binned, out_dir_app, z_suffix());
    }

    /* Compute power spectrum and bin it */
    if (sim.out_opt.get_pwr){
        fftw_execute_dft_r2c(p_F_pwr, power_aux[0]);
        pwr_spec_k(power_aux[0], &power_aux[0]);
        gen_pow_spec_binned(sim, power_aux[0], &pwr_spec_binned);
    }

    /* Printing power spectrum */
    if (sim.out_opt.print_pwr){
        print_pow_spec(pwr_spec_binned, out_dir_app, "_par" + z_suffix());
        if (!is_init_pwr_spec_0){
            pwr_spec_binned_0 = pwr_spec_binned;
            D_init = growth_factor(b, sim.cosmo);
            is_init_pwr_spec_0 = true;
        }
        FTYPE D_now = growth_factor(b, sim.cosmo);
        print_pow_spec_diff(pwr_spec_binned, pwr_spec_binned_0, D_now / D_init, out_dir_app, "_par" + z_suffix());
        print_pow_spec_diff(pwr_spec_binned, pwr_spec_input, D_now, out_dir_app, "_input" + z_suffix());
        print_pow_spec_diff(pwr_spec_binned, pwr_spec_binned_0, pwr_spec_input, D_now, D_init,
                            out_dir_app, "_hybrid" + z_suffix());
    }

    /* Extrapolate power spectrum beyond range of simulation box */
    if (sim.out_opt.get_pk_extrap){
        Extrap_Pk<FTYPE, 2> P_k(pwr_spec_binned, sim);
    /* Print extrapolated power spectrum */
        if (sim.out_opt.print_extrap_pwr){
            gen_pow_spec_binned_from_extrap(sim, P_k, &pwr_spec_binned);
            print_pow_spec(pwr_spec_binned, out_dir_app, "_extrap" + z_suffix());
        }
    /* Printing correlation function */
        if (sim.out_opt.print_corr){
            gen_corr_func_binned_gsl_qawf(sim, P_k, corr_func_binned);
            print_corr_func(corr_func_binned, out_dir_app, "_gsl_qawf_par" + z_suffix());
            gen_corr_func_binned_gsl_qawf_lin(sim, b, corr_func_binned);
            print_corr_func(corr_func_binned, out_dir_app, "_gsl_qawf_lin" + z_suffix());
        }
    }

    /* Velocity power spectrum */
    if (sim.out_opt.print_vel_pwr && get_vel_from_par(particles, &power_aux, sim)){
        fftw_execute_dft_r2c_triple(p_F_pwr, power_aux);
        vel_pwr_spec_k(power_aux, &power_aux[0]);
        gen_pow_spec_binned(sim, power_aux[0], &pwr_spec_binned);
        print_vel_pow_spec(pwr_spec_binned, out_dir_app, z_suffix());
        if (!is_init_vel_pwr_spec_0){
            vel_pwr_spec_binned_0 = pwr_spec_binned;
            is_init_vel_pwr_spec_0 = true;
            dDda_init = growth_change(b, sim.cosmo);
        }
        print_vel_pow_spec_diff(pwr_spec_binned, vel_pwr_spec_binned_0, growth_change(b, sim.cosmo) / dDda_init, out_dir_app, z_suffix());
    }
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
    App_Var<Particle_v>(sim, app_str), expotential (sim.box_opt.mesh_num)
{
    memory_alloc += sizeof(FTYPE)*expotential.length;
}

/**
 * @class:	App_Var_FP_mod
 * @brief:	class containing variables for modified Frozen-potential approximation
 */
 
 App_Var_FP_mod::App_Var_FP_mod(const Sim_Param &sim, string app_str):
    App_Var<Particle_v>(sim, app_str), linked_list(sim.box_opt.par_num, sim.app_opt.M, sim.app_opt.Hc)
{
    memory_alloc += sizeof(int)*(linked_list.HOC.length+linked_list.par_num);

    // precompute short range force
    size_t res = size_t(sim.app_opt.rs/0.05)+1; // force resolution 5% of mesh cell
    const FTYPE r0 = sim.app_opt.rs / (res-1);
    Data_Vec<FTYPE, 2> data(res);
    FTYPE r;
    const FTYPE e2 = pow_(sim.box_opt.Ng*0.1, 2); // softening of 10% of average interparticle length

    #pragma omp parallel for private(r)
    for(unsigned i = 0; i < res; i++)
    {
        r = i*r0;
        data[0][i] = pow_(r, 2); // store square of r
        data[1][i] = (force_tot(r, e2) - force_ref(r, sim.app_opt.a))/(4*PI);
    }
    fs_interp.init(data);
}

/**
 * @class LinkedList
 * @brief class handling linked lists
 */


LinkedList::LinkedList(unsigned par_num, int m, FTYPE hc):
	par_num(par_num), Hc(hc), LL(par_num), HOC(m, m, m) {}
	
void LinkedList::get_linked_list(Particle_v* particles)
{
	HOC.assign(-1);
	for (unsigned i = 0; i < par_num; i++)
	{
		LL[i] = HOC(particles[i].position/Hc);
		HOC(particles[i].position/Hc) = i;
	}
}

template void Tracking::update_track_par(Particle_x* particles);
template void Tracking::update_track_par(Particle_v* particles);
template class App_Var<Particle_x>;
template class App_Var<Particle_v>;

#ifdef TEST
#include "test_core.cpp"
#endif