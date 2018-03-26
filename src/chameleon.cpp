#include "core_app.h"
#include "core_power.h"
#include "core_mesh.h"
#include "core_out.h"
#include "chameleon.hpp"

using namespace std;

namespace{
constexpr FTYPE_t MPL = (FTYPE_t)1; // chi in units of Planck mass
// (FTYPE_t)2.435E18; // reduced Planck mass, [GeV/c^2]
// constexpr FTYPE_t fm_to_Mpc = (FTYPE_t)3.2407792896664E-38; // 1 fm = ? Mpc
// constexpr FTYPE_t hbarc = (FTYPE_t)197.327053; // reduced Planck constant times speed of light, [MeV fm]
// constexpr FTYPE_t hbarc_cosmo = hbarc*FTYPE_t(1E-9)*fm_to_Mpc; // [GeV Mpc]
// constexpr FTYPE_t G_N = hbarc_cosmo*FTYPE_t(6.70711*1E-39); // gravitational constant, [GeV Mpc / (GeV/c^2)^2]
constexpr FTYPE_t c_kms = (FTYPE_t)299792.458; // speed of light [km / s]

// compute chi in units of chi_a
#define CHI_A_UNITS

template<typename T>
void transform_Mesh_to_Grid(const Mesh& mesh, Grid<3, T> &grid)
{/* copy data in Mesh 'N*N*(N+2)' onto MultiGrid 'N*N*N' */
    unsigned int ix, iy, iz;
    const unsigned N_tot = grid.get_Ntot();
    const unsigned N = grid.get_N();

    if (mesh.N != N) throw std::range_error("Mesh of a different size than Grid!");

    #pragma omp parallel for private(ix, iy, iz)
    for (unsigned i = 0; i < N_tot; ++i)
    {
        ix = i % N;
        iy = i / N % N;
        iz = i / (N*N);
        grid[i] = mesh(ix, iy, iz);
    }
}

template<typename T>
void transform_Mesh_to_Grid(const Mesh& mesh, MultiGrid<3, T> &grid)
{
    transform_Mesh_to_Grid(mesh, grid.get_grid());
    grid.restrict_down_all();
}

template<typename T>
void transform_Grid_to_Mesh(Mesh& mesh, const Grid<3, T> &grid)
{/* copy data in MultiGrid 'N*N*N' onto Mesh 'N*N*(N+2)' */
    unsigned int ix, iy, iz;
    const unsigned N_tot = grid.get_Ntot();
    const unsigned N = grid.get_N();

    if (mesh.N != N) throw std::range_error("Mesh of a different size than Grid!");

    #pragma omp parallel for private(ix, iy, iz)
    for (unsigned i = 0; i < N_tot; ++i)
    {
        ix = i % N;
        iy = i / N % N;
        iz = i / (N*N);
        mesh(ix, iy, iz) = grid[i];
    }
}

template<typename T>
void transform_Grid_to_Mesh(Mesh& mesh, const MultiGrid<3, T> &grid){ transform_Grid_to_Mesh(mesh, grid.get_grid()); }

template<typename T>
void transform_Grid_to_Mesh(Mesh& mesh, const MultiGridSolver<3, T> &sol){ transform_Grid_to_Mesh(mesh, sol.get_grid()); }

template<typename T>
T min(const std::vector<T>& data){ return *std::min_element(data.begin(), data.end()); }

FTYPE_t min(const Mesh& data){ return min(data.data); }

template<typename T>
FTYPE_t min(const Grid<3, T> &grid){ return min(grid.get_vec()); }

template<typename T>
FTYPE_t min(const MultiGrid<3, T> &grid){ return min(grid.get_grid()); }

} // namespace <unique> end

template<typename T>
ChiSolver<T>::ChiSolver(unsigned int N, int Nmin, const Sim_Param& sim, bool verbose):
    MultiGridSolver<3, T>(N, Nmin, verbose), n(sim.chi_opt.n), chi_0(2*sim.chi_opt.beta*MPL*sim.chi_opt.phi),
    chi_prefactor_0( // dimensionless prefactor to poisson equation
    #ifndef CHI_A_UNITS // beta*rho_m,0 / Mpl^2 + computing units
    3*sim.chi_opt.beta*sim.cosmo.Omega_m*pow(sim.cosmo.H0 // beta*rho_m,0 / Mpl^2
    * sim.cosmo.h / c_kms // units factor for 'c = 1' and [L] = Mpc / h
    * sim.box_opt.box_size // dimension factor for laplacian
    ,2))
    #else // beta*rho_m,0 / Mpl^2 / (a^3 chi_a) + computing units
    3/2.*sim.cosmo.Omega_m*pow(sim.cosmo.H0 * sim.cosmo.h / c_kms * sim.box_opt.box_size ,2) / sim.chi_opt.phi)
    #endif
{
    if ((n <= 0) || (n >= 1) || (chi_0 <= 0)) throw out_of_range("invalid values of chameleon power-law potential parameters");
}

template<typename T>
void ChiSolver<T>::set_time(T a_, const Cosmo_Param& cosmo)
{
    a = a_; // set new time
    a_3 = pow(a_, 3);
    D = growth_factor(a_, cosmo);

    #ifdef CHI_A_UNITS
    chi_prefactor = chi_prefactor_0*pow(a, -3.*(2.-n)/(1.-n));
    #endif
}

template<typename T>
void ChiSolver<T>::set_initial_guess()
{
    if (!a_3) throw out_of_range("invalid value of scale factor");
    if (!MultiGridSolver<3, T>::get_external_field_size()) throw out_of_range("overdensity not set"); 
    cout << "Setting initial guess for chameleon field...\n";

    T* const f = MultiGridSolver<3, T>::get_y(); // initial guess
    T const* const rho = MultiGridSolver<3,T>::get_external_field(0, 0); // overdensity
    const unsigned N_tot = MultiGridSolver<3, T>::get_Ntot();

    #pragma omp parallel for
    for (unsigned i = 0; i < N_tot; ++i)
    {
        f[i] = chi_min(rho[i]);
    }
}

// The dicretized equation L(phi)
template<typename T>
T  ChiSolver<T>::l_operator(unsigned int level, std::vector<unsigned int>& index_list, bool addsource)
{ 
    const unsigned int i = index_list[0];
    
    // Gridspacing
    const T h = 1.0/T( MultiGridSolver<3, T>::get_N(level) );

    // Solution and pde-source grid at current level
    T const* const chi = MultiGridSolver<3, T>::get_y(level); // solution

    const T rho = MultiGridSolver<3,T>::get_external_field(level, 0)[i];          

    // The right hand side of the PDE 
    if(chi[i] <= 0)
    { // if the solution is overshot, try bulk field
        MultiGridSolver<3, T>::get_y(level)[i] = rho > -1 ? chi_min(rho) : chi_min(0);
    }
    #ifndef CHI_A_UNITS
    T source = ((1+rho)/a_3 - pow(chi_0/chi[i], 1-n)) * chi_prefactor_0;
    #else
    T source = (1 + rho - pow(chi[i], n - 1)) * chi_prefactor;
    #endif

    // Compute the standard kinetic term [D^2 phi] (in 1D this is phi''_i =  phi_{i+1} + phi_{i-1} - 2 phi_{i} )
    T kinetic = -2*chi[i]*3; // chi, '-2*3' is factor in 3D discrete laplacian
    // go through all surrounding points
    for(auto it = index_list.begin() + 1; it < index_list.end(); ++it) kinetic += chi[*it];

    // source term arising rom restricting the equation down to the lower level
    if( level > 0 && addsource ) source += MultiGridSolver<3, T>::get_multigrid_source(level, i);

    // The discretized equation of motion L_{ijk...}(phi) = 0
    return kinetic/(h*h) - source;
}

// Differential of the L operator: dL_{ijk...}/dphi_{ijk...}
template<typename T>
T  ChiSolver<T>::dl_operator(unsigned int level, std::vector<unsigned int>& index_list){
    // solution
    const T chi = MultiGridSolver<3, T>::get_y(level)[ index_list[0] ];

    // Gridspacing
    const T h = 1.0/T( MultiGridSolver<3, T>::get_N(level) );

    // Derivative of kinetic term
    const T dkinetic = -2.0*3;
    
    // Derivative of source
    #ifndef CHI_A_UNITS // when using initial overdensity, current otherwise
    const T dsource = chi_prefactor_0*(1-n)/chi_0*pow(chi_0/chi, 2-n);
    #else
    const T dsource = chi_prefactor*(1-n)*pow(chi, n-2);
    #endif

    return dkinetic/(h*h) - dsource;
}

template<typename T>
bool ChiSolver<T>::check_convergence(){
    // bring variables from MultiGridSolver<3, T>:: namespace here
    const double& _rms_res_i = MultiGridSolver<3, T>::_rms_res_i;
    const double& _rms_res = MultiGridSolver<3, T>::_rms_res;
    const double& _rms_res_old = MultiGridSolver<3, T>::_rms_res_old;
    const double& _verbose = MultiGridSolver<3, T>::_verbose;
    const double& _istep_vcycle = MultiGridSolver<3, T>::_istep_vcycle;
    const double& _eps_converge = MultiGridSolver<3, T>::_eps_converge;
    const double& _maxsteps = MultiGridSolver<3, T>::_maxsteps;

    // Compute ratio of residual to previous residual
    double err = _rms_res_old != 0.0 ? _rms_res/_rms_res_old : 0.0;
    bool converged = false;

    // Print out some information
    if(_verbose){
        std::cout << "    Checking for convergence at step = " << _istep_vcycle << std::endl;
        std::cout << "        Residual = " << _rms_res << "  Residual_old = " <<  _rms_res_old << std::endl;
        std::cout << "        Residual_i = " << _rms_res_i << "  Err = " << err << std::endl;
    }

    // Convergence criterion -- till the residuals are getting smaller
        if(err > 1){ // reject solution, wait for better one
            if(_verbose) std::cout << "    The solution stopped converging err = " << err << " > 1" << std::endl;
            conv_stop = true;
        } else if((err > _eps_converge) && (err <= 1)){
            std::cout << std::endl;
            std::cout << "    The solution stopped converging fast enough err = " << err << " > " << _eps_converge << " ( res = " << _rms_res << " ) istep = " << _istep_vcycle << "\n" << std::endl;
            converged = true;
        } else {
            if (conv_stop){
                std::cout << std::endl;
                std::cout << "    The solution has converged err = " << err << " > " << _eps_converge << " ( res = " << _rms_res << " ) istep = " << _istep_vcycle << "\n" << std::endl;
                converged = true;
            } else {
                if(_verbose) std::cout << "    The solution is still converging err = " << err << " !> " << _eps_converge << std::endl;
            }
        }

    // Define converged if istep exceeds maxsteps to avoid infinite loop...
    if(_istep_vcycle >= _maxsteps){
        std::cout << "    WARNING: MultigridSolver failed to converge! Reached istep = maxsteps = " << _maxsteps << std::endl;
        std::cout << "    res = " << _rms_res << " res_old = " << _rms_res_old << " res_i = " << _rms_res_i << std::endl;
        converged  = true;
    }

    if (converged) conv_stop = false;
    return converged;
}

template<typename T>
T  ChiSolver<T>::chi_min(T delta) const
{
    #ifndef CHI_A_UNITS
    return chi_0*std::pow(a_3/(1+delta), 1/(1-n));
    #else
    return std::pow(1+delta, 1/(n-1));
    #endif
}

App_Var_chi::App_Var_chi(const Sim_Param &sim, std::string app_str):
    App_Var<Particle_v<FTYPE_t>>(sim, app_str), sol(sim.box_opt.mesh_num, sim, false), drho(sim.box_opt.mesh_num)
{
    // EFFICIENTLY ALLOCATE VECTOR OF MESHES
    chi_force.reserve(3);
    for(size_t i = 0; i < 3; i++){
        chi_force.emplace_back(sim.box_opt.mesh_num);
    }
    memory_alloc += sizeof(FTYPE_t)*app_field[0].length*app_field.size();

    // SET CHI SOLVER
    sol.add_external_grid(&drho);
    sol.set_maxsteps(20);
}

static void w_k_correction(Mesh& rho_k)
{
	FTYPE_t w_k;
    Vec_3D<int> k_vec;
    const unsigned NM = rho_k.N;
    const unsigned half_length = rho_k.length / 2;

	#pragma omp parallel for private(w_k, k_vec)
	for(unsigned i=0; i < half_length;i++)
	{
		w_k = 1.;
		get_k_vec(NM, i, k_vec);
		for (int j = 0; j < 3; j++) if (k_vec[j] != 0) w_k *= pow(sin(k_vec[j]*PI/NM)/(k_vec[j]*PI/NM), 2); //< ORDER = 1 (CIC)
        rho_k[2*i] /= w_k;
		rho_k[2*i+1] /= w_k;
	}
}

void App_Var_chi::save_drho_from_particles(Mesh& aux_field)
{
    cout << "Storing density distribution...\n";
    get_rho_from_par(particles, aux_field, sim);
    #ifdef CHI_W_CORR
    fftw_execute_dft_r2c(p_F, aux_field);
    w_k_correction(aux_field);
    fftw_execute_dft_c2r(p_B, aux_field);
    #endif
    transform_Mesh_to_Grid(aux_field, drho);
}

void App_Var_chi::solve(FTYPE_t a)
{
    sol.set_time(a, sim.cosmo);
    sol.set_epsilon(0.98);
    save_drho_from_particles(chi_force[0]);
    cout << "Solving equations of motion for chameleon field...\n";
    sol.solve();
}

void App_Var_chi::print_output()
{/* Print standard output */
    App_Var<Particle_v<FTYPE_t>>::print_output();

    /* Chameleon power spectrum */
    if (sim.out_opt.print_pwr)
    {
        solve(b);

        transform_Grid_to_Mesh(chi_force[0], sol); // get solution
        fftw_execute_dft_r2c(p_F, chi_force[0]); // get chi(k)
        pwr_spec_k_init(chi_force[0], chi_force[0]); // get chi(k)^2, NO w_k correction
        gen_pow_spec_binned(sim, chi_force[0], pwr_spec_binned); // get average Pk
        print_pow_spec(pwr_spec_binned, out_dir_app, "_chi" + z_suffix()); // print
    }
}
template class ChiSolver<CHI_PREC_t>;

#ifdef TEST
#include "test_chameleon.cpp"
#endif