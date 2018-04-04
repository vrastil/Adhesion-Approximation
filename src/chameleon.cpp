/**
 * @file:	chameleon.cpp
 * @brief:	chameleon solver
 * @note:   chameleon field is in units of 'chi_a = 2*beta*Mpl*Phi_s*a^3/1-n)
 */


#include "core_app.h"
#include "core_power.h"
#include "core_mesh.h"
#include "core_out.h"
#include "chameleon.hpp"

using namespace std;
namespace{
constexpr FTYPE_t MPL = (FTYPE_t)1; // mass & chi in units of Planck mass
constexpr FTYPE_t c_kms = (FTYPE_t)299792.458; // speed of light [km / s]

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
    // beta*rho_m,0 / (Mpl*chi_0) + computing units; additional a^(-3 -3/(1-n)) at each timestep
    3/2.*sim.cosmo.Omega_m*pow(sim.cosmo.H0 * sim.cosmo.h / c_kms * sim.box_opt.box_size ,2) / sim.chi_opt.phi)
{
    if ((n <= 0) || (n >= 1) || (chi_0 <= 0)) throw out_of_range("invalid values of chameleon power-law potential parameters");
}

template<typename T>
void ChiSolver<T>::set_time(T a_, const Cosmo_Param& cosmo)
{
    a = a_; // set new time
    a_3 = pow(a_, 3);
    D = growth_factor(a_, cosmo);
    chi_prefactor = chi_prefactor_0*pow(a, -3.*(2.-n)/(1.-n));
}

template<typename T>
void ChiSolver<T>::set_bulk_field()
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

template<typename T>
void ChiSolver<T>::get_chi_k(Mesh& rho_k, const T h)
{/* transform input density in k-space into linear prediction for chameleon field
    includes w(k) corrections for interpolation of particles */

    const unsigned N = rho_k.N;
    const unsigned l_half = rho_k.length/2;
    const T mass_sq = (1-n)*chi_prefactor*pow(h, 2); // dimensionless square mass, with derivative factor
    const T chi_a_n = -1/(1-n); // prefactor for chi(k), in chi_a units
    
    T k2, g_k;

	#pragma omp parallel for private(k2, g_k)
	for(unsigned i=0; i < l_half;i++){
		k2 = get_k_sq(N, i);
		if (k2 == 0)
        {
            rho_k[2*i] = 0;
            rho_k[2*i+1] = 0;
        }
		else
        {
            g_k = chi_a_n/(k2+mass_sq)*mass_sq; // Green function
			rho_k[2*i] *= g_k;
			rho_k[2*i+1] *= g_k;
		}
	}
}

template<typename T>
void ChiSolver<T>::get_chi_x()
{// transform dchi into chi, in chi_a units this means only 'chi = 1 + dchi'
    Grid<3, T>& chi = MultiGridSolver<3, T>::get_grid(); // guess
    const unsigned N_tot = MultiGridSolver<3, T>::get_Ntot();

    #pragma omp parallel for
    for (unsigned i = 0; i < N_tot; ++i) ++chi[i];
}

template<typename T>
void ChiSolver<T>::screen_corr()
{// check for invalid values and/or screening regime, try to improve guess
    return;

    // Grid<3, T>& chi = MultiGridSolver<3, T>::get_grid(); // guess
    // const unsigned N_tot = MultiGridSolver<3, T>::get_Ntot();
    // T const* const rho_grid = MultiGridSolver<3,T>::get_external_field(0, 0); // overdensity
    // std::vector<unsigned int> index_list;
    // bool converged = true;

    // // check for zero values (those in screened regime but outside dense objects)
    // unsigned num_loop = 0;
    // while (!converged){
    //     num_loop++;
    //     converged = true;
    //     #pragma omp parallel for private(index_list)
    //     for (unsigned i = 0; i < N_tot; ++i)
    //     {
    //         if (chi[i] == 0) 
    //         { // assign average value
    //             MultiGridSolver<3, T>::get_neighbor_gridindex(index_list, i, N);
    //             for(auto it = index_list.begin() + 1; it < index_list.end(); ++it) chi[i] += chi[*it];
    //             chi[i] /= 6;
    //             if (chi[i] == 0) converged = false; // wait until next loop
    //         } 
    //     }
    // }
    // cout << "No zero values after " << num_loop << " loops.\n";
}

template<typename T>
void ChiSolver<T>::set_linear(Mesh& rho, const FFTW_PLAN_TYPE& p_F, const FFTW_PLAN_TYPE& p_B, const T h)
{
    // get delta(k)
    fftw_execute_dft_r2c(p_F, rho);

    // get dchi(k)
    get_chi_k(rho, h);

    // get dchi(x)
    fftw_execute_dft_c2r(p_B, rho);

    // copy dchi(x) onto Grid
    transform_Mesh_to_Grid(rho, MultiGridSolver<3, T>::get_grid());

    // get chi(x)
    get_chi_x();

    // check for invalid values
    screen_corr();
}

// The dicretized equation L(phi)
template<typename T>
T  ChiSolver<T>::l_operator(unsigned int level, std::vector<unsigned int>& index_list, bool addsource, const T h)
{ 
    const unsigned int i = index_list[0];

    // Solution and pde-source grid at current level
    T const* const chi = MultiGridSolver<3, T>::get_y(level); // solution

    const T rho = MultiGridSolver<3,T>::get_external_field(level, 0)[i];          

    // The right hand side of the PDE 
    if(chi[i] <= 0)
    { // if the solution is overshot, try bulk field
        MultiGridSolver<3, T>::get_y(level)[i] = chi_min(rho);
        return 0;
    }
    T source = (1 + rho - pow(chi[i], n - 1)) * chi_prefactor;

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
T  ChiSolver<T>::dl_operator(unsigned int level, std::vector<unsigned int>& index_list, const T h){
    // solution
    const T chi = MultiGridSolver<3, T>::get_y(level)[ index_list[0] ];

    // Derivative of kinetic term
    const T dkinetic = -2.0*3;
    
    // Derivative of source
    const T dsource = chi_prefactor*(1-n)*pow(chi, n-2);

    return dkinetic/(h*h) - dsource;
}

template<typename T>
bool ChiSolver<T>::check_convergence(){
    // bring variables from MultiGridSolver<3, T>:: namespace here
    const double _rms_res_i = MultiGridSolver<3, T>::_rms_res_i;
    const double _rms_res = MultiGridSolver<3, T>::_rms_res;
    const double _rms_res_old = MultiGridSolver<3, T>::_rms_res_old;
    const double _verbose = MultiGridSolver<3, T>::_verbose;
    const double _istep_vcycle = MultiGridSolver<3, T>::_istep_vcycle;
    const double _eps_converge = MultiGridSolver<3, T>::_eps_converge;
    const double _maxsteps = MultiGridSolver<3, T>::_maxsteps;

    // Compute ratio of residual to previous residual
    double err = _rms_res_old != 0.0 ? _rms_res/_rms_res_old : 0.0;
    bool converged = false;

    // Print out some information
    if(_verbose){
        std::cout << "    Checking for convergence at step = " << _istep_vcycle << "\n";
        std::cout << "        Residual = " << _rms_res << "  Residual_old = " <<  _rms_res_old << "\n";
        std::cout << "        Residual_i = " << _rms_res_i << "  Err = " << err << "\n";
    }

    // Convergence criterion
    if((_rms_res < _eps_converge) && (err > m_err_stop_min)){ // residuals below treshold
        std::cout << "\n    The solution has converged res = " << _rms_res << " < " << _eps_converge 
                  << " ( err = " << err << " ) istep = " << _istep_vcycle << "\n\n";
        converged = true;
    } else if(err > 1){ // reject solution, wait for better one
        if(_verbose) std::cout << "    The solution stopped converging res = " << _rms_res << " !< " << _eps_converge
                               << " err = " << err << " > 1" << " num(err > 1) = " << m_conv_stop << "\n";
        m_conv_stop++;
    } else if(err > m_err_stop){ // convergence too slow
        std::cout << "\n    The solution stopped converging fast enough err = " << err << " > " << m_err_stop
                  << " ( res = " << _rms_res << " ) istep = " << _istep_vcycle << "\n\n";
        converged = true;
    } else {
        if (m_conv_stop){
            std::cout << "\n    The solution stopped converging num(err > 1)  = " << m_conv_stop
                      << " res = " << _rms_res << " > " << _eps_converge  << " ) istep = " << _istep_vcycle << "\n\n";
            converged = true;
        } else {
            if(_verbose) std::cout << "    The solution is still converging err = " << err << " !> " << _eps_converge
                                   << " res = " << _rms_res << " > " << _eps_converge << "\n";
        }
    }

    // Define converged if istep exceeds maxsteps to avoid infinite loop...
    if(_istep_vcycle >= _maxsteps){
        std::cout << "    WARNING: MultigridSolver failed to converge! Reached istep = maxsteps = " << _maxsteps << "\n";
        std::cout << "    res = " << _rms_res << " res_old = " << _rms_res_old << " res_i = " << _rms_res_i << "\n";
        converged  = true;
    }

    if (converged) m_conv_stop = 0;
    return converged;
}

template<typename T>
void ChiSolver<T>::set_convergence(double eps, double err_stop, double err_stop_min)
{
    MultiGridSolver<3, T>::set_epsilon(eps);
    m_err_stop = err_stop;
    m_err_stop_min = err_stop_min;
}

template<typename T>
T  ChiSolver<T>::chi_min(T delta) const
{
    return delta > -1 ? std::pow(1+delta, 1/(n-1)) : 1;
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
    sol.set_convergence(1e-6, 0.95, 0.1);
    sol.set_maxsteps(15);
}

void App_Var_chi::save_drho_from_particles(Mesh& aux_field)
{
    cout << "Storing density distribution...\n";
    get_rho_from_par(particles, aux_field, sim);
    transform_Mesh_to_Grid(aux_field, drho);
}

void App_Var_chi::solve(FTYPE_t a)
{
    sol.set_time(a, sim.cosmo);
    save_drho_from_particles(chi_force[0]);
    cout << "Setting guess for chameleon field...\n";
    sol.set_linear(chi_force[0], p_F, p_B, 1/(2*PI));
    // cout << "Solving equations of motion for chameleon field...\n";
    // sol.solve();
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