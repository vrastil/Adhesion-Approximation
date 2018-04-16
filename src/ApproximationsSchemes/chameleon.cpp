/**
 * @file:	chameleon.cpp
 * @brief:	chameleon solver
 * @note:   chameleon field is in units of 'chi_a = 2*beta*Mpl*Phi_s*a^3/1-n)
 */

#include "ApproximationsSchemes/chameleon.hpp"
#include "core_app.h"
#include "core_power.h"
#include "core_mesh.h"
#include "core_out.h"
#include "integration.hpp"
#include "params.hpp"
#include "MultiGridSolver/multigrid_solver.h"

using namespace std;
namespace{

// accuracy of chameleon solver
typedef long double CHI_PREC_t;

// mass & chi in units of Planck mass
constexpr FTYPE_t MPL = (FTYPE_t)1;

// speed of light [km / s]
constexpr FTYPE_t c_kms = (FTYPE_t)299792.458;

// unphysical value of overdensity (< -1) used to indicate that chameleon filed at this point should not be changed
constexpr CHI_PREC_t MARK_CHI_BOUND_COND =  (CHI_PREC_t)-2;

// correction multiplier for chameleon field when the solution is overshot
constexpr CHI_PREC_t CHI_CORR_MULT = (CHI_PREC_t)0.8;

// slow-down multiplier for Newton`s method
constexpr CHI_PREC_t CHI_SLOW_MULT = (CHI_PREC_t)0.05;

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

/**
 * @class:	ChiSolver
 * @brief:	class to solve chameleon equations of motion
 */

template<typename T>
class ChiSolver : public MultiGridSolver<3, T>
{
private:
    // Parameters
    const T n;       // Hu-Sawicki paramater
    const T chi_0;   // 2*beta*Mpl*phi_scr
    const T chi_prefactor_0; // dimensionless prefactor to poisson equation at a = 1
    T chi_prefactor; // time-dependent prefactor

    // convergence parameters
    unsigned m_conv_stop = 0; // number of unsuccessful sweeps
    double m_rms_stop_min;      // iterate at least until _rms_res > m_conv_eps_min
    double m_err_stop;     // stop iteration when: 1 >  err > m_err_stop
    double m_err_stop_min; // iterate at least until: err > m_err_stop_min
    double m_num_fail;     // give up converging if number of failed iteration (err > 1) is > m_num_fail

    void get_chi_k(Mesh& rho_k, const T h)
    {/* transform input density in k-space into linear prediction for chameleon field,
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

    void get_chi_x()
    {/* transform dchi into chi, in chi_a units this means only 'chi = 1 + dchi' */
        Grid<3, T>& chi = MultiGridSolver<3, T>::get_grid(); // guess
        const unsigned N_tot = MultiGridSolver<3, T>::get_Ntot();

        #pragma omp parallel for
        for (unsigned i = 0; i < N_tot; ++i) ++chi[i];
    }

    
    bool check_surr_dens(T const* const rho_grid, std::vector<unsigned int> index_list, unsigned i, unsigned N)
    {/* internal method for finding highest density in nearby points */
        // never fix bulk field in under-dense region
        if (rho_grid[i] <= 0) return false;

        // check surrounding points if theres is higher density
        MultiGridSolver<3, T>::get_neighbor_gridindex(index_list, i, N);
        for(unsigned i_s : index_list) if (rho_grid[i_s] > rho_grid[i]) return false;

        // if current point has the highest density, fix chameleon value
        return true;
    }

public:
    ChiSolver(unsigned int N, int Nmin, const Sim_Param& sim, bool verbose = true) :
        MultiGridSolver<3, T>(N, Nmin, verbose), n(sim.chi_opt.n), chi_0(2*sim.chi_opt.beta*MPL*sim.chi_opt.phi),
        chi_prefactor_0( // dimensionless prefactor to poisson equation
        // beta*rho_m,0 / (Mpl*chi_0) + computing units; additional a^(-3 -3/(1-n)) at each timestep
        3/2.*sim.cosmo.Omega_m*pow(sim.cosmo.H0 * sim.cosmo.h / c_kms * sim.box_opt.box_size ,2) / sim.chi_opt.phi)
        {
            if ((n <= 0) || (n >= 1) || (chi_0 <= 0)) throw out_of_range("invalid values of chameleon power-law potential parameters");
        }

    ChiSolver(unsigned int N, const Sim_Param& sim, bool verbose = true) : ChiSolver(N, 2, sim, verbose) {}

    void set_time(T a, const Cosmo_Param& cosmo)
    {
        chi_prefactor = chi_prefactor_0*pow(a, -3.*(2.-n)/(1.-n));
    }

    T  l_operator(unsigned int level, std::vector<unsigned int>& index_list, bool addsource, const T h) override
    {/* The dicretized equation L(phi) */
        // Solution and pde-source grid at current level
        const unsigned int i = index_list[0];
        T const* const chi = MultiGridSolver<3, T>::get_y(level); // solution
        const T chi_i = chi[i];
        const T rho = MultiGridSolver<3,T>::get_external_field(level, 0)[i];

        // do not change values in screened regions
        if (rho == MARK_CHI_BOUND_COND) return 0;

        // The right hand side of the PDE 
        T source = (1 + rho - pow(chi[i], n - 1)) * chi_prefactor;

        // Compute the standard kinetic term [D^2 phi] (in 1D this is phi''_i =  phi_{i+1} + phi_{i-1} - 2 phi_{i} )
        T kinetic = -2*chi_i*3; // chi, '-2*3' is factor in 3D discrete laplacian
        // go through all surrounding points
        for(auto it = index_list.begin() + 1; it < index_list.end(); ++it) kinetic += chi[*it];

        // source term arising rom restricting the equation down to the lower level
        if( level > 0 && addsource ) source += MultiGridSolver<3, T>::get_multigrid_source(level, i);

        // The discretized equation of motion L_{ijk...}(phi) = 0
        return kinetic/(h*h) - source;
    }

    // Differential of the L operator: dL_{ijk...}/dphi_{ijk...}
    T dl_operator(unsigned int level, std::vector<unsigned int>& index_list, const T h) override
    {
        // solution
        const T chi = MultiGridSolver<3, T>::get_y(level)[ index_list[0] ];

        // Derivative of source
        const T dsource = chi_prefactor*(1-n)*pow(chi, n-2);

        // Derivative of kinetic term
        const T dkinetic = -2.0*3;

        return dkinetic/(h*h) - dsource;
    }

    
    T upd_operator(T f, T l, T dl) override
    {/* Method for updating solution:
        try Newton`s method and check for unphysical values */
        const T f_new = f - CHI_SLOW_MULT*l/dl;
        return f_new > 0 ? f_new : f*CHI_CORR_MULT;
    }

    // 
    void correct_sol(Grid<3,T>& f, const Grid<3,T>& corr) override
    {/* Method for correcting solution when going up,
        check for unphysical values */
        const unsigned Ntot  = f.get_Ntot();
        T f_new;
        #pragma omp parallel for private(f_new)
        for(unsigned i = 0; i < Ntot; i++)
        {
            f_new = f[i] + CHI_SLOW_MULT*corr[i];
            f[i] = f_new > 0 ? f_new : f[i]*CHI_CORR_MULT;
        }
    }

    bool check_convergence() override
    {/* Criterion for defining convergence */
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
            if (_rms_res > m_rms_stop_min) m_conv_stop++;
            else{
                std::cout << "\n    The solution stopped converging fast enough err = " << err << " > " << m_err_stop
                        << " ( res = " << _rms_res << " ) istep = " << _istep_vcycle << "\n\n";
                converged = true;
            }
        } else {
            if (m_conv_stop >= m_num_fail){
                std::cout << "\n    The solution stopped converging num(err > 1)  = " << m_conv_stop
                        << " res = " << _rms_res << " > " << _eps_converge  << " ) istep = " << _istep_vcycle << "\n\n";
                converged = true;
            } else {
                if(_verbose) std::cout << "    The solution is still converging err = " << err << " !> " << m_err_stop
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

    void set_convergence(double eps, double err_stop, double err_stop_min, double rms_stop_min, double num_fail)
    {
        MultiGridSolver<3, T>::set_epsilon(eps);
        m_err_stop = err_stop;
        m_err_stop_min = err_stop_min;
        m_rms_stop_min = rms_stop_min;
        m_num_fail = num_fail;
    }

    // set chameleon guess to bulk value, need to set time and add rho before call to this function
    void set_bulk_field()
    {
        if (!chi_prefactor) throw out_of_range("invalid value of scale factor");
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

    void set_linear(Mesh& rho, const FFTW_PLAN_TYPE& p_F, const FFTW_PLAN_TYPE& p_B, const T h)
    {/* set chameleon guess to liear prediction */
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
    }

    void set_screened()
    {/* check solution for invalid values (non-linear regime), fix values in high density regions, try to improve guess in others */
        Grid<3, T>& chi = MultiGridSolver<3, T>::get_grid(); // guess
        const unsigned N_tot = MultiGridSolver<3, T>::get_Ntot();
        const unsigned N = MultiGridSolver<3, T>::get_N();
        T* const rho_grid = MultiGridSolver<3,T>::get_external_field(0, 0); // overdensity
        std::vector<unsigned int> index_list;
        unsigned num_high_density = 0;

        #pragma omp parallel for private(index_list)
        for (unsigned i = 0; i < N_tot; ++i)
        {
            if (chi[i] <= 0) // non-linear regime
            {
                // '0' to indicate high-density region
                chi[i] = check_surr_dens(rho_grid, index_list, i, N) ? 0 : chi_min(rho_grid[i]);
            }
        }
        #pragma omp parallel for reduction(+:num_high_density)
        for (unsigned i = 0; i < N_tot; ++i)
        {
            if (chi[i] == 0) // fix chameleon to bulk value, set unphysical density to indicate screened regime
            {
                ++num_high_density;
                chi[i] = chi_min(rho_grid[i]);
                rho_grid[i] = MARK_CHI_BOUND_COND;
            }
        }

        cout << "Identified and fixed " << num_high_density << "(" << std::setprecision(2) << num_high_density*100.0/N_tot <<  "%) points\n";
    }

    T chi_min(T delta) const
    {/* get chi_bulk for given overdensity */
        return delta > -1 ? std::pow(1+delta, 1/(n-1)) : 1;
    }

}; // class ChiSolver end

} // namespace <unique> end

/**
 * @class:	App_Var_Chi
 * @brief:	class containing variables and methods for chameleon gravity
 */

class App_Var_Chi::ChiImpl
{
public:
    // CONSTRUCTOR
    ChiImpl(const Sim_Param &sim):
        sol(sim.box_opt.mesh_num, sim, false), drho(sim.box_opt.mesh_num)
    {
        // EFFICIENTLY ALLOCATE VECTOR OF MESHES
        chi_force.reserve(3);
        for(size_t i = 0; i < 3; i++){
            chi_force.emplace_back(sim.box_opt.mesh_num);
        }

        // ALLOCATED MEMORY
        memory_alloc  = sizeof(FTYPE_t)*chi_force[0].length*chi_force.size();
        memory_alloc += sizeof(CHI_PREC_t)*8*(sol.get_Ntot()-1)/7 // MultiGrid<3, CHI_PREC_t>
                                          *3; // _f, _res, _source
        memory_alloc += sizeof(CHI_PREC_t)*8*(sol.get_Ntot()-1)/7;// MultiGrid<3, CHI_PREC_t> drho

        // SET CHI SOLVER
        sol.add_external_grid(&drho);
        sol.set_convergence(1e-6, 0.95, 0.1, 1e-3, 5);
        sol.set_maxsteps(40);
    }

    // VARIABLES
    ChiSolver<CHI_PREC_t> sol;
    MultiGrid<3, CHI_PREC_t> drho;
    std::vector<Mesh> chi_force;
    uint64_t memory_alloc;

    // METHODS
    void solve(FTYPE_t a, const std::vector<Particle_v<FTYPE_t>>& particles, const Sim_Param &sim, const FFTW_PLAN_TYPE& p_F, const FFTW_PLAN_TYPE& p_B)
    {
        // set prefactor
        sol.set_time(a, sim.cosmo);

        // save (over)density from particles
        cout << "Storing density distribution...\n";
        get_rho_from_par(particles, chi_force[0], sim);
        transform_Mesh_to_Grid(chi_force[0], drho);

        // set guess from linear theory and correct unphysical values
        cout << "Setting linear guess for chameleon field...\n";
        sol.set_linear(chi_force[0], p_F, p_B, 1/(2*PI));
        sol.set_screened();

        // get multigrid_solver runnig
        cout << "Solving equations of motion for chameleon field...\n";
        sol.solve();
    }

    void gen_pow_spec_binned(const Sim_Param &sim, Data_Vec<FTYPE_t,2>& pwr_spec_binned, const FFTW_PLAN_TYPE& p_F)
    {
        transform_Grid_to_Mesh(chi_force[0], sol); // get solution
        fftw_execute_dft_r2c(p_F, chi_force[0]); // get chi(k)
        pwr_spec_k_init(chi_force[0], chi_force[0]); // get chi(k)^2, NO w_k correction
        ::gen_pow_spec_binned(sim, chi_force[0], pwr_spec_binned); // get average Pk
    }
};

App_Var_Chi::App_Var_Chi(const Sim_Param &sim):
    App_Var<Particle_v<FTYPE_t>>(sim, "CHI", "Chameleon gravity"), m_impl(new ChiImpl(sim))
{
    memory_alloc += m_impl->memory_alloc;
}

App_Var_Chi::~App_Var_Chi() = default;

void App_Var_Chi::print_output()
{
    /* Print standard output */
    App_Var<Particle_v<FTYPE_t>>::print_output();

    /* Chameleon power spectrum */
    if (sim.out_opt.print_pwr)
    {
        m_impl->solve(a(), particles, sim, p_F, p_B); // get solution for current time
        m_impl->gen_pow_spec_binned(sim, pwr_spec_binned, p_F); // get chameleon power spectrum
        print_pow_spec(pwr_spec_binned, get_out_dir(), "_chi" + get_z_suffix()); // print
    }
}

void App_Var_Chi::upd_pos()
{// Leapfrog method for chameleon gravity (frozen-potential)
    auto kick_step = [&]()
    {
        m_impl->solve(a_half(), particles, sim, p_F, p_B);
        kick_step_w_momentum(sim.cosmo, a_half(), da(), particles, app_field);
    };
    stream_kick_stream(da(), particles, kick_step, sim.box_opt.mesh_num);
}

#ifdef TEST
#include "test_chameleon.cpp"
#endif