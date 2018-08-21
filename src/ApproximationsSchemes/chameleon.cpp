/**
 * @brief chameleon model of gravity implementation
 * 
 * @file chameleon.cpp
 * @author Michal Vrastil
 * @date 2018-07-08
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

#include <algorithm>

/*****************************//**
 * PRIVATE FUNCTIONS DEFINITIONS *
 *********************************/
namespace{

/**
 * @brief accuracy of chameleon solver
 * @typedef CHI_PREC_t
 * 
 */
typedef double CHI_PREC_t;

/**
 * @brief mass & chi in units of Planck mass
 * @var MPL
 * 
 */
constexpr FTYPE_t MPL = (FTYPE_t)1;

/**
 * @brief speed of light [km / s]
 * @var c_kms
 * 
 */
constexpr FTYPE_t c_kms = (FTYPE_t)299792.458;

/**
 * @brief unphysical value of overdensity (< -1) used to indicate that chameleon filed at this point should not be changed
 * @var MARK_CHI_BOUND_COND
 * 
 */
constexpr CHI_PREC_t MARK_CHI_BOUND_COND =  (CHI_PREC_t)-2;

/**
 * @brief when the relative change in solution is less than this, use Newton`s method. Bisection method otherwise
 * @var SWITCH_BIS_NEW
 * 
 */
constexpr CHI_PREC_t SWITCH_BIS_NEW = (CHI_PREC_t)0.1;

/**
 * @brief minimum value of chameleon field, in chi_a units it is '0', in phi units it is '-1'
 * @var CHI_MIN
 * 
 */
constexpr CHI_PREC_t CHI_MIN = (CHI_PREC_t)-1;

/**
 * \defgroup CHICON Chameleon convergence criteria
 * @{
 */
constexpr double CONVERGENCE_RES = 0; ///< stop when total (rms) residual below
constexpr double CONVERGENCE_RES_MIN = 0; ///< do not stop if solution didn`t converge below
constexpr double CONVERGENCE_ERR = 0.97; ///< stop when improvements between steps slow below
constexpr double CONVERGENCE_ERR_MIN = 0.8; ///< do not stop if solution is still improving
constexpr unsigned CONVERGENCE_NUM_FAIL = 5; ///< stop when number of failed steps is over
constexpr unsigned CONVERGENCE_BI_STEPS = 5; ///< maximal number of steps inside bisection rootfindg method
constexpr unsigned CONVERGENCE_BI_STEPS_INIT = 3; ///< maximal number of steps inside bisection initialization method
constexpr CHI_PREC_t CONVERGENCE_BI_DCHI = (CHI_PREC_t)1e-2; ///< stop bisection method when chi doesn`t chanege
constexpr CHI_PREC_t CONVERGENCE_BI_L = (CHI_PREC_t)1e-2; ///< stop bisection method when residual below
/**@}*/

/**
 * @brief acces \p Exit_Status namespace
 */
using ES = MultiGridSolver<3, CHI_PREC_t>::Exit_Status;

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
void transform_Mesh_to_MultiGrid(const Mesh& mesh, MultiGrid<3, T> &mltgrid)
{
    transform_Mesh_to_Grid(mesh, mltgrid.get_grid());
    mltgrid.restrict_down_all();
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
void transform_MultiGrid_to_Mesh(Mesh& mesh, const MultiGrid<3, T> &mltgrid){ transform_Grid_to_Mesh(mesh, mltgrid.get_grid()); }

template<typename T>
void transform_MultiGridSolver_to_Mesh(Mesh& mesh, const MultiGridSolver<3, T> &sol){ transform_Grid_to_Mesh(mesh, sol.get_grid()); }

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
public:
    // Parameters
    const T n;       ///< Hu-Sawicki paramater
    const T beta;   ///< Chameleon coupling constant
    const T chi_0;   ///< 2*beta*Mpl*phi_scr
    const T phi_prefactor;   ///< prefactor for Poisson equation for gravitational potential
    const T chi_prefactor_0; ///< dimensionless prefactor to poisson equation at a = 1
    T chi_prefactor; ///< time-dependent prefactor

    // convergence parameters
    unsigned m_conv_stop = 0; // number of unsuccessful sweeps
    double m_rms_stop_min;      // iterate at least until _rms_res < m_conv_eps_min
    double m_err_stop;     // stop iteration when: 1 >  err > m_err_stop
    double m_err_stop_min; // iterate at least until: err > m_err_stop_min
    unsigned m_num_fail;     // give up converging if number of failed iteration (err > 1) is > m_num_fail

    // bisection convergence parameters
    unsigned m_max_bisection_steps; // at given point perfom max this number of inteval halving
    T m_dchi_stop;                  // if change in chi is small, stop halving
    T m_l_stop;                     // if residuum is small, stop halving

    // variables for checking solution in deep-screened regime, for each level
    std::vector<std::map<unsigned, T>> fix_vals; ///< <index, value>


    ChiSolver(unsigned int N, int Nmin, const Sim_Param& sim, bool verbose = true) :
        MultiGridSolver<3, T>(N, Nmin, verbose), n(sim.chi_opt.n), beta(sim.chi_opt.beta), chi_0(2*beta*MPL*sim.chi_opt.phi),
        phi_prefactor( // prefactor for Poisson equation for gravitational potential, [] = (h/Mpc)^2
            // 4*pi*G*rho_m,0 + computing units [Mpc/h]
            3/2.*sim.cosmo.Omega_m*pow(sim.cosmo.H0 * sim.cosmo.h / c_kms ,2)
        ),
        chi_prefactor_0( // dimensionless prefactor to poisson equation
            // beta*rho_m,0 / (Mpl*chi_0) + computing units [Mpc/h]; additional a^(-3 -3/(1-n)) at each timestep
            phi_prefactor*pow(sim.box_opt.box_size ,2) / sim.chi_opt.phi
        ),
        fix_vals(this->get_Nlevel())
        {
            if ((n <= 0) || (n >= 1) || (chi_0 <= 0)) throw std::out_of_range("invalid values of chameleon power-law potential parameters");
        }

    ChiSolver(unsigned int N, const Sim_Param& sim, bool verbose = true) : ChiSolver(N, 2, sim, verbose) {}

    T chi_a(T a) const
    {
        return chi_0*pow(a, 3/(1-n));
    }

    T chi_force_units(T a) const
    {
        return beta/MPL*chi_a(a)/phi_prefactor;
    }

    void get_chi_k(Mesh& rho_k)
    {/* transform input density in k-space into linear prediction for chameleon field,
        includes w(k) corrections for interpolation of particles */
        const unsigned N = rho_k.N;
        const unsigned l_half = rho_k.length/2;
        const T mass_sq = (1-n)*chi_prefactor/pow(2*PI, 2); // dimensionless square mass, with derivative factor k* = 2*PI / L
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
    
    bool check_surr_dens(T const* const rho_grid, std::vector<unsigned int> index_list, unsigned i, unsigned N)
    {/* internal method for finding highest density in nearby points */
        // never fix bulk field in under-dense region
        if (rho_grid[i] <= 0) return false;

        // check surrounding points if theres is higher density
        this->get_neighbor_gridindex(index_list, i, N);
        for(unsigned i_s : index_list) if (rho_grid[i_s] > rho_grid[i]) return false;

        // if current point has the highest density, fix chameleon value
        return true;
    }

    void set_time(T a, const Cosmo_Param& cosmo)
    {
        chi_prefactor = chi_prefactor_0*pow(a, -3.*(2.-n)/(1.-n));
    }

    T get_chi_prefactor() const { return chi_prefactor; }
    T get_phi_prefactor() const { return phi_prefactor; }

    T  l_operator(const T chi_i, const unsigned int level, const std::vector<unsigned int>& index_list, const bool addsource, const T h) const 
    {/* The dicretized equation L(phi) */
        // Solution and pde-source grid at current level
        const unsigned int i = index_list[0];
        T const* const chi = this->get_y(level); // solution
        const T rho = this->get_external_field(level, 0)[i];
        
        // do not change values in screened regions
        if (rho == MARK_CHI_BOUND_COND) return 0;

        // The right hand side of the PDE 
        T source = (1 + rho - pow(chi_i - CHI_MIN, n - 1)) * chi_prefactor;

        // Compute the standard kinetic term [D^2 phi] (in 1D this is phi''_i =  phi_{i+1} + phi_{i-1} - 2 phi_{i} )
        T kinetic = -2*chi_i*3; // chi, '-2*3' is factor in 3D discrete laplacian
        // go through all surrounding points
        for(auto it = index_list.begin() + 1; it < index_list.end(); ++it) kinetic += chi[*it];

        // source term arising rom restricting the equation down to the lower level
        if( level > 0 && addsource ) source += this->get_multigrid_source(level, i);

        // The discretized equation of motion L_{ijk...}(phi) = 0
        return kinetic/(h*h) - source;
    }

    T  l_operator(const unsigned int level, const std::vector<unsigned int>& index_list, const bool addsource, const T h) const override
    {/* The dicretized equation L(phi) */
        // Solution and pde-source grid at current level
        const unsigned int i = index_list[0];
        T const* const chi = this->get_y(level); // solution
        const T chi_i = chi[i];
        return l_operator(chi_i, level, index_list, addsource, h);
    }

    // Differential of the L operator: dL_{ijk...}/dphi_{ijk...}
    T dl_operator(const unsigned int level, const std::vector<unsigned int>& index_list, const T h) const override
    {
        // solution
        const T chi = this->get_y(level)[ index_list[0] ];

        // Derivative of source
        const T dsource = chi_prefactor*(1-n)*pow(chi - CHI_MIN, n-2);

        // Derivative of kinetic term
        const T dkinetic = -2.0*3;

        return dkinetic/(h*h) - dsource;
    }

    bool find_opposite_l_sign(const T f1, const T l1, T df, T& f2, T& l2, const unsigned int level, const std::vector<unsigned int>& index_list, const T h) const
    {/* find such 'f2' that 'l_operator(f2)' has opposite sign than l1
        use df as a guess in which direction to be looking */
        f2 = f1;
        for (unsigned j = 0; j < CONVERGENCE_BI_STEPS_INIT; ++j)
        {
            f2 += df;
            if (f2 <= CHI_MIN)
            {
                const unsigned i = index_list[0];
                const T rho = this->get_external_field(level, 0)[i];

                f2 = chi_min(rho);
                j = CONVERGENCE_BI_STEPS_INIT;//< end of loop but first check for an improvement
            }
            l2 = l_operator(f2, level, index_list, true, h);
            if (l1*l2 <= 0) return true;
        }
        return false;
    }

    bool check_bisection_convergence(const T df_new, const T l_new) const
    {
        return ((std::abs(l_new) < m_l_stop) || (std::abs(df_new) < m_dchi_stop));
    }

    T bisection_step(T& f1, T& l1, T& f2, T& l2, const unsigned int level, const std::vector<unsigned int>& index_list, const T h) const
    {/* given 'f1' and 'f2' with different signs of l_operator(f_i) perform one step of bisection:
        f_new = (f1 + f2) / 2
        change whichever l_operator(f_i) has the same sign as l_operator(f_new)
        return value indicates convergence -- 0 (unphysical for chameleon) not, otherwise yes*/

        const T f_new = (f1 + f2) / 2;
        const T l_new = l_operator(f_new, level, index_list, true, h);

        if (check_bisection_convergence(f2 - f_new, l_new)) return f_new;

        if (l1*l_new > 0) {
            f1 = f_new;
            l1 = l_new;
            
        } else {
            f2 = f_new;
            l2 = l_new;
        }

        return CHI_MIN;
    }

    T bisection(T f1, T l1, const T df, const unsigned int level, const std::vector<unsigned int>& index_list, const T h) const
    {/* initialize bisection solver -- find two values with opposite value of l_operator -- and start iterating */
        T f_new, f2, l2;

        // get initial guess -- return when failed
        if (!find_opposite_l_sign(f1, l1, df, f2, l2, level, index_list, h)) return f2;

        // iterate
        for (unsigned i = 0; i < m_max_bisection_steps; ++i){
            f_new = bisection_step(f1, l1, f2, l2, level, index_list, h);
            if (f_new != CHI_MIN)
            {
                return f_new;
            }
        }

        // failed iteration
        return f1;
    }

    
    T upd_operator(const T f, const unsigned int level, const std::vector<unsigned int>& index_list, const T h) const override
    {/* Method for updating solution:
        if df is large, try bisection, otherwise Newton`s method
        try Newton`s method and check for unphysical values */
        T l  =  l_operator(level, index_list, true, h);
        T dl = dl_operator(level, index_list, h);
        T df = -l/dl;
        T df_rel = std::abs(df/(f - CHI_MIN));

        static_assert((SWITCH_BIS_NEW < 1), "Newton`s method cannot be allowed with negative values. Adjust 'SWITCH_BIS_NEW < 1'.");

        return df_rel < SWITCH_BIS_NEW ? f + df : bisection(f, l, df / 2, level, index_list, h);
    }

    void correct_sol(Grid<3,T>& f, const Grid<3,T>& corr, const unsigned int level) override
    {/* Method for correcting solution when going up,
        check for unphysical values */

        const unsigned Ntot  = f.get_Ntot();
        const T * const rho = this->get_external_field(level, 0);

        // bisection variables
        const unsigned int N = this->get_N(level);
        const T h = 1.0/T( N );
        std::vector<unsigned int> index_list;
        
        #pragma omp parallel for private(index_list)
        for(unsigned i = 0; i < Ntot; i++)
        {
            // do not change values in screened regions
            if (rho[i] == MARK_CHI_BOUND_COND) continue;
            else if (std::abs(corr[i]/(f[i] - CHI_MIN)) < SWITCH_BIS_NEW) f[i] += corr[i];
            else{
                this->get_neighbor_gridindex(index_list, i, N);
                T l =  l_operator(level, index_list, true, h);
                f[i] = bisection(f[i], l, corr[i] / 2, level, index_list, h);
            }   
        }
    }

    /**
     * @brief check if solution already converged
     * 
     * @return ES exit status (SUCCESS, FAILURE, SLOW, ITERATE, MAX_STEPS)
     */
    ES check_convergence() override
    {
        /// - bring variables from this-> namespace here
        const double _rms_res_i = this->_rms_res_i;
        const double _rms_res = this->_rms_res;
        const double _rms_res_old = this->_rms_res_old;
        const double _verbose = this->_verbose;
        const double _istep_vcycle = this->_istep_vcycle;
        const double _eps_converge = this->_eps_converge;
        const double _maxsteps = this->_maxsteps;

        /// - Compute ratio of residual to previous residual
        const double err = _rms_res_old != 0.0 ? _rms_res/_rms_res_old : 0.0;
        ES converged = ES::ITERATE;

        /// - convergence criteria
        bool res_below_optimal_tresh = (_rms_res < _eps_converge) || (_eps_converge == 0.0);
        bool res_below_minimal_tresh = (_rms_res < m_rms_stop_min) || (m_rms_stop_min == 0.0);
        bool fast_convergence = err < m_err_stop_min;
        bool slow_convergence = err > m_err_stop;
        bool worse_solution = err > 1;
        bool over_max_steps = _istep_vcycle >= _maxsteps;

        /// - print some information
        auto print_success = [=](unsigned m_conv_stop){
            std::cout << "\n\tSUCCESS: res = " << _rms_res << ", err = " << err << ", num_err = " << m_conv_stop << " (istep = " << _istep_vcycle << ")\n";
            };
        auto print_failure = [=](unsigned m_conv_stop){
            std::cout << "\tFAILURE: res = " << _rms_res << ", err = " << err << ", num_err = " << m_conv_stop << " (istep = " << _istep_vcycle << ")\n";
            };

        auto print_iterate = [=](unsigned m_conv_stop){
            std::cout << "\tITERATE: res = " << _rms_res << ", err = " << err << ", num_err = " << m_conv_stop << " (istep = " << _istep_vcycle << ")\n";
            };

        /// - Print out some information
        if (_verbose){
            std::cout << "    Checking for convergence at step = " << _istep_vcycle << "\n";
            std::cout << "        Residual = " << _rms_res << "  Residual_old = " <<  _rms_res_old << "\n";
            std::cout << "        Residual_i = " << _rms_res_i << "  Err = " << err << "\n";
        }

        /// - check individual conditions
        if (res_below_optimal_tresh && !fast_convergence){
            print_success(m_conv_stop);
            converged = ES::SUCCESS;
        }
        else if (worse_solution)
        {
            if(_verbose) print_failure(m_conv_stop);
            m_conv_stop++;
            converged = ES::ITERATE;
        }
        else if (slow_convergence)
        {
            if (res_below_minimal_tresh)
            {
                print_success(m_conv_stop);
                converged = ES::SLOW;
            }
            else
            {
                m_conv_stop++;
                if (m_conv_stop >= m_num_fail)
                {
                    print_failure(m_conv_stop);
                    converged = ES::FAILURE;
                }
                else
                {
                     if(_verbose) print_iterate(m_conv_stop);
                    converged = ES::ITERATE;
                }
            }
        }
        else
        {
            if (m_conv_stop >= m_num_fail)
            {
                print_failure(m_conv_stop);
                converged = ES::FAILURE;
            }
            else
            {
                if(_verbose) print_iterate(m_conv_stop);
               converged = ES::ITERATE;
            }
        }

        if(over_max_steps)
        {
            print_failure(m_conv_stop);
            converged = ES::MAX_STEPS;
        }
        
        if (converged != ES::ITERATE)
        {
            m_conv_stop = 0;
        }
        return converged;
    }

    void check_solution(unsigned level, Grid<3,T>& chi) override
    {
        for (auto fv : fix_vals[level]) chi[fv.first] = fv.second;
    }

    void set_convergence(double eps, double err_stop, double err_stop_min, double rms_stop_min, unsigned num_fail)
    {
        this->set_epsilon(eps);
        m_err_stop = err_stop;
        m_err_stop_min = err_stop_min;
        m_rms_stop_min = rms_stop_min;
        m_num_fail = num_fail;
    }

    void set_bisection_convergence(unsigned max_bi_step, T dchi_stop, T l_stop)
    {
        m_max_bisection_steps = max_bi_step;
        m_dchi_stop = dchi_stop;
        m_l_stop = l_stop;
    }

    void set_def_convergence()
    {// set convergence -- compensate for units in which we compute laplace operator
        const double err_mod = get_chi_prefactor();
        set_convergence(err_mod*CONVERGENCE_RES, CONVERGENCE_ERR, CONVERGENCE_ERR_MIN, err_mod*CONVERGENCE_RES_MIN, CONVERGENCE_NUM_FAIL);
        set_bisection_convergence(CONVERGENCE_BI_STEPS_INIT, CONVERGENCE_BI_DCHI, CONVERGENCE_BI_L);
    }

    // set chameleon guess to bulk value, need to set time and add rho before call to this function
    void set_bulk_field()
    {
        if (!chi_prefactor) throw std::out_of_range("invalid value of scale factor");
        if (!this->get_external_field_size()) throw std::out_of_range("overdensity not set"); 
        std::cout << "Setting initial guess for chameleon field...\n";

        T* const f = this->get_y(); // initial guess
        T const* const rho = this->get_external_field(0, 0); // overdensity
        const unsigned N_tot = this->get_Ntot();

        #pragma omp parallel for
        for (unsigned i = 0; i < N_tot; ++i)
        {
            f[i] = chi_min(rho[i]);
        }
    }

    void set_linear_sol_at_level(Mesh& rho, const FFTW_PLAN_TYPE& p_F, const FFTW_PLAN_TYPE& p_B, unsigned level)
    {/* set chameleon guess to linear prediction at specific level */
        // check we have the right grids
        if (this->get_N(level) != rho.N) throw std::range_error("Mesh of a different size than Grid at current level!");

        // get delta(k)
        fftw_execute_dft_r2c(p_F, rho);

        // get dchi(k)
        get_chi_k(rho);

        // get dchi(x)
        fftw_execute_dft_c2r(p_B, rho);

        // transform dchi into chi, in chi_a units this means only 'chi = 1 + dchi' */
        rho += 1 + CHI_MIN;

        // copy dchi(x) onto Grid at level
        transform_Mesh_to_Grid(rho, this->get_grid(level));
    }

    void set_linear_recursively(unsigned level)
    {
        // we are at the bottom level
        if (level >= this->get_Nlevel()) return;

        // extract parameters
        const unsigned N = this->get_N(level);

        // create temporary Mesh to compute linear potential, copy density
        Mesh rho(N);
        transform_Grid_to_Mesh(rho, this->get_external_grid(level, 0));

        // initialize FFTW
        const FFTW_PLAN_TYPE p_F = FFTW_PLAN_R2C(N, N, N, rho.real(), rho.complex(), FFTW_ESTIMATE);
        const FFTW_PLAN_TYPE p_B = FFTW_PLAN_C2R(N, N, N, rho.complex(), rho.real(), FFTW_ESTIMATE);

        // set linear prediction
        set_linear_sol_at_level(rho, p_F, p_B, level);

        // solve linear prediction for the next level
        set_linear_recursively(level + 1);
    }

    void set_linear(Mesh& rho, const FFTW_PLAN_TYPE& p_F, const FFTW_PLAN_TYPE& p_B)
    {/* set chameleon guess to linear prediction */

        // solve level = 0, use already allocated space and created plans
        set_linear_sol_at_level(rho, p_F, p_B, 0);

        // recursively solve at level > 0, create new grids and plans (small)
        set_linear_recursively(1);
    }

    void set_screened(unsigned level = 0)
    {/* check solution for invalid values (non-linear regime), fix values in high density regions, try to improve guess in others */
        if (level >= this->get_Nlevel()) return; ///< we are at the bottom level

        Grid<3, T>& chi = this->get_grid(level); // guess
        const unsigned N_tot = this->get_Ntot(level);
        const unsigned N = this->get_N(level);
        T* const rho_grid = this->get_external_field(level, 0); // overdensity
        std::vector<unsigned int> index_list;

        #pragma omp parallel for private(index_list)
        for (unsigned i = 0; i < N_tot; ++i)
        {
            if (chi[i] <= CHI_MIN) // non-linear regime
            {
                if (check_surr_dens(rho_grid, index_list, i, N)) chi[i] = CHI_MIN; ///< CHI_MIN to indicate high-density region
                else chi[i] = rho_grid[i] > 0 ? chi_min(rho_grid[i]) : 1/2. + CHI_MIN;//< phi_s / 2 in underdense region as starting point
            }
        }

        fix_vals[level].clear();

        // writing into map, do not use omp
        for (unsigned i = 0; i < N_tot; ++i)
        {
            if (chi[i] == CHI_MIN) // fix chameleon to bulk value, set unphysical density to indicate screened regime
            {

                chi[i] = chi_min(rho_grid[i]);
                rho_grid[i] = MARK_CHI_BOUND_COND;
                fix_vals[level].emplace(i, chi[i]);
            }
        }

        unsigned num_high_density = fix_vals[level].size();
        std::cout << "Identified and fixed " << num_high_density << "(" << std::setprecision(2) << num_high_density*100.0/N_tot <<  "%) points at level " << level << "\n";

        set_screened(level + 1); ///< recursive call to fix all levels
    }

    T chi_min(T delta) const
    {/* get chi_bulk for given overdensity */
        if (delta > -1) return std::pow(1+delta, 1/(n-1)) + CHI_MIN;
        else if (delta == -1) return 1 + CHI_MIN;
        else throw std::range_error("invalid values of density: " + std::to_string(delta) + " < -1)");        
    }

}; // class ChiSolver end
} ///< end of anonymous namespace (private definitions)

/****************************//**
 * PUBLIC FUNCTIONS DEFINITIONS *
 ********************************/

/**
 * @class:	App_Var_Chi
 * @brief:	class containing variables and methods for chameleon gravity
 */

class App_Var_Chi::ChiImpl
{
public:
    // CONSTRUCTOR
    ChiImpl(const Sim_Param &sim):
        sol(sim.box_opt.mesh_num, sim, false), drho(sim.box_opt.mesh_num), N_level_orig(sol.get_Nlevel()), x_0(sim.x_0())
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
        sol.set_bisection_convergence(CONVERGENCE_BI_STEPS_INIT, CONVERGENCE_BI_DCHI, CONVERGENCE_BI_L);
    }

    // VARIABLES
    ChiSolver<CHI_PREC_t> sol;
    MultiGrid<3, CHI_PREC_t> drho;
    std::vector<Mesh> chi_force;
    uint64_t memory_alloc;

    // METHODS
    void solve(FTYPE_t a, const std::vector<Particle_v<FTYPE_t>>& particles, const Sim_Param &sim, const FFTW_PLAN_TYPE& p_F, const FFTW_PLAN_TYPE& p_B)
    {
        /// - set prefactor
        sol.set_time(a, sim.cosmo);

        /// - set convergence -- compensate for units in which we compute laplace operator
        sol.set_def_convergence();

        /// - save (over)density from particles
        std::cout << "Storing density distribution...\n";
        get_rho_from_par(particles, chi_force[0], sim);
        transform_Mesh_to_MultiGrid(chi_force[0], drho);

        /// - set guess from linear theory and correct unphysical values
        std::cout << "Setting linear guess for chameleon field...\n";
        sol.set_linear(chi_force[0], p_F, p_B);
        sol.set_screened();

        /// - get multigrid_solver runnig
        std::cout << "Solving equations of motion for chameleon field...\n";
        solve_multigrid(); ///< solve using multigrid teqniques
        solve_finest(); ///< solve only on the finest mesh using NGS sweeps
    }

    void gen_pow_spec_binned(const Sim_Param &sim, Data_Vec<FTYPE_t,2>& pwr_spec_binned, const FFTW_PLAN_TYPE& p_F)
    {
        transform_MultiGridSolver_to_Mesh(chi_force[0], sol); // - get solution
        fftw_execute_dft_r2c(p_F, chi_force[0]); // - get chi(k)
        pwr_spec_k_init(chi_force[0], chi_force[0]); // - get chi(k)^2, NO w_k correction
        ::gen_pow_spec_binned(sim, chi_force[0], pwr_spec_binned); // - get average Pk
    }

    void get_chi_force(const FFTW_PLAN_TYPE& p_F, const FFTW_PLAN_TYPE& p_B)
    {
        transform_MultiGridSolver_to_Mesh(chi_force[0], sol); // - get solution
        fftw_execute_dft_r2c(p_F, chi_force[0]); // - get chi(k)
        gen_displ_k_cic(chi_force, chi_force[0]); // - get -k*chi(k)
        fftw_execute_dft_c2r_triple(p_B, chi_force);// - get chi force
    }

    void kick_step_w_chi(const Cosmo_Param &cosmo, const FTYPE_t a, const FTYPE_t da, std::vector<Particle_v<FTYPE_t>>& particles, const std::vector< Mesh> &force_field)
    {
        const unsigned Np = particles.size();
        Vec_3D<FTYPE_t> force;
        const FTYPE_t D = growth_factor(a, cosmo);
        const FTYPE_t OL = cosmo.Omega_L()*pow(a,3);
        const FTYPE_t Om = cosmo.Omega_m;
        const FTYPE_t OLa = OL/(Om+OL);

        /// - -3/2a represents usual EOM, the rest are LCDM corrections
        const FTYPE_t f1 = 3/(2*a)*(1 + OLa);
        const FTYPE_t f2 = 3/(2*a)*(1 - OLa)*D/a;
        /// - chameleon force factor + units
        const FTYPE_t f3 = a/D*sol.chi_force_units(a)/pow2(x_0);

        // force.fill(0.);
        // assign_from(force_field, particles[0].position, force);
        // std::cout << "\n=======> Gr = " << force.norm() << "\tGr_x = " << force[0];
        // force.fill(0.);
        // assign_from(chi_force, particles[0].position, force, f3);
        // std::cout << "\tChi = " << force.norm() << "\tChi_x = " << force[0] << "\n";        
        
        #pragma omp parallel for private(force)
        for (unsigned i = 0; i < Np; i++)
        {
            force.fill(0.);
            assign_from(force_field, particles[i].position, force);
            assign_from(chi_force, particles[i].position, force, f3);
            force = force*f2 - particles[i].velocity*f1;		
            particles[i].velocity += force*da;
        }
    }

private:
    const unsigned N_level_orig;
    const FTYPE_t x_0;

    void solve_multigrid()
    {
        sol.set_ngs_sweeps(3, 6); ///< fine, coarse
        sol.set_maxsteps(30);
        ES status = sol.solve();
    }

    void solve_finest()
    {
        sol.set_ngs_sweeps(5, 0); ///< fine, coarse
        sol.set_Nlevel(1);
        sol.set_maxsteps(30);
        ES status = sol.solve();
        sol.set_Nlevel(N_level_orig);
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
        m_impl->get_chi_force(p_F, p_B);
        //kick_step_w_momentum(sim.cosmo, a_half(), da(), particles, app_field);
        m_impl->kick_step_w_chi(sim.cosmo, a_half(), da(), particles, app_field);
    };
    stream_kick_stream(da(), particles, kick_step, sim.box_opt.mesh_num);
}