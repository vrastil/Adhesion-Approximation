#pragma once

#include "stdafx.h"
#include "params.hpp"
#include "app_var.hpp"
#include "MultiGridSolver/multigrid_solver.h"

template<typename T>
class ChiSolver : public MultiGridSolver<3, T>
{
private:
    // Parameters
    T a;    // scale factor
    T a_3;  // prec-compute a^3
    T D;    // pre-compute D(a)

    const T n;        // Hu-Sawicki paramater
    const T chi_0;    // 2*beta*Mpl*phi_scr
    const T chi_prefactor_0; // dimensionless prefactor to poisson equation at a = 1
    T chi_prefactor;

    // convergence parameters
    unsigned m_conv_stop = 0; // number of unsuccessful sweeps
    double m_err_stop;        // stop iteration when: 1 >  err > m_err_stop
    double m_err_stop_min;    // iterate at least until: err > m_err_stop_min

public:
    // Constructors
    ChiSolver(unsigned int N, const Sim_Param& sim, bool verbose = true) : ChiSolver(N, 2, sim, verbose) {}
    ChiSolver(unsigned int N, int Nmin, const Sim_Param& sim, bool verbose = true);
    void set_time(T a, const Cosmo_Param& cosmo);

    // The dicretized equation L(phi)
    T  l_operator(unsigned int level, std::vector<unsigned int>& index_list, bool addsource, const T h) override;

    // Differential of the L operator: dL_{ijk...}/dphi_{ijk...}
    T dl_operator(unsigned int level, std::vector<unsigned int>& index_list, const T h) override;

    // Criterion for defining convergence
    bool check_convergence() override;
    void set_convergence(double eps, double err_stop, double err_stop_min);

    // set chameleon guess to bulk value, need to set time and add rho before call to this function
    void set_bulk_field();

    // set chameleon guess to liear prediction
    void set_linear(Mesh& rho, const FFTW_PLAN_TYPE& p_F, const FFTW_PLAN_TYPE& p_B, const T x_0);

    // get chi_bulk for given overdensity
    T chi_min(T delta) const;
};

/**
 * @class:	App_Var_chi
 * @brief:	class containing variables for chameleon gravity
 */

typedef long double CHI_PREC_t;
 
class App_Var_chi: public App_Var<Particle_v<FTYPE_t>>
{
public:
	// CONSTRUCTORS & DESTRUCTOR
	App_Var_chi(const Sim_Param &sim, std::string app_str);

	// VARIABLES
    ChiSolver<CHI_PREC_t> sol;
    MultiGrid<3, CHI_PREC_t> drho;
    std::vector<Mesh> chi_force;

    // METHODS
    void print_output();
    void solve(FTYPE_t a);

private:
    void save_drho_from_particles(Mesh& aux_field);
};