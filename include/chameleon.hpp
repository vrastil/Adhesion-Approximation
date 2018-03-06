#pragma once

#include "stdafx.h"
#include "params.hpp"
#include "app_var.hpp"
#include "MultiGridSolver/multigrid_solver.h"

// use rho = D*rho_0 ifdef LINEAR_CHI_SOLVER, assign particles onto Mesh at each timestep otherwise
#define LINEAR_CHI_SOLVER

template<typename T>
class ChiSolver : public MultiGridSolver<3, T>
{
private:
    // Parameters
    T a;    // scale factor
    T a_0;  // scale factor from previous step
    T a_3;  // prec-compute a^3
    T D;    // pre-compute D(a)

public:
    // Constructors
    ChiSolver(unsigned int N, const Sim_Param& sim, bool verbose = true) : ChiSolver(N, 2, sim, verbose) {}
    ChiSolver(unsigned int N, int Nmin, const Sim_Param& sim, bool verbose = true);
    void set_time(T a, const Cosmo_Param& cosmo);

    // Parameters
    const T n;        // Hu-Sawicki paramater
    const T chi_0;    // 2*beta*Mpl*phi_scr
    const T chi_prefactor; // beta*rho_m,0 / Mpl, [dimensionless]

    // The dicretized equation L(phi)
    T  l_operator(unsigned int level, std::vector<unsigned int>& index_list, bool addsource);

    // Differential of the L operator: dL_{ijk...}/dphi_{ijk...}
    T dl_operator(unsigned int level, std::vector<unsigned int>& index_list);

    // set initial guess to bulk value, need to set time and add rho before call to this function
    void set_initial_guess();

    // multiply guess from previus step by factor corresponding to evolution of bulk field
    void set_next_guess(const Cosmo_Param& cosmo);

    // get chi_bulk for given overdensity
    T chi_min(T delta) const { return chi_0*std::pow(a_3/(1+delta), 1/(1-n)); }
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
    void save_init_drho_k(const Mesh& dro, Mesh& aux_field);
    void save_drho_from_particles(Mesh& aux_field);
    void print_output();
};