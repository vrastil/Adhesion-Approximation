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
    const T n;        // Hu-Sawicki paramater
    const T chi_0;    // 2*beta*Mpl*phi_scr
    const T chi_prefactor; // beta*rho_m,0 / Mpl, [1]
    T a;    // scale factor
    T a_3;  // prec-compute a^3
    T D;    // pre-compute D(a)

public:
    // Constructors
    ChiSolver() {}
    ChiSolver(unsigned int N, const Sim_Param& sim, bool verbose = true) : ChiSolver(N, 2, sim, verbose) {}
    ChiSolver(unsigned int N, int Nmin, const Sim_Param& sim, bool verbose = true);
    void set_time(T a, const Cosmo_Param& cosmo);

    // The dicretized equation L(phi)
    T  l_operator(unsigned int level, std::vector<unsigned int>& index_list, bool addsource);

    // Differential of the L operator: dL_{ijk...}/dphi_{ijk...}
    T dl_operator(unsigned int level, std::vector<unsigned int>& index_list);

    // set initial guess to bulk value, need to set time and add rho before call to this function
    void set_initial_guess();

    // get chi_bulk for given initial overdensity
    T chi_min(T delta_0) const { return chi_0*std::pow(a_3/(1+D*delta_0), 1/(1-n)); }
};

/**
 * @class:	App_Var_chi
 * @brief:	class containing variables for chameleon gravity
 */
 
 class App_Var_chi: public App_Var<Particle_v<FTYPE>>
{
public:
	// CONSTRUCTORS & DESTRUCTOR
	App_Var_chi(const Sim_Param &sim, std::string app_str);

	// VARIABLES
    ChiSolver<FTYPE> sol;
    MultiGrid<3, FTYPE> drho;

    // METHODS
    void save_init_drho_k(const Mesh& dro, Mesh& aux_field);
};