/**
 * @brief modified frozen-potential approximation implementation
 * 
 * @file mod_frozen_potential.cpp
 * @author Michal Vrastil
 * @date 2018-07-11
 */
#include "mod_frozen_potential.hpp"
#include "core_app.h"
#include "core_mesh.h"
#include "integration.hpp"
#include "params.hpp"

namespace {

/**
 * @class LinkedList
 * @brief class handling linked lists
 */

class LinkedList
{
public:
	// CONSTRUCTORS & DESTRUCTOR
	LinkedList(size_t par_num, size_t m, FTYPE_t hc):
	    par_num(par_num), Hc(hc), LL(par_num), HOC(m, m, m) {}
	
	// VARIABLES
	size_t par_num;
	FTYPE_t Hc;
	std::vector<size_t> LL;
	Mesh_base<size_t> HOC;
	
	// METHODS
	void get_linked_list(const std::vector<Particle_v<FTYPE_t>>& particles)
    {
        HOC.assign(-1);
        for (size_t i = 0; i < par_num; i++)
        {
            LL[i] = HOC(particles[i].position/Hc);
            HOC(particles[i].position/Hc) = i;
        }
    }
};

FTYPE_t force_ref(const FTYPE_t r, const FTYPE_t a){
	// Reference force for an S_2-shaped particle
	FTYPE_t z = 2 * r / a;
	if (z > 2) return 1 / (r*r);
	else if (z > 1) return (12 / (z*z) - 224 + 896 * z - 840 * z*z + 224 * pow(z, 3) +
							70 * pow(z, 4) - 48 * pow(z, 5) + 7 * pow(z, 6)) / (35 * a*a);
	else return (224 * z - 224 * pow(z, 3) + 70 * pow(z, 4) + 48 * pow(z, 5) - 21 * pow(z, 7)) / (35 * a*a);
}

FTYPE_t force_tot(const FTYPE_t r, const FTYPE_t e2){
	return 1 / (r*r+e2);
}

void force_short(const Sim_Param &sim, const FTYPE_t D, const LinkedList& linked_list, const  std::vector<Particle_v<FTYPE_t>>& particles,
				 const Vec_3D<FTYPE_t>& position, Vec_3D<FTYPE_t>& force, Interp_obj& fs_interp)
{	// Calculate short range force in position, force is added
    #define FORCE_SHORT_NO_INTER
	size_t p;
	Vec_3D<FTYPE_t> dr_vec;
    FTYPE_t dr2;
    FTYPE_t dr; // <-- #ifdef FORCE_SHORT_NO_INTER
    const FTYPE_t m = pow((FTYPE_t)sim.box_opt.Ng, 3) / D;
    const size_t Nm = sim.box_opt.mesh_num;
    const FTYPE_t rs2 = pow2(sim.app_opt.rs);
    const FTYPE_t e2 = pow2(sim.box_opt.Ng*0.1); // <-- #ifdef FORCE_SHORT_NO_INTER

    IT<3> it(position, sim.app_opt.Hc);
    do{
        p = linked_list.HOC(it.vec);
        while (p != -1){
            dr_vec = get_sgn_distance(particles[p].position, position, Nm);
            dr2 = dr_vec.norm2();
            if ((dr2 < rs2) && (dr2 != 0)) // Short range force is set 0 for separation larger than cutoff radius
            {
                #ifdef FORCE_SHORT_NO_INTER
                dr = sqrt(dr2);
                force += dr_vec*(force_tot(dr, e2) - force_ref(dr, sim.app_opt.a))*m/(dr*4*PI);
                #else
                force += dr_vec*(m/sqrt(dr2)*fs_interp.eval(dr2));
                #endif
            }
            p = linked_list.LL[p];
        }
    } while( it.iter() );
}

void kick_step_w_pp(const Sim_Param &sim, const FTYPE_t a, const FTYPE_t da,  std::vector<Particle_v<FTYPE_t>>& particles, const  std::vector< Mesh> &force_field,
                    LinkedList& linked_list, Interp_obj& fs_interp)
{    // 2nd order ODE with long & short range potential
    const size_t Np = particles.size();
    Vec_3D<FTYPE_t> force;
    const FTYPE_t D = growth_factor(a, sim.cosmo);
    const FTYPE_t OL = sim.cosmo.Omega_L()*pow(a,3);
    const FTYPE_t Om = sim.cosmo.Omega_m;
    // -3/2a represents usual EOM, the rest are LCDM corrections
    const FTYPE_t f1 = 3/(2*a)*(Om+2*OL)/(Om+OL);
    const FTYPE_t f2 = 3/(2*a)*Om/(Om+OL)*D/a;
    
    BOOST_LOG_TRIVIAL(debug) << "Creating linked list...";
	linked_list.get_linked_list(particles);

    BOOST_LOG_TRIVIAL(debug) << "Computing short and long range parts of the potential...";
    #pragma omp parallel for private(force)
    for (size_t i = 0; i < Np; i++)
	{
        force.fill(0.);
        assign_from(force_field, particles[i].position, force); // long-range force
        force_short(sim, D, linked_list, particles, particles[i].position, force, fs_interp); // short range force

        force = force*f2 - particles[i].velocity*f1;		
        particles[i].velocity += force*da;
    }
}

}// end of anonymous namespace

/**
 * @class:	App_Var_FP_mod
 * @brief:	class containing variables and methods for modified Frozen-potential approximation
 */

class App_Var_FP_mod::FP_ppImpl
{
public:
    FP_ppImpl(const Sim_Param &sim): linked_list(sim.box_opt.par_num, sim.app_opt.M, sim.app_opt.Hc)
    {
        memory_alloc = sizeof(size_t)*(linked_list.HOC.length+linked_list.par_num);

        // precompute short range force
        size_t res = size_t(sim.app_opt.rs/0.05)+1; // force resolution 5% of mesh cell
        const FTYPE_t r0 = sim.app_opt.rs / (res-1);
        Data_Vec<FTYPE_t, 2> data(res);
        FTYPE_t r;
        const FTYPE_t e2 = pow2(sim.box_opt.Ng*0.1); // softening of 10% of average interparticle length

        #pragma omp parallel for private(r)
        for(size_t i = 0; i < res; i++)
        {
            r = i*r0;
            data[0][i] = pow2(r); // store square of r
            data[1][i] = (force_tot(r, e2) - force_ref(r, sim.app_opt.a))/(4*PI);
        }
        fs_interp.init(data);
    }

	// VARIABLES
    LinkedList linked_list;
    Interp_obj fs_interp;
    uint64_t memory_alloc;
};

 
App_Var_FP_mod::App_Var_FP_mod(const Sim_Param &sim):
    App_Var<Particle_v<FTYPE_t>>(sim, "FP_pp", "Modified Frozen-potential approximation"), m_impl(new FP_ppImpl(sim))
{
    memory_alloc += m_impl->memory_alloc;
}

App_Var_FP_mod::~App_Var_FP_mod() = default;

void App_Var_FP_mod::pot_corr(std::vector<Mesh>& vel_field, Mesh& pot_k)
{
    /* Computing displacement in k-space with S2 shaped particles */
	gen_displ_k_S2(vel_field, pot_k, sim.app_opt.a);
    
    /* Computing force in q-space */
    BOOST_LOG_TRIVIAL(debug) << "Computing force in q-space...";
    fftw_execute_dft_c2r_triple(p_B, vel_field);
}

void App_Var_FP_mod::upd_pos()
{// Leapfrog method for modified frozen-potential
    auto kick_step = [&](){ kick_step_w_pp(sim, a_half(), da(), particles, app_field, m_impl->linked_list, m_impl->fs_interp); };
    stream_kick_stream(da(), particles, kick_step, sim.box_opt.mesh_num);
}