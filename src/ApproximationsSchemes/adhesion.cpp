/**
 * @brief ahesion approximation implementation
 * 
 * @file adhesion.cpp
 * @author Michal Vrastil
 * @date 2018-07-11
 */

#include "ApproximationsSchemes/adhesion.hpp"
#include "core_app.h"
#include "core_mesh.h"
#include "integration.hpp"
#include "params.hpp"

#include <algorithm>

namespace {
const FTYPE_t ACC = 1e-10;
const FTYPE_t log_acc = log(ACC);

void gen_init_expot(const Mesh& potential, Mesh& expotential, FTYPE_t nu)
{
	printf("Storing initial expotenital in q-space...\n");
    // store exponent only
    // *expotential = potential; !!! <- do not use this in case potential and expotential are meshes of different size
    #pragma omp parallel for
    for (size_t i = 0; i < expotential.length; i++) expotential[i] = -potential[i] / (2*nu);
}

FTYPE_t get_summation(const std::vector<FTYPE_t>& exp_aux)
{
    FTYPE_t max_exp = *max_element(exp_aux.begin(), exp_aux.end());
    FTYPE_t sum = 0;
    for(auto const& a_exp: exp_aux) {
        if ((a_exp - max_exp) > log_acc) sum+= exp(a_exp - max_exp);
    }
    return max_exp + log(sum);
}

void convolution_y1(Mesh& potential, const std::vector<FTYPE_t>& gaussian, const Mesh& expotential_0){
	// multi-thread index is y3
    // compute f1 (x1, y2, y3)

    const int N = potential.N;
    std::vector<FTYPE_t> exp_aux;
    
	#pragma omp parallel for private(exp_aux)
	for (int x1 = 0; x1 < N; x1++){
		for (int y2 = 0; y2 < N; y2++){
			for (int y3 = 0; y3 < N; y3++){
                exp_aux.reserve(N);
                // fill in exponents
                for (int y1 = 0; y1 < N; y1++){
                    exp_aux.push_back(expotential_0(y1, y2, y3)+gaussian[std::abs(x1-y1)]);
				}
				potential(x1, y2, y3) = get_summation(exp_aux); // potential is now f1
                exp_aux.clear();
			}
		}
	}
}

void convolution_y2(Mesh& potential, const std::vector<FTYPE_t>& gaussian){
    // compute f2 (x1, x2, y3)

    const int N = potential.N;
	std::vector<FTYPE_t> sum_aux;
	std::vector<FTYPE_t> exp_aux;

	#pragma omp parallel for private(sum_aux, exp_aux)
	for (int x1 = 0; x1 < N; x1++){
		for (int y3 = 0; y3 < N; y3++){
            sum_aux.reserve(N);
			for (int x2 = 0; x2 < N; x2++){
                exp_aux.reserve(N);
				// fill in exponents
                for (int y2 = 0; y2 < N; y2++){
                    exp_aux.push_back(potential(x1, y2, y3) + gaussian[abs(x2-y2)]);
				}
				sum_aux.push_back(get_summation(exp_aux));
                exp_aux.clear();
			}

			for (int x2 = 0; x2 < N; x2++){
				potential(x1, x2, y3) = sum_aux[x2]; // potential is now f2
			}
			sum_aux.clear();
		}
	}
}

void convolution_y3(Mesh& potential, const std::vector<FTYPE_t>& gaussian){
    // compute f3 (x1, x2, x3) == expotential(x, b)

    const int N = potential.N;
	std::vector<FTYPE_t> sum_aux;
    std::vector<FTYPE_t> exp_aux;

	#pragma omp parallel for private(sum_aux, exp_aux)
	for (int x1 = 0; x1 < N; x1++){
		for (int x2 = 0; x2 < N; x2++){
            sum_aux.reserve(N);
			for (int x3 = 0; x3 < N; x3++){
                exp_aux.reserve(N);
				// fill in exponents
                for (int y3 = 0; y3 < N; y3++){
                    exp_aux.push_back(potential(x1, x2, y3) + gaussian[abs(x3-y3)]);
				}
				sum_aux.push_back(get_summation(exp_aux));
                exp_aux.clear();
			}
			for (int x3 = 0; x3 < N; x3++){
				potential(x1, x2, x3) = sum_aux[x3]; // potential is now f3
			}
			sum_aux.clear();
		}
	}
}

void gen_expot(Mesh& potential,  const Mesh& expotential_0, FTYPE_t nu, FTYPE_t b)
{
	/* Computing convolution using direct sum */
	printf("Computing expotential in q-space...\n");
	/*
	f(x1, x2, x3) = \int dy^3 { g(y1, y2, y3) * h(x1 - y1) * h(x2 - y2) * h(x3 - y3)}
	..
	f1 (x1, y2, y3) = \int dy1 { g  (y1, y2, y3) * h(x1 - y1)}	:: N^3 sums of length N
	f2 (x1, x2, y3) = \int dy2 { f1 (x1, y2, y3) * h(x2 - y2)}	:: N^3 sums of length N
	f3 (x1, x2, x3) = \int dy3 { f2 (x1, x2, y3) * h(x3 - y3)}	:: N^3 sums of length N
	*/

	// store values of exponential - every convolution uses the same exp(-r^2/4bv)
	std::vector<FTYPE_t> gaussian(expotential_0.N);

	#pragma omp parallel for
	for (size_t i = 0; i < expotential_0.N; i++){
		gaussian[i]=-i*i/(4*b*nu);
	}

	convolution_y1(potential, gaussian, expotential_0);
	convolution_y2(potential, gaussian);
	convolution_y3(potential, gaussian);
}

}// end of anonymous namespace

/**
 * @class:	App_Var_AA
 * @brief:	class containing variables and methods for adhesion approximation
 */

class App_Var_AA::AAImpl
{
public:
    AAImpl(const Sim_Param &sim): expotential(sim.box_opt.mesh_num)
    {
        memory_alloc = sizeof(FTYPE_t)*expotential.length;
    }

    void aa_convolution(App_Var_AA& APP)
    {
        printf("Computing potential...\n");	
        gen_expot(APP.app_field[0], expotential, APP.sim.app_opt.nu, APP.a_half());
        APP.app_field[0] *= -2*APP.sim.app_opt.nu;
                    
        printf("Computing velocity field via FFT...\n");
        fftw_execute_dft_r2c(APP.p_F, APP.app_field[0]);
        gen_displ_k(APP.app_field, APP.app_field[0]);
        fftw_execute_dft_c2r_triple(APP.p_B, APP.app_field);
    }

	// VARIABLES
	Mesh expotential;
    uint64_t memory_alloc;
};
 
App_Var_AA::App_Var_AA(const Sim_Param &sim):
    App_Var<Particle_v<FTYPE_t>>(sim, "AA", "Adhesion approximation"), m_impl(new AAImpl(sim))
{
    memory_alloc += m_impl->memory_alloc;
}

App_Var_AA::~App_Var_AA() = default;

void App_Var_AA::pot_corr()
{
    /* Computing initial expotential */
	fftw_execute_dft_c2r(p_B, power_aux[0]);
    gen_init_expot(power_aux[0], m_impl->expotential, sim.app_opt.nu);
}

void App_Var_AA::upd_pos()
{// Leapfrog method for adhesion
    m_impl->aa_convolution(*this);
    auto kick_step = [&](){ kick_step_w_momentum(sim.cosmo, a_half(), da(), particles, app_field); };
    stream_kick_stream(da(), particles, kick_step, sim.box_opt.mesh_num);
}