#include "stdafx.h"
#include "core.h"
#include "core_power.h"
#include "chameleon.hpp"

using namespace std;

template<typename T>
inline T chi_min(T chi_0, T a, T n, T delta){ return chi_0*pow(pow(a, 3)/(1+delta), 1/(1-n)); }

template<typename T>
void transform_Mesh_to_Grid(const Mesh& mesh, MultiGrid<3, T> &grid)
{/* copy data in Mesh 'N*N*(N+2)' onto MultiGrid 'N*N*N' */
    unsigned int ix, iy, iz;
    const unsigned N_tot = grid.get_N_tot(0);
    const unsigned N = grid.get_N(0);

    #pragma omp parallel for private(ix, iy, iz)
    for (unsigned i = 0; i < N_tot; ++i)
    {
        ix = i % N;
        iy = i / N % N;
        iz = i / (N*N) % N;
        grid[i] = mesh(ix, iy, iz);
    }
    grid.restrict_down_all();
}

template<typename T>
void transform_Grid_to_Mesh(Mesh& mesh, const MultiGrid<3, T> &grid)
{/* copy data in MultiGrid 'N*N*N' onto Mesh 'N*N*(N+2)' */
    unsigned int ix, iy, iz;
    const unsigned N_tot = grid.get_N_tot(0);
    const unsigned N = grid.get_N(0);

    #pragma omp parallel for private(ix, iy, iz)
    for (unsigned i = 0; i < N_tot; ++i)
    {
        ix = i % N;
        iy = i / N % N;
        iz = i / (N*N) % N;
        mesh(ix, iy, iz) = grid[i];
    }
}

template<typename T>
ChiSolver<T>::ChiSolver(unsigned int N, int Nmin, const Sim_Param& sim, bool verbose):
    MultiGridSolver<3, T>(N, Nmin, verbose), n(sim.chi_opt.n), chi_0(sim.chi_opt.chi_0),
    chi_prefactor( // beta*rho_m,0 / Mpl, [(h/Mpc)^3]
        3*sim.chi_opt.beta*sim.cosmo.H0*sim.cosmo.H0*sim.cosmo.Omega_m/(8*PI*G_N*MPL)
        /(c_kms*c_kms*pow(sim.cosmo.h, 3))) // units factor
{       
        
}

template<typename T>
void ChiSolver<T>::set_time(T a, const Cosmo_Param& cosmo)
{
    a_3 = pow(a, 3);
    D = growth_factor(a, cosmo);
}

// The dicretized equation L(phi)
template<typename T>
T  ChiSolver<T>::l_operator(unsigned int level, std::vector<unsigned int>& index_list, bool addsource)
{ 
    const unsigned int i = index_list[0];
    
    // Gridspacing
    const T h = 1.0/T( MultiGridSolver<3, T>::get_N(level) );

    // Solution and pde-source grid at current level
    T *const dchi = MultiGridSolver<3, T>::get_y(level); // chi = dchi + chi_A
    T *const chi_A = MultiGridSolver<3,T>::get_external_field(level, 0); // solution from previus time-step

    const T rho_0 = MultiGridSolver<3,T>::get_external_field(level, 1)[i]; // initial overdensity
    const T chi = dchi[i] + chi_A[i]; // full solution

    // Compute the standard kinetic term [D^2 phi] (in 1D this is phi''_i =  phi_{i+1} + phi_{i-1} - 2 phi_{i} )
    T kinetic = -2*chi*3; // chi, '-2*3' is factor in 3D discrete laplacian
    // go through all surrounding points
    for(auto it = index_list.begin() + 1; it < index_list.end(); ++it) kinetic += dchi[*it] + chi_A[*it];
        

    // The right hand side of the PDE 
    T source = (1+D*rho_0)/a_3 - pow(chi_0/chi, 1-n);
    source *= chi_prefactor; // beta*rho_m,0 / Mpl, [(h/Mpc)^3]

    // source term arising rom restricting the equation down to the lower level
    if( level > 0 && addsource ) source += MultiGridSolver<3, T>::get_multigrid_source(level, i);

    // The discretized equation of motion L_{ijk...}(phi) = 0
    return kinetic/(h*h) - source;
}

// Differential of the L operator: dL_{ijk...}/dphi_{ijk...}
template<typename T>
T  ChiSolver<T>::dl_operator(unsigned int level, std::vector<unsigned int>& index_list){
    const unsigned int i = index_list[0];
    T *const dchi = MultiGridSolver<3, T>::get_y(level); // chi = dchi + chi_A
    T *const chi_A = MultiGridSolver<3,T>::get_external_field(level, 0); // solution from previus time-step
    const T chi = dchi[i] + chi_A[i]; // full solution

    // Gridspacing
    const T h = 1.0/T( MultiGridSolver<3, T>::get_N(level) );

    // Derivative of kinetic term
    const T dkinetic = -2.0*3;
    
    // Derivative of source
    const T dsource = chi_prefactor*(1-n)/chi_0*pow(chi_0/chi, 2-n);

    return dkinetic/(h*h) - dsource;
}