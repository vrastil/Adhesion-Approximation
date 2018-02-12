#include "stdafx.h"
#include "core.h"
#include "core_power.h"
#include "chameleon.hpp"

using namespace std;

constexpr FTYPE MPL = (FTYPE)1; // chi in units of Planck mass
// (FTYPE)2.435E18; // reduced Planck mass, [GeV/c^2]
// constexpr FTYPE fm_to_Mpc = (FTYPE)3.2407792896664E-38; // 1 fm = ? Mpc
// constexpr FTYPE hbarc = (FTYPE)197.327053; // reduced Planck constant times speed of light, [MeV fm]
// constexpr FTYPE hbarc_cosmo = hbarc*FTYPE(1E-9)*fm_to_Mpc; // [GeV Mpc]
// constexpr FTYPE G_N = hbarc_cosmo*FTYPE(6.70711*1E-39); // gravitational constant, [GeV Mpc / (GeV/c^2)^2]
 constexpr FTYPE c_kms = (FTYPE)299792.458; // speed of light [km / s]

template<typename T>
inline T chi_min(T chi_0, T a, T n, T delta){ return chi_0*pow(pow(a, 3)/(1+delta), 1/(1-n)); }

template<typename T>
void transform_Mesh_to_Grid(const Mesh& mesh, MultiGrid<3, T> &grid)
{/* copy data in Mesh 'N*N*(N+2)' onto MultiGrid 'N*N*N' */
    unsigned int ix, iy, iz;
    const unsigned N_tot = grid.get_Ntot();
    const unsigned N = grid.get_N();

    if (mesh.N != N) throw std::range_error("Mesh of a different size than MultiGrid!");

    #pragma omp parallel for private(ix, iy, iz)
    for (unsigned i = 0; i < N_tot; ++i)
    {
        ix = i % N;
        iy = i / N % N;
        iz = i / (N*N) % N;
        grid[0][i] = mesh(ix, iy, iz);
    }
    grid.restrict_down_all();
}

template<typename T>
void transform_Grid_to_Mesh(Mesh& mesh, const MultiGrid<3, T> &grid)
{/* copy data in MultiGrid 'N*N*N' onto Mesh 'N*N*(N+2)' */
    unsigned int ix, iy, iz;
    const unsigned N_tot = grid.get_Ntot();
    const unsigned N = grid.get_N();

    if (mesh.N != N) throw std::range_error("Mesh of a different size than MultiGrid!");

    #pragma omp parallel for private(ix, iy, iz)
    for (unsigned i = 0; i < N_tot; ++i)
    {
        ix = i % N;
        iy = i / N % N;
        iz = i / (N*N) % N;
        mesh(ix, iy, iz) = grid[0][i];
    }
}

template<typename T>
ChiSolver<T>::ChiSolver(unsigned int N, int Nmin, const Sim_Param& sim, bool verbose):
    MultiGridSolver<3, T>(N, Nmin, verbose), n(sim.chi_opt.n), chi_0(2*sim.chi_opt.beta*MPL*sim.chi_opt.phi),
    chi_prefactor( // beta*rho_m,0 / Mpl^2 + computing units, e.g. dimensionless
        3*sim.chi_opt.beta*sim.cosmo.Omega_m
        *pow(sim.cosmo.H0 // Hubble constant
        / (sim.cosmo.h * c_kms) // units factor for 'c = 1' and [L] = Mpc / h
        * sim.x_0() // dimension factor for laplacian
        ,2))
{}

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
    source *= chi_prefactor; // beta*rho_m,0 / Mpl^2, [dimensionless]

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

template<typename T>
void set_bulk(MultiGrid<3, T>& chi_A, const MultiGrid<3, T>& rho, const T a, const Chi_Opt& chi_opt)
{
    const unsigned N_tot = chi_A.get_Ntot();
    const T chi_0 = 2*chi_opt.beta*MPL*chi_opt.phi;
    const T n = chi_opt.n;

    #pragma omp parallel for
    for (unsigned i = 0; i < N_tot; ++i)
    {
        chi_A[0][i] = chi_min(chi_0, a, n, rho[0][i]);
    }

    chi_A.restrict_down_all();
}


template void set_bulk(MultiGrid<3, FTYPE>&, const MultiGrid<3, FTYPE>&, const FTYPE, const Chi_Opt&);

#ifdef TEST
#include "test_chameleon.cpp"
#endif