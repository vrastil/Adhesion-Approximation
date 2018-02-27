#include "stdafx.h"
#include "params.hpp"
#include "core_power.h"
#include "core_mesh.h"
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
        iz = i / (N*N) % N;
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

    if (mesh.N != N) throw std::range_error("Mesh of a different size than MultiGrid!");

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
void transform_Grid_to_Mesh(Mesh& mesh, const MultiGrid<3, T> &grid)
{
    transform_Grid_to_Mesh(mesh, grid.get_grid());
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
{
    if ((n <= 0) || (n >= 1)) throw out_of_range("invalid value of chameleon power-law potential exponent");
}

template<typename T>
void ChiSolver<T>::set_time(T a_, const Cosmo_Param& cosmo)
{
    a = a_;
    a_3 = pow(a_, 3);
    D = growth_factor(a_, cosmo);
}

template<typename T>
void ChiSolver<T>::set_initial_guess()
{
    if (!a_3) throw out_of_range("invalid value of scale factor");
    if (!MultiGridSolver<3, T>::get_external_field_size()) throw out_of_range("initial overdensity not set"); 

    T* const f = MultiGridSolver<3, T>::get_y(); // initial guess
    T const* const rho = MultiGridSolver<3,T>::get_external_field(0, 0); // overdensity
    const unsigned N_tot = MultiGridSolver<3, T>::get_Ntot();

    #pragma omp parallel for
    for (unsigned i = 0; i < N_tot; ++i)
    {
        f[i] = chi_min(rho[i]);
    }
}

// The dicretized equation L(phi)
template<typename T>
T  ChiSolver<T>::l_operator(unsigned int level, std::vector<unsigned int>& index_list, bool addsource)
{ 
    const unsigned int i = index_list[0];
    
    // Gridspacing
    const T h = 1.0/T( MultiGridSolver<3, T>::get_N(level) );

    // Solution and pde-source grid at current level
    T const* const chi = MultiGridSolver<3, T>::get_y(level); // solution
    const T rho_0 = MultiGridSolver<3,T>::get_external_field(level, 0)[i]; // initial overdensity

    // Compute the standard kinetic term [D^2 phi] (in 1D this is phi''_i =  phi_{i+1} + phi_{i-1} - 2 phi_{i} )
    T kinetic = -2*chi[i]*3; // chi, '-2*3' is factor in 3D discrete laplacian
    // go through all surrounding points
    for(auto it = index_list.begin() + 1; it < index_list.end(); ++it) kinetic += chi[*it];
        
    // The right hand side of the PDE 
    T source = (1+D*rho_0)/a_3 - pow(chi_0/chi[i], 1-n);
    source *= chi_prefactor; // beta*rho_m,0 / Mpl^2, [dimensionless]

    // source term arising rom restricting the equation down to the lower level
    if( level > 0 && addsource ) source += MultiGridSolver<3, T>::get_multigrid_source(level, i);

    // The discretized equation of motion L_{ijk...}(phi) = 0
    return kinetic/(h*h) - source;
}

// Differential of the L operator: dL_{ijk...}/dphi_{ijk...}
template<typename T>
T  ChiSolver<T>::dl_operator(unsigned int level, std::vector<unsigned int>& index_list){
    // solution
    const T chi = MultiGridSolver<3, T>::get_y(level)[ index_list[0] ];

    // Gridspacing
    const T h = 1.0/T( MultiGridSolver<3, T>::get_N(level) );

    // Derivative of kinetic term
    const T dkinetic = -2.0*3;
    
    // Derivative of source
    const T dsource = chi_prefactor*(1-n)/chi_0*pow(chi_0/chi, 2-n);

    return dkinetic/(h*h) - dsource;
}

App_Var_chi::App_Var_chi(const Sim_Param &sim, std::string app_str):
    App_Var<Particle_v<FTYPE>>(sim, app_str), sol(sim.box_opt.mesh_num, sim), drho(sim.box_opt.mesh_num)
{
    sol.set_epsilon(2e-5*sim.chi_opt.phi);
    sol.add_external_grid(&drho);
}

void  App_Var_chi::save_init_drho_k(const Mesh& dro_k, Mesh& aux_field)
{
    // do not overwrite aux_field if Mesh of different type
    if (dro_k.N != aux_field.N) throw std::range_error("Meshes of a different sizes!");

    // save initial overdensity
    cout << "Storing initial density distribution...\n";
    aux_field = dro_k;
    fftw_execute_dft_c2r(p_B, aux_field);
    transform_Mesh_to_Grid(aux_field, drho);

    // set initial guess to bulk field
    cout << "Setting initial guess for chameleon field...\n";
    sol.set_time(b, sim.cosmo);
    sol.set_initial_guess();
}

#ifdef TEST
#include "test_chameleon.cpp"
#endif