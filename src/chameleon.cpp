#include "core_app.h"
#include "core_power.h"
#include "core_mesh.h"
#include "core_out.h"
#include "chameleon.hpp"

using namespace std;

constexpr FTYPE_t MPL = (FTYPE_t)1; // chi in units of Planck mass
// (FTYPE_t)2.435E18; // reduced Planck mass, [GeV/c^2]
// constexpr FTYPE_t fm_to_Mpc = (FTYPE_t)3.2407792896664E-38; // 1 fm = ? Mpc
// constexpr FTYPE_t hbarc = (FTYPE_t)197.327053; // reduced Planck constant times speed of light, [MeV fm]
// constexpr FTYPE_t hbarc_cosmo = hbarc*FTYPE_t(1E-9)*fm_to_Mpc; // [GeV Mpc]
// constexpr FTYPE_t G_N = hbarc_cosmo*FTYPE_t(6.70711*1E-39); // gravitational constant, [GeV Mpc / (GeV/c^2)^2]
 constexpr FTYPE_t c_kms = 1;
 // (FTYPE_t)299792.458; // speed of light [km / s]

template<typename T>
static void transform_Mesh_to_Grid(const Mesh& mesh, Grid<3, T> &grid)
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
static void transform_Mesh_to_Grid(const Mesh& mesh, MultiGrid<3, T> &grid)
{
    transform_Mesh_to_Grid(mesh, grid.get_grid());
    grid.restrict_down_all();
}

template<typename T>
static void transform_Grid_to_Mesh(Mesh& mesh, const Grid<3, T> &grid)
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
static void transform_Grid_to_Mesh(Mesh& mesh, const MultiGrid<3, T> &grid)
{
    transform_Grid_to_Mesh(mesh, grid.get_grid());
}

template<typename T>
static void transform_Grid_to_Mesh(Mesh& mesh, const MultiGridSolver<3, T> &sol)
{
    transform_Grid_to_Mesh(mesh, sol.get_grid());
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
    a_0 = a; // store old time
    a = a_; // set new time
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

template<typename T>
void ChiSolver<T>::set_next_guess(const Cosmo_Param& cosmo)
{
    if(!a_0) return; // do not set for initial conditions

    const T a_0_3 = pow(a_0, 3);
    const T D_0 = growth_factor(a_0, cosmo);

    T* const f = MultiGridSolver<3, T>::get_y(); // initial guess
    T const* const rho = MultiGridSolver<3,T>::get_external_field(0, 0); // overdensity
    const unsigned N_tot = MultiGridSolver<3, T>::get_Ntot();

    #pragma omp parallel for
    for (unsigned i = 0; i < N_tot; ++i)
    {
        f[i] *= pow(a_3*(1+D_0*rho[i])/(a_0_3*(1+D*rho[i])), 1/(1-n));
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

    // The right hand side of the PDE 
    T source;
    if(chi[i] <= 0)
    { // if the solution is overshot, try bulk field
        MultiGridSolver<3, T>::get_y(level)[i] = chi_min(rho_0);
        source = 0;
    } else {
        source = (1+D*rho_0)/a_3 - pow(chi_0/chi[i], 1-n);
        source *= chi_prefactor; // beta*rho_m,0 / Mpl^2, [dimensionless]
    }

    // Compute the standard kinetic term [D^2 phi] (in 1D this is phi''_i =  phi_{i+1} + phi_{i-1} - 2 phi_{i} )
    T kinetic = -2*chi[i]*3; // chi, '-2*3' is factor in 3D discrete laplacian
    // go through all surrounding points
    for(auto it = index_list.begin() + 1; it < index_list.end(); ++it) kinetic += chi[*it];

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
    App_Var<Particle_v<FTYPE_t>>(sim, app_str), sol(sim.box_opt.mesh_num, sim, false), drho(sim.box_opt.mesh_num)
{
    // EFFICIENTLY ALLOCATE VECTOR OF MESHES
    chi_force.reserve(3);
    for(size_t i = 0; i < 3; i++){
        chi_force.emplace_back(sim.box_opt.mesh_num);
    }
    memory_alloc += sizeof(FTYPE_t)*app_field[0].length*app_field.size();

    // SET CHI SOLVER
    sol.add_external_grid(&drho);
    sol.set_maxsteps(8);
}

template<typename T>
static T min(const std::vector<T>& data)
{
    return *std::min_element(data.begin(), data.end());
}

static FTYPE_t min(const Mesh& data){
    return min(data.data);
}

void App_Var_chi::save_init_drho_k(const Mesh& dro_k, Mesh& aux_field)
{
    // do not overwrite aux_field if Mesh of different type
    if (dro_k.N != aux_field.N) throw std::range_error("Meshes of a different sizes!");

    // save initial overdensity
    cout << "Storing initial density distribution...\n";
    aux_field = dro_k;
    fftw_execute_dft_c2r(p_B, aux_field);
    transform_Mesh_to_Grid(aux_field, drho);

    cout << "Minimal value of initial overdensity:\t" << min(aux_field) << "\n";

    // set initial guess to bulk field
    cout << "Setting initial guess for chameleon field...\n";
    sol.set_time(b, sim.cosmo);
    sol.set_initial_guess();
}

static void pwr_spec_chi_k(const Mesh &chi_k, Mesh& power_aux, const FTYPE_t prefactor = 1)
{
    Vec_3D<int> k_vec;
    FTYPE_t k;
    const unsigned NM = chi_k.N;
    const unsigned half_length = chi_k.length / 2;

	#pragma omp parallel for private(k_vec, k)
	for(unsigned i=0; i < half_length;i++)
	{
		get_k_vec(NM, i, k_vec);
        k = k_vec.norm();
        power_aux[2*i] = (pow2(chi_k[2*i]) + pow2(chi_k[2*i+1]))*pow(k, 4)*prefactor;
		power_aux[2*i+1] = k;
	}
}

void App_Var_chi::print_output()
{/* Print standard output */
    App_Var<Particle_v<FTYPE_t>>::print_output();

    /* Chameleon power spectrum */
    if (sim.out_opt.print_pwr)
    {
        transform_Grid_to_Mesh(chi_force[0], sol); // get solution
        fftw_execute_dft_r2c(p_F, chi_force[0]); // get chi(k)
        pwr_spec_chi_k(chi_force[0], chi_force[0], 1/sol.chi_prefactor); // get chi(k)^2 * k^4
        gen_pow_spec_binned(sim, chi_force[0], pwr_spec_binned); // get average Pk
        print_pow_spec(pwr_spec_binned, out_dir_app, "_chi" + z_suffix()); // print
    }
}
template class ChiSolver<CHI_PREC_t>;

#ifdef TEST
#include "test_chameleon.cpp"
#endif