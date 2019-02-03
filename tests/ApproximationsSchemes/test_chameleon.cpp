#include <catch.hpp>
#include "../test.hpp"
#include "chameleon.cpp"

namespace{

template<typename T>
T mean(const std::vector<T>& data)
{
    T tmp(0);
	
	#pragma omp parallel for reduction(+:tmp)
	for (auto it = data.begin(); it < data.end(); ++it) tmp += *it;
	
	return tmp / data.size();
}

FTYPE_t mean(const Mesh& data){ return mean(data.data); }

void get_neighbor_gridindex(std::vector<size_t>& index_list, size_t i, size_t N)
{
    index_list = std::vector<size_t>(7);
    index_list[0] = i;
    for(size_t j = 0, n = 1; j < 3; j++, n *= N){
        size_t ii = i/n % N;
        size_t iminus = ii >= 1   ? ii - 1 : N - 1;
        size_t iplus  = ii <= N-2 ? ii + 1 : 0;
        index_list[2*j+1] = i + (iminus - ii) * n;
        index_list[2*j+2] = i + (iplus - ii) * n;
    }
}

void init_overdensity(Sim_Param& sim, Mesh& rho, MultiGrid<3, CHI_PREC_t>& rho_grid)
{
    const size_t N = sim.test_opt.N_grid;
    const FTYPE_t R = sim.test_opt.R_sphere;
    const FTYPE_t rho_0 = sim.test_opt.rho_sphere;
    const size_t N_tot = rho.length;

    FTYPE_t mean_rho;
    FTYPE_t R2_;
    const FTYPE_t R2 = R*R;
    size_t ix0 = N/2, iy0 = N/2, iz0 = N/2;
    
    #pragma omp parallel for private(R2_)
    for (size_t ix = 0; ix < N; ++ix)
    {
        for (size_t iy = 0; iy < N; ++iy)
        {
            for (size_t iz = 0; iz < N; ++iz)
            {
                R2_ = pow2(ix - ix0) + pow2(iy - iy0) + pow2(iz - iz0);
                rho(ix, iy, iz) = R2_ <= R2 ? rho_0 : 0;
            }
        }
    }
    mean_rho = mean(rho);

    if (mean_rho > 1) throw std::out_of_range("invalid values of sphere parameters (overdensity: " + std::to_string(-mean_rho) + " < -1)");

    rho -= mean_rho;
    transform_Mesh_to_MultiGrid(rho, rho_grid);

    // check density
    REQUIRE( mean(rho) == Approx(0.) );
    REQUIRE( rho(ix0, iy0, iz0) == Approx(rho_0 - mean_rho) );
    REQUIRE( rho(0, 0, 0) == Approx(-mean_rho) );

    // update background underdensity
    sim.test_opt.rho_b = -mean_rho;
}

void get_grav_pot(Mesh& rho, const FFTW_PLAN_TYPE& p_F, const FFTW_PLAN_TYPE& p_B, FTYPE_t box_size, FTYPE_t phi_prefactor)
{
    fftw_execute_dft_r2c(p_F, rho);
    gen_pot_k(rho);
    fftw_execute_dft_c2r(p_B, rho);
    rho *= pow(box_size / rho.N, 2) // computing units (laplace)
        *  phi_prefactor; // prefactor of poisson equation
}

template <typename T> static int sgn(T val)
{
	return (T(0) < val) - (val < T(0));
}

void print_mesh(const std::string& file_name, const Mesh& pot, const FTYPE_t mod = -1 - CHI_MIN)
{
    Ofstream File(file_name);
    size_t N = pot.N;
    int ix0 = N/2, iy0 = N/2, iz0 = N/2;

    auto print_r_chi = [&File, ix0, iy0, iz0,&pot, mod]
                       (int i, int j, int k)
                        {
                            Vec_3D<size_t> pos(i, j, k); //< to call periodic position on Mesh
                            FTYPE_t r = sqrt(pow2(i - ix0) + pow2(j - iy0) + pow2(k - iz0));
                            File << r << "\t" << pot(pos) + mod  << "\n";
                        };

    for (size_t i = 0; i < N; ++i)
    {
        print_r_chi(i, i, i);
        print_r_chi(i, i, i + 1);
        print_r_chi(i, i + 1, i + 1);
    }
}

} // namespace <unique> end

TEST_CASE( "UNIT TEST: create Multigrid and copy data to/from Mesh", "[multigrid]" )
{
    print_unit_msg("create Multigrid and copy data to/from Mesh");

    constexpr size_t N = 32;
    srand(time(0));
    MultiGrid<3, FTYPE_t> grid(N);
    Mesh mesh_from(N);
    Mesh mesh_to(N);

    const std::vector<size_t> some_indices = // Mesh indices = N*N*(N+2)
                                            {4, // (0, 0, 4)
                                            (N+2)*2+5, // (0, 2, 5)
                                            N*(N+2)*4+8*(N+2)+5}; // (4, 8, 5)
    for(size_t i : some_indices) mesh_from[i] = 2.0*rand()/RAND_MAX - 1;

    transform_Mesh_to_MultiGrid(mesh_from, grid);
    CHECK(min(mesh_from) == min(grid));
    transform_MultiGrid_to_Mesh(mesh_to, grid);

    for(size_t i : some_indices) CHECK( mesh_from[i] == mesh_to[i] );
    CHECK(min(mesh_from) == min(grid));
    CHECK(min(mesh_to) == min(grid));
}

TEST_CASE( "UNIT TEST: create and initialize ChiSolver, check bulk field", "[chameleon]" )
{
    print_unit_msg("create and initialize ChiSolver, check bulk field");

    constexpr size_t N = 32;
    constexpr FTYPE_t a = 0.5;
    constexpr FTYPE_t rho_0 = 0.3;
    const std::vector<size_t> some_indices = // Mesh indices = N*N*(N+2)
                                            {4, // (0, 0, 4)
                                            (N+2)*2+5, // (0, 2, 5)
                                            N*(N+2)*4+8*(N+2)+5}; // (4, 8, 5)

    // initialize Sim_Param
    const int argc = 1;
    const char* const argv[1] = {"test"};
    Sim_Param sim(argc, argv);

    // check thowing of exception with unphysical values of chameleon parameters
    auto constructor = [&sim](){ ChiSolver<FTYPE_t> _sol(N, sim, false); };
    sim.chi_opt.n = -1;
    CHECK_THROWS_AS( constructor(), std::out_of_range );
    sim.chi_opt.n = 0.5;
    sim.chi_opt.beta = -1;
    CHECK_THROWS_AS( constructor(), std::out_of_range );

    // initialize ChiSolver
    sim.chi_opt.beta = 1/sqrt(6.);
    ChiSolver<FTYPE_t> sol(N, sim, false);

    // check thowing of exceptions in uninitialised state
    CHECK_THROWS_AS( sol.set_bulk_field(), std::out_of_range );

    // initialize overdensity -- constant density
    MultiGrid<3, FTYPE_t> rho_grid(N);
    {
        Mesh rho(N);
        rho.assign(rho_0);
        transform_Mesh_to_MultiGrid(rho, rho_grid);
    }

    // check density
    for(size_t i : some_indices) CHECK( rho_grid[0][i] == rho_0 );

    // set ChiSolver -- bulk field
    sol.set_time(a, sim.cosmo);
    sol.add_external_grid(&rho_grid);
    sol.set_bulk_field();

    // check bulk field
    {
        const FTYPE_t chi_bulk = sol.chi_min(rho_0);
        FTYPE_t const* const chi = sol.get_y();
        for(size_t i : some_indices) REQUIRE( chi[i] == Approx(chi_bulk));
    }

    // check that EOM is satisfied
    std::vector<size_t> index_list;
    size_t level = 0;
    const FTYPE_t h = 1.0/FTYPE_t( sol.get_N(level) );
    for(size_t i : some_indices)
    {
        get_neighbor_gridindex(index_list, i, N);
        CHECK( sol.l_operator(level, index_list, true, h) == Approx(0.));
    }
}

TEST_CASE( "UNIT TEST: create and initialize ChiSolver, solve sphere", "[chameleon]" )
{
    print_unit_msg("create and initialize ChiSolver, solve sphere");

    // initialize Sim_Param
    const char* const argv[1] = {"test"};
    Sim_Param sim(1, argv);
    const size_t N = sim.test_opt.N_grid;
    const size_t N_min = sim.test_opt.N_min;
    const bool verbose = sim.test_opt.verbose;

    // initialize ChiSolver
    ChiSolver<CHI_PREC_t> sol(N, N_min, sim, verbose);

    // initialize overdensity -- constant density in sphere of radius R, center at x0, y0, z0
    MultiGrid<3, CHI_PREC_t> rho_grid(N);
    Mesh rho(N);
    init_overdensity(sim, rho, rho_grid);

    // FFTW preparation
    if (!FFTW_PLAN_OMP_INIT()){
        throw std::runtime_error("Errors during multi-thread initialization");
    }
    FFTW_PLAN_OMP(sim.run_opt.nt);
    const FFTW_PLAN_TYPE p_F = FFTW_PLAN_R2C(N, N, N, rho.real(), rho.complex(), FFTW_ESTIMATE);
    const FFTW_PLAN_TYPE p_B = FFTW_PLAN_C2R(N, N, N, rho.complex(), rho.real(), FFTW_ESTIMATE);

    // set ChiSolver
    sol.set_time(1, sim.cosmo);
    sol.set_def_convergence();
    sol.add_external_grid(&rho_grid);
    sol.set_ngs_sweeps(sim.test_opt.fine_sweeps, sim.test_opt.coarse_sweeps); //< fine, coarse

    // compute gravitational potential
    Mesh phi_pot(rho); //< copy density
    get_grav_pot(phi_pot, p_F, p_B, sim.box_opt.box_size, sol.get_phi_prefactor());

    // get linear prediction
    sol.set_linear(rho, p_F, p_B);
    sol.set_screened();

    // full solution on Mesh
    Mesh chi_full(N);

    // create directory structure
    std::string out_dir = sim.out_opt.out_dir + "test_ChiSolver/";
    remove_all_files(out_dir);
    create_dir(out_dir);
    sim.print_info(out_dir, "test");

    // print gravitational potential
    print_mesh(out_dir + "grav_pot.dat", phi_pot, 0);

    // Solve the equation -- full V-cycles
    BOOST_LOG_TRIVIAL(info) << "Starting V-cycles with intermediate output...\n";
    
    int istep = 0;
    const size_t step_per_iter = sim.test_opt.step_per_iter;
    while(1)
    {
        // istep, max_step (use int for the last iteration -- negative max_step)
        istep += sol.get_istep();
        int max_step = istep + step_per_iter > sim.test_opt.max_steps ? sim.test_opt.max_steps - istep : step_per_iter;

        // print chi_full
        transform_MultiGridSolver_to_Mesh(chi_full, sol);
        print_mesh(out_dir + "chi_istep_" + std::to_string(istep) + ".dat", chi_full);

        // check max_step and convergence
        if ((max_step <= 0) || ((istep != 0) && (sol.get_istep() < step_per_iter))) break;

        // set & solve
        sol.set_maxsteps(max_step);
        sol.solve();
    }

    // FFTW CLEANUP
	FFTW_DEST_PLAN(p_F);
    FFTW_DEST_PLAN(p_B);
	FFTW_PLAN_OMP_CLEAN();
}