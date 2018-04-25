#include <catch.hpp>
#include "../test.hpp"
#include "ApproximationsSchemes/chameleon.cpp"

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

void get_neighbor_gridindex(std::vector<unsigned int>& index_list, unsigned int i, unsigned int N)
{
    index_list = std::vector<unsigned int>(7);
    index_list[0] = i;
    for(unsigned int j = 0, n = 1; j < 3; j++, n *= N){
        unsigned int ii = i/n % N;
        unsigned int iminus = ii >= 1   ? ii - 1 : N - 1;
        unsigned int iplus  = ii <= N-2 ? ii + 1 : 0;
        index_list[2*j+1] = i + (iminus - ii) * n;
        index_list[2*j+2] = i + (iplus - ii) * n;
    }
}

} // namespace <unique> end

TEST_CASE( "UNIT TEST: create Multigrid and copy data to/from Mesh", "[multigrid]" )
{
    print_unit_msg("create Multigrid and copy data to/from Mesh");

    constexpr unsigned N = 32;
    srand(time(0));
    MultiGrid<3, FTYPE_t> grid(N);
    Mesh mesh_from(N);
    Mesh mesh_to(N);

    const std::vector<unsigned> some_indices = // Mesh indices = N*N*(N+2)
                                            {4, // (0, 0, 4)
                                            (N+2)*2+5, // (0, 2, 5)
                                            N*(N+2)*4+8*(N+2)+5}; // (4, 8, 5)
    for(unsigned i : some_indices) mesh_from[i] = 2.0*rand()/RAND_MAX - 1;

    transform_Mesh_to_Grid(mesh_from, grid);
    CHECK(min(mesh_from) == min(grid));
    transform_Grid_to_Mesh(mesh_to, grid);

    for(unsigned i : some_indices) CHECK( mesh_from[i] == mesh_to[i] );
    CHECK(min(mesh_from) == min(grid));
    CHECK(min(mesh_to) == min(grid));
}

TEST_CASE( "UNIT TEST: create and initialize ChiSolver, check bulk field", "[chameleon]" )
{
    print_unit_msg("create and initialize ChiSolver, check bulk field");

    constexpr unsigned N = 32;
    constexpr FTYPE_t a = 0.5;
    constexpr FTYPE_t rho_0 = 0.3;
    const std::vector<unsigned> some_indices = // Mesh indices = N*N*(N+2)
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
        transform_Mesh_to_Grid(rho, rho_grid);
    }

    // check density
    for(unsigned i : some_indices) CHECK( rho_grid[0][i] == rho_0 );

    // set ChiSolver -- bulk field
    sol.set_time(a, sim.cosmo);
    sol.add_external_grid(&rho_grid);
    sol.set_bulk_field();

    // check bulk field
    {
        const FTYPE_t chi_bulk = sol.chi_min(rho_0);
        FTYPE_t const* const chi = sol.get_y();
        for(unsigned i : some_indices) REQUIRE( chi[i] == Approx(chi_bulk));
    }

    // check that EOM is satisfied
    std::vector<unsigned int> index_list;
    unsigned level = 0;
    const FTYPE_t h = 1.0/FTYPE_t( sol.get_N(level) );
    for(unsigned i : some_indices)
    {
        get_neighbor_gridindex(index_list, i, N);
        CHECK( sol.l_operator(level, index_list, true, h) == Approx(0.));
    }
}

TEST_CASE( "UNIT TEST: create and initialize ChiSolver, solve sphere", "[chameleon]" )
{
    print_unit_msg("create and initialize ChiSolver, solve sphere");

    constexpr int N = 128;
    constexpr FTYPE_t R = 2;
    constexpr FTYPE_t R2 = R*R;
    int ix0 = N/2, iy0 = N/2, iz0 = N/2;

    constexpr FTYPE_t rho_0 = 0.5;

    // initialize Sim_Param
    const int argc = 1;
    const char* const argv[1] = {"test"};
    Sim_Param sim(argc, argv);

    // initialize ChiSolver
    ChiSolver<CHI_PREC_t> sol(N, sim, true);

    // initialize overdensity -- constant density in sphere of radius R, center at x0, y0, z0
    MultiGrid<3, CHI_PREC_t> rho_grid(N);
    FTYPE_t mean_rho;
    Mesh rho(N);
    {
        const int N_tot = rho.length;
        FTYPE_t R2_;
        #pragma omp parallel for private(R2_)
        for (int ix = 0; ix < N; ++ix)
        {
            for (int iy = 0; iy < N; ++iy)
            {
                for (int iz = 0; iz < N; ++iz)
                {
                    R2_ = pow(ix - ix0, 2) + pow(iy - iy0, 2) + pow(iz - iz0, 2);
                    rho(ix, iy, iz) = R2_ < R2 ? rho_0 : 0;
                }
            }
        }
        mean_rho = mean(rho);
        rho -= mean_rho;
        transform_Mesh_to_Grid(rho, rho_grid);

        // check density
        REQUIRE( mean(rho) == Approx(0.) );
        REQUIRE( rho(ix0, iy0, iz0) == Approx(rho_0 - mean_rho) );
        REQUIRE( rho(0, 0, 0) == Approx(-mean_rho) );
    }

    // FFTW preparation
    if (!FFTW_PLAN_OMP_INIT()){
        throw runtime_error("Errors during multi-thread initialization");
    }
    FFTW_PLAN_OMP(sim.run_opt.nt);
    const FFTW_PLAN_TYPE p_F = FFTW_PLAN_R2C(N, N, N, rho.real(), rho.complex(), FFTW_ESTIMATE);
    const FFTW_PLAN_TYPE p_B = FFTW_PLAN_C2R(N, N, N, rho.complex(), rho.real(), FFTW_ESTIMATE);

    // set ChiSolver
    const double err_mod = pow(sim.box_opt.box_size, 2);
    sol.set_convergence(err_mod*1e-6, 0.95, 0.1, err_mod*1e-3, 5);
    sol.set_bisection_convergence(5, 0, 1e-3);
    sol.set_time(1, sim.cosmo);
    sol.add_external_grid(&rho_grid);

    // compute gravitational potential
    Mesh phi_pot(rho); //< copy density
    fftw_execute_dft_r2c(p_F, phi_pot);
    gen_pot_k(phi_pot);
    fftw_execute_dft_c2r(p_B, phi_pot);
    phi_pot *= pow(sim.box_opt.box_size/sim.box_opt.mesh_num, 2);

    // get linear prediction
    sol.set_linear(rho, p_F, p_B, 1/(2*PI));
    sol.set_screened();


    // Solve the equation
    sol.solve();
    
    // copy onto Mesh
    Mesh chi_full(N);
    transform_Grid_to_Mesh(chi_full, sol.get_grid());

    // create directory structure and open file
    std::string out_dir = sim.out_opt.out_dir + "test_ChiSolver/";
    string file_name = out_dir + "chi_values.dat";
    create_dir(out_dir);
    Ofstream File(file_name);

    // print chi_full
    File << setprecision(18);
    File << "# N :=\t" << N << "\n";
    File << "# R :=\t" << sqrt(R2) << "\n";
    File << "# rho_0 :=\t" << rho_0 << "\n";
    File << "# phi screening :=\t" << sim.chi_opt.phi << "\n";
    File << "# phi gravitational :=\t" << 3./2.*R2*rho_0 << "\n";
    File << "# r\t(chi(r)-chi_0)\n";
    const CHI_PREC_t chi_0 = sol.chi_min(-mean_rho);
    // const FTYPE_t chi_0 = max(chi_full);

    auto print_r_chi = [&File, ix0, iy0, iz0,&chi_full, chi_0, &phi_pot]
                       (int i, int j, int k){
        File << sqrt(pow(i - ix0, 2) + pow(j - iy0, 2) + pow(k - iz0, 2)) << "\t" << chi_full(i, j, k) -1
             << "\t" << phi_pot(i, j, k) << "\n";
    };

    for (int i = 0; i < N - 1; ++i)
    {
        print_r_chi(i, i, i);
        print_r_chi(i, i, i + 1);
        print_r_chi(i, i + 1, i + 1);
    }

    // FFTW CLEANUP
	FFTW_DEST_PLAN(p_F);
    FFTW_DEST_PLAN(p_B);
	FFTW_PLAN_OMP_CLEAN();
}