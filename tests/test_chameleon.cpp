#include <catch.hpp>
#include "core_mesh.h"
#include "core_out.h"

TEST_CASE( "UNIT TEST: create Multigrid and copy data to/from Mesh", "[multigrid]" )
{
    print_unit_msg("create Multigrid and copy data to/from Mesh");

    constexpr unsigned N = 32;
    srand(time(0));
    MultiGrid<3, FTYPE> grid(N);
    Mesh mesh_from(N);
    Mesh mesh_to(N);

    const std::vector<unsigned> some_indices = {4, N*2+5, N*N*4+8*N+5};
    for(unsigned i : some_indices) mesh_from[i] = rand();

    transform_Mesh_to_Grid(mesh_from, grid);
    transform_Grid_to_Mesh(mesh_to, grid);

    for(unsigned i : some_indices) CHECK( mesh_from[i] == mesh_to[i] );
}

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

TEST_CASE( "UNIT TEST: create and initialize ChiSolver, check bulk field", "[chameleon]" )
{
    print_unit_msg("create and initialize ChiSolver, check bulk field");

    constexpr unsigned N = 32;
    constexpr FTYPE a = 0.5;
    const std::vector<unsigned> some_indices = {4, N*2+5, N*N*4+8*N+5};

    // initialize Sim_Param
    const int argc = 1;
    const char* const argv[1] = {"test"};
    Sim_Param sim(argc, argv);

    // initialize ChiSolver
    ChiSolver<FTYPE> sol(N, sim);

    // initialize overdensity -- constant density
    MultiGrid<3, FTYPE> rho_grid(N);
    {
        Mesh rho(N);
        rho.assign(0.);
        transform_Mesh_to_Grid(rho, rho_grid);
    }

    // check density
    for(unsigned i : some_indices) CHECK( rho_grid[0][i] == 0. );

    // initialize chi_A -- bulk field
    MultiGrid<3, FTYPE> chi_A(N);
    set_bulk(chi_A, rho_grid, a, sim.chi_opt);

    // check bulk field
    const FTYPE chi_bulk = chi_min(2*sim.chi_opt.beta*MPL*sim.chi_opt.phi, a, sim.chi_opt.n, FTYPE(0));
    for(unsigned i : some_indices) REQUIRE( chi_A[0][i] == Approx(chi_bulk));

    // set ChiSolver
    sol.set_time(a, sim.cosmo);
    sol.add_external_grid(&chi_A);
    sol.add_external_grid(&rho_grid);

    // check that EOM is satisfied
    std::vector<unsigned int> index_list;
    unsigned level = 0;
    for(unsigned i : some_indices)
    {
        get_neighbor_gridindex(index_list, i, N);
        REQUIRE( sol.l_operator(level, index_list, true) == Approx(0.));
    }
}

TEST_CASE( "UNIT TEST: create and initialize ChiSolver, solve sphere", "[chameleon]" )
{
    print_unit_msg("create and initialize ChiSolver, solve sphere");

    constexpr int N = 64;
    constexpr FTYPE R2 = (N/8)*(N/8);
    int ix0 = N/2, iy0 = N/2, iz0 = N/2;

    constexpr FTYPE rho_0 = 1E-8;

    // initialize Sim_Param
    const int argc = 1;
    const char* const argv[1] = {"test"};
    Sim_Param sim(argc, argv);

    // initialize ChiSolver
    ChiSolver<FTYPE> sol(N, sim);

    // initialize overdensity -- constant density in sphere of radius R, center at x0, y0, z0
    MultiGrid<3, FTYPE> rho_grid(N);
    FTYPE mean_rho;
    {
        Mesh rho(N);
        const int N_tot = rho.length;
        FTYPE R2_;
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
        REQUIRE( rho(N/2, N/2, N/2) == Approx(rho_0 - mean_rho) );
        REQUIRE( rho(0, 0, 0) == Approx(-mean_rho) );
    }

    // initialize chi_A -- bulk field
    MultiGrid<3, FTYPE> chi_A(N);
    set_bulk(chi_A, rho_grid, FTYPE(1), sim.chi_opt);

    // set ChiSolver
    sol.set_epsilon(1e-16);
    sol.set_time(1, sim.cosmo);
    sol.add_external_grid(&chi_A);
    sol.add_external_grid(&rho_grid);

    // Solve the equation
    sol.solve();

    // add solution to bulk field
    {
        FTYPE* dchi = sol.get_y();
        const unsigned N_tot = chi_A.get_Ntot();
        #pragma omp parallel for
        for (unsigned i = 0; i < N_tot; ++i) chi_A[0][i] += dchi[i];
    }
    
    // copy onto Mesh
    Mesh chi_full(N);
    transform_Grid_to_Mesh(chi_full, chi_A);

    // create directory structure and open file
    std::string out_dir = sim.out_opt.out_dir + "test_ChiSolver/";
    string file_name = out_dir + "chi_values.dat";
    create_dir(out_dir);
    Ofstream File(file_name);

    // print chi_full
    File << setprecision(12);
    File << "# N :=\t" << N << "\n";
    File << "# R :=\t" << sqrt(R2) << "\n";
    File << "# rho_0 :=\t" << rho_0 << "\n";
    File << "# phi screening :=\t" << sim.chi_opt.phi << "\n";
    File << "# phi gravitational :=\t" << R2*rho_0/(4*MPL*MPL) << "\n";
    File << "# r\tchi(r)-chi_0\n";
    FTYPE r;
    const FTYPE chi_0 = chi_min(2*sim.chi_opt.beta*MPL*sim.chi_opt.phi, 1.0, sim.chi_opt.n, -mean_rho);
    //const FTYPE chi_0 = 2*sim.chi_opt.beta*MPL*sim.chi_opt.phi;

    for (int i = 0; i < N; ++i)
    {
        for (int j = i; j < i + 2; ++j)
        {
            for (int k = i; k < i + 2; ++k)
            {
                r = sqrt(pow(i - ix0, 2) + pow(j - iy0, 2) + pow(k - iz0, 2));
                File << r << "\t" << chi_full(i, j, k) - chi_0 << "\n";
            }
        }   
    }
}