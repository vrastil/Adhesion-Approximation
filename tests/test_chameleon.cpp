#include <catch.hpp>

TEST_CASE( "UNIT TEST: create Multigrid and copy data to/from Mesh", "[multigrid]" )
{
    unsigned N = 32;
    srand(time(0));
    MultiGrid<3, FTYPE> grid(N);
    Mesh mesh_from(N);
    Mesh mesh_to(N);

    std::vector<unsigned> some_indices = {4, N*2+5, N*N*4+8*N+5};
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

TEST_CASE( "UNIT TEST: create and initialize ChiSolver, check bulk field ", "[chameleon]" )
{
    const unsigned N = 32;
    const FTYPE a = 0.5;
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

    for(unsigned i : some_indices) CHECK( rho_grid[0][i] == 0. );

    // initialize chi_A -- bulk field
    MultiGrid<3, FTYPE> chi_A(N);
    const unsigned N_tot = chi_A.get_Ntot();
    const FTYPE chi_0 = 2*sim.chi_opt.beta*MPL*sim.chi_opt.phi;
    const FTYPE n = sim.chi_opt.n;
    auto chi_min_ = [chi_0, a, n](FTYPE delta){ return chi_min<FTYPE>(chi_0, a, n, delta); };
    #pragma omp parallel for
    for (unsigned i = 0; i < N_tot; ++i)
    {
        chi_A[0][i] = chi_min_(rho_grid[0][i]);
    }
    chi_A.restrict_down_all();

    const FTYPE chi_bulk = chi_min(chi_0, a, n, FTYPE(0));
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