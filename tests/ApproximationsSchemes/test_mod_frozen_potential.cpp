#include <catch.hpp>
#include "../test.hpp"
#include "ApproximationsSchemes/mod_frozen_potential.cpp"

// void force_test(Sim_Param& sim)
// {
//     // 1 particle prep
//     cout << "Adjust simulation parameters to 1 particle.\n";
//     sim.par_num = 1;
//     sim.mesh_num_pwr = sim.mesh_num;
//     sim.Ng = sim.Ng_pwr = sim.mesh_num / pow(sim.par_num, 1/3.);
//     sim.print_info();
//     App_Var_FP_mod APP(sim, "test_pp");
//     APP.print_mem();

//     cout << "Place particle in the middle of the box.\n";
//     APP.particles[0] = Particle_v<FTYPE_t>(sim.mesh_num/2, sim.mesh_num/2., sim.mesh_num/2., 0, 0, 0); // middle, no velocity
//     get_rho_from_par(APP.particles, &APP.app_field[0], sim); // assign density
//     printf("Transforming density into k-sapce...\n");
//     fftw_execute_dft_r2c(APP.p_F_pwr, APP.app_field[0]); // get \rho(k)
//     gen_pot_k(APP.app_field[0], &APP.power_aux[0]); // get \phi(k)
//     gen_displ_k_S2(&APP.app_field, APP.power_aux[0], APP.sim.a);
//     printf("Computing force in q-space...\n");
//     fftw_execute_dft_c2r_triple(APP.p_B, APP.app_field);

//     printf("Creating linked list...\n");
//     APP.linked_list.get_linked_list(APP.particles);

//     for (unsigned i = 0; i <sim.par_num; i++)
//     {
//         cout << "LL[" << i << "] = " << APP.linked_list.LL[i] << "\n";
//     }
//     for (unsigned i = 0; i < APP.linked_list.HOC.length; i++){
//         if (APP.linked_list.HOC[i] != -1){
//             cout << "HOC [" << i << "] = " << APP.linked_list.HOC[i]
//             << "\tChain position = "
//             << i / (APP.linked_list.HOC.N2*APP.linked_list.HOC.N3) << "  "
//             << (i / APP.linked_list.HOC.N3)  % APP.linked_list.HOC.N2 << "  "
//             << i % APP.linked_list.HOC.N3 <<  "\n";
//         }
//     }

//     Vec_3D<FTYPE_t> f_long, f_short, f_total, dr_vec, cur_pos;
//     FTYPE_t dr;
//     FTYPE_t m = pow(sim.Ng, 3);

//     // cout << "\n\nr_vec\t\tr\t|\tshort\t\t\t\tlong\t\t\t\ttotal\t\t\t|\ts\tl\tt\tm/4PIr2\n";
//     string file_name = "/home/vrastil/Documents/GIT/Adhesion-Approximation/output/test_runs/test_pp_run/data_rs_" + to_string(sim.rs) + ".dat";
//     ofstream File(file_name);
//     File << "#r\ts\tl\tt\tm/4PIr2\n";
//     File << setprecision(8);
//     FTYPE_t inc = sim.rs/3.;
//     for (FTYPE_t i = 5; i < sim.mesh_num-5; i += inc)
//     {
//         f_long.fill(0.);
//         f_short.fill(0.);
//         cur_pos = Vec_3D<FTYPE_t>(i, i, i);
//         dr_vec = get_sgn_distance(APP.particles[0].position, cur_pos, sim.mesh_num);
//         dr = dr_vec.norm();
//         if (dr > 5*sim.rs) inc = sim.rs/3.;
//         else if (dr > sim.rs) inc = sim.rs/10.;
//         else inc = sim.rs/40.;

//         assign_from(APP.app_field, cur_pos, &f_long);
//         force_short(sim, 1, APP.linked_list, APP.particles, cur_pos, &f_short);
//         f_total = f_short + f_long;

//         if (dr != 0){
//             File << scientific;
//             //File << dr_vec[0] << ", ...\t";
//             //File << dr << "\t|\t";
//             File << dr << "\t";
//         //    File << f_short[0] << " " << f_short[1] << " " << f_short[2] << "\t\t";
//         //    File << f_long[0] << " " << f_long[1] << " " << f_long[2] << "\t\t";
//         //    File << f_total[0] << " " << f_total[1] << " " << f_total[2] << "\t|\t";
//             File << f_short.norm() << "\t" << f_long.norm() << "\t" << f_total.norm() << "\t";
//             File << m/(4*PI*dr*dr) << "\n";
//         }
//     }
//     File.close();
// }