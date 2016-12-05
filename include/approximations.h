/*
Header for the approximations.cpp file
*/

typedef double(*t_power)(double, double*);

void zel_app(int, int, int, t_power, double*, double, double, int, std::string);

void frozen_flow(int, int, int, t_power, double*, double, double, int, std::string);

void truncated_zeldovich(int, int, int, t_power, double*, double, double, int, std::string, double);

void adhesiom_approximation(int, int, int, t_power, double*, double, double, double, int, std::string);

void frozen_potential(int mesh_num, int Ng, int L, t_power power_spectrum, double* parameters, double z_in, double z_out, int nt, std::string out_dir);

void mod_frozen_flow(int mesh_num, int Ng, int L, t_power power_spectrum, double* parameters, double z_in, double z_out, int nt, std::string out_dir);