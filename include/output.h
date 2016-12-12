/*
Header for the output.cpp file
*/
typedef double(*t_power)(double, double*);

//std::string work_dir(std::string);

// REWRITE TO MULTITHREAD!

void work_dir_over(std::string);
void print_par_pos(int, int, int, double**, std::string);
void print_par_pos_cut(int, int, int, double**, std::string, std::string);
void print_par_pos_cut_small(int, int, int, double**, std::string, std::string);
void print_pow_spec(int, fftw_complex*, std::string, std::string, double, double, int);
void print_pow_spec_diff(fftw_complex*, fftw_complex*, std::string, std::string, int, double);
void gen_pow_spec_binned(int, fftw_complex*, fftw_complex*, double, double, int);

void get_track_par_id(int, int*, int);
void update_track_par(double***, double**, int, int*, int);
void print_track_par(double***, int, int, int, int, std::string, std::string);
void print_rho_map(int, int,  double*, std::string, std::string);

typedef double t_supp[2];
void upd_suppresion(t_supp* supp, fftw_complex* pwr_spec_binned, fftw_complex* pwr_spec_binned_0, double b, int step);
void print_suppression(t_supp* supp, std::string out_dir, int num_step, std::string suffix);
void gen_dens_binned(int mesh_num, int Ng, double* rho, int* dens_binned, int bin_num);
void print_dens_bin(int* dens_binned, std::string out_dir, std::string suffix, int bin_num, int mesh_num);