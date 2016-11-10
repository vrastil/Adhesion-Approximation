/*
Header for the grid_fce.cpp file
*/

typedef double(*t_power)(double, double*);

/* POWER SPECTRUM FUNCTIONS */
void power_spectrum_norm(double*);
double power_spectrum_T(double, double*);
double power_spectrum(double, double*);
double flat_power_spectrum(double, double*);
double single_power_spectrum(double, double*);

void gen_sin_rho_k(int, int, double*, fftw_plan, double, int);

void pwr_spec (int, int, int, double**, double*, fftw_plan, int, const int);
void pwr_spec(int, int, int, double*, fftw_plan, int, const int);
void comp_power_spec(int, fftw_complex*, int, int, int);
void pwr_spec_0_nt(int, int, double*, fftw_complex*, int, const int);
void get_rho_par(int, int, double**, double*, int, const int);

/* FFT FUNCTIONS */
void fftw_execute_dft_c2r_triple(fftw_plan, fftw_complex**, double**);
void normalize_FFT_FORWARD_nt(int, fftw_complex*, int);
void normalize_FFT_BACKWARD_nt(int, double*, int);
void normalize_FFT_FORWARD_nt(int, fftw_complex**, int);
void normalize_FFT_BACKWARD_nt(int, double**, int);

/* GENERAL FUNCTIONS */
int isPowerOfTwo (unsigned int);
const char *humanSize(uint64_t);
void get_coord_fftw_complex(int, int, int*);
void get_coord_double(int, int, int*);
int get_per(int, int);
double get_per(double, int);
double get_k_sq(int , int);
double get_k_sq(int , int*);
int get_pos(int, int);
void assign_fc(double*, double*, double &, int, const int, bool);

/* RANDOM FIELD FUNCTIONS */
double genrand_real(int*);
double generateGaussianNoise(int*, const double&, const double &);
void gen_gauss_white_noise_nt(int, double*, const double &, const double &, int);

/* APPROXIMATION FUNCTIONS */
void gen_rho_dist_k(int , int , double*, fftw_plan, t_power , double*, double, int);
void gen_displ_k_w_nt(int, double**, double*, int);
void gen_displ_k_nt(int, double**, double*, int);
void gen_no_displ_pos_vel_nt(int, int, double***, int);
void gen_no_displ_pos_nt(int, int, double**, int);
void gen_init_expot(int, double*, double*, double, bool, fftw_plan, int);
void gen_expot(int, double*, double*, double, double, bool, fftw_plan); // <-- REWRITE TO MULTITHREAD!
void fin_diff(int, double**, double*); // <-- REWRITE TO MULTITHREAD!

/* MULTI-THREAD FUNCTIONS */
void exp2pot_nt(int, double*, double*, double, int);
void exp2pot_nt(int, double*, double, int);
void pot2exp_nt(int, double*, double*, double, int);
void pot2exp_nt(int, double*, double, int);
void init_dens_field(int, double*, int);

void gen_rho_k_nt(int, int, fftw_complex*, t_power, double*, double, int);

void upd_pos_nt(int, int, double**, double**, double, int, int);
void upd_pos_vel0_nt(int, int, double**, double**, double, int, int);
void upd_pos_mid_nt(int, int, double***, double**, double, int, int);
void upd_pos_leapfrog_nt(int par_num, int mesh_num, double*** par_pos, double** vec_field, double db, int order, int nt);
void gen_pot_k_nt(int, double*, int);
void gen_displ_pos_ng_nt(int, int, double**, double**, double, int);


/* ALLOCATION FUNCTIONS */
void** alloc_zeldovich(int, int, int, int);
void dealloc_zeldovich(void**);
void** alloc_adhesion(int, int, int, int);
void dealloc_adhesion(void**);
void** alloc_frozen_pot(int mesh_num, int par_num, int bin_num, int nt);
void dealloc_frozen_pot(void** arrays);