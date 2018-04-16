#include <boost/filesystem.hpp>
#include "core_out.h"
#include "templates/class_particles.hpp"

namespace fs = boost::filesystem;
using namespace std;

#define BUFF_SIZE 1024 * 1024 * 16 // 16 MB buffer


Ofstream::Ofstream(string file_name) : ofstream(file_name), buf(new char[BUFF_SIZE])
{
    if (!this->is_open()) throw runtime_error("Error while opening '" + file_name + "'");
    this->rdbuf()->pubsetbuf(buf, sizeof(buf));
}

Ofstream::~Ofstream()
{
    delete[] buf;
    if (this->is_open()) this->close();
}

Ifstream::Ifstream(string file_name) : ifstream(file_name)
{
    if (!this->is_open()) throw runtime_error("Error while opening '" + file_name + "'");
}

Ifstream::~Ifstream()
{
    if (this->is_open()) this->close();
}

// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
string currentDateTime()
{
	time_t     now = time(0);
	struct tm  tstruct;
	char       buf[80];
	tstruct = *gmtime(&now);
	strftime(buf, sizeof(buf), "%y%m%d_%H%M%S", &tstruct);
	
	string returnval(buf);
    return returnval;
}

string std_out_dir(string pre_subdir, const Sim_Param &sim)
{
    return sim.out_opt.out_dir + pre_subdir + currentDateTime() + "_" + to_string(sim.box_opt.mesh_num) +"m_" +
           to_string(sim.box_opt.Ng) + "p_" + to_string(sim.box_opt.mesh_num_pwr) +"M_" + to_string((int)sim.box_opt.box_size) + "b/";
}

void create_dir(string out_dir)
{
	fs::path dir(out_dir.c_str());
	if(fs::create_directories(dir)){
        cout << "Directory Created: "<< out_dir << "\n";
    }
}

template <class T>
void print_par_pos_cut_small(const std::vector<T>& particles, const Sim_Param &sim, string out_dir, string suffix)
{
   out_dir += "par_cut/";
   string file_name = out_dir + "par_cut" + suffix + ".dat";
   Ofstream File(file_name);
   
   std::cout << "Writing small cut through the box of particles into file " << file_name << "\n";
   File << "# This file contains positions of particles in units [Mpc/h].\n";
   FTYPE_t x, y, z, dx;
   const FTYPE_t x_0 = sim.x_0();
   const unsigned N = sim.box_opt.par_num;
   for(unsigned i=0; i < N; i++)
   {
       x = particles[i].position[0];
       y = particles[i].position[1];
       z = particles[i].position[2];			
       dx = abs(y - sim.box_opt.mesh_num/2);
       if ((dx < 0.5) && (x < sim.box_opt.mesh_num/4.) && (z < sim.box_opt.mesh_num/4.))
       {
           // cut (L/4 x L/4 x 0.5)
           File << x*x_0 << "\t" << z*x_0 << "\t" << y*x_0 << "\n";
       }
   }
}

void print_pow_spec(const Data_Vec<FTYPE_t, 2> &pwr_spec_binned, string out_dir, string suffix)
{
	out_dir += "pwr_spec/";
	string file_name = out_dir + "pwr_spec" + suffix + ".dat";
	Ofstream File(file_name);
	
	cout << "Writing power spectrum into file " << file_name << "\n";
	File << "# This file contains power spectrum P(k) in units [(Mpc/h)^3] depending on wavenumber k in units [h/Mpc].\n"
	        "# k [h/Mpc]\tP(k) [(Mpc/h)^3]\n";

    File << scientific;
    const unsigned size = pwr_spec_binned.size();
	for (unsigned j = 0; j < size; j++){
        for (unsigned i = 0; i < 2; i++){
            File << pwr_spec_binned[i][j] << "\t";
        }
        File << "\n";
	}
}

void print_vel_pow_spec(const Data_Vec<FTYPE_t, 2> &pwr_spec_binned, string out_dir, string suffix)
{
	out_dir += "vel_pwr_spec/";
	string file_name = out_dir + "vel_pwr_spec" + suffix + ".dat";
	Ofstream File(file_name);
	
	cout << "Writing velocity divergence power spectrum into file " << file_name << "\n";
	File << "# This file contains velocity divergence power spectrum P(k) in units [(Mpc/h)^3] depending on wavenumber k in units [h/Mpc].\n"
	        "# k [h/Mpc]\tP(k) [(Mpc/h)^3]\n";

    File << scientific;
    const unsigned size = pwr_spec_binned.size();
	for (unsigned j = 0; j < size; j++){
        for (unsigned i = 0; i < 2; i++){
            File << pwr_spec_binned[i][j] << "\t";
        }
        File << "\n";
	}
}

void print_corr_func(const Data_Vec<FTYPE_t, 2> &pwr_spec_binned, string out_dir, string suffix)
{
	out_dir += "corr_func/";
	string file_name = out_dir + "corr_func" + suffix + ".dat";
	Ofstream File(file_name);
	
	cout << "Writing correlation function into file " << file_name << "\n";
	File << "# This file contains correlation function depending on distance r in units [Mpc/h].\n"
	        "# x [Mpc/h]\txsi(r)\n";
    
    const unsigned N = pwr_spec_binned.size();
	for (unsigned j = 0; j < N; j++){
		if (pwr_spec_binned[1][j]) File << pwr_spec_binned[0][j] << "\t" << pwr_spec_binned[1][j] << "\n";
	}
}

template<typename T>
T rel_error(const T& a, const T& b)
{
    return a ? fabs((a-b)/a) : fabs(a-b);
}

template<typename T>
bool is_err(const vector<T>& vec1, const vector<T>& vec2, unsigned bin)
{
    const T err = rel_error( vec1[bin], vec2[bin]);
    constexpr T prec_err = is_same<T, float>::value ? 1e-3f : 1e-7;
    constexpr T prec_war = is_same<T, float>::value ? 1e-5f : 1e-12;

    if (err > prec_err){
        cout << "ERROR! Different values of k in bin " << bin << "! Relative error = " << err << "\n";
        return true;
    }
    else if (err > prec_war) cout << "WARNING! Different values of k in bin " << bin << "! Relative error = " << err << "\n";
    return false;
}

void print_pow_spec_diff(const Data_Vec<FTYPE_t, 2> &pwr_spec_binned, const Data_Vec<FTYPE_t, 2> &pwr_spec_binned_0,
	FTYPE_t growth, string out_dir, string suffix)
{
    out_dir += "pwr_diff/";
    string file_name = out_dir + "pwr_spec_diff" + suffix + ".dat";
    Ofstream File(file_name);
    
	cout << "Writing power spectrum difference into file " << file_name << "\n";
    File << "# This file contains relative difference between power spectrum P(k)\n"
            "# and lineary extrapolated power spectrum of initial position of particles\n"
            "# depending on wavenumber k in units [h/Mpc].\n"
	        "# k [h/Mpc]\t(P(k, z)-P_lin(k, z))/P_lin(k, z)\n";

	FTYPE_t P_k, P_lin;
    cout.precision(10);
    const unsigned size = pwr_spec_binned.size();
	for (unsigned j = 0; j < size; j++){
        if (is_err(pwr_spec_binned[0], pwr_spec_binned_0[0], j)) continue;
        P_k = pwr_spec_binned[1][j];
        P_lin = pwr_spec_binned_0[1][j] * pow2(growth);
        File << scientific << pwr_spec_binned_0[0][j] << "\t" << fixed << (P_k-P_lin)/P_lin << "\n";
	}
}

void print_pow_spec_diff(const Data_Vec<FTYPE_t, 2> &pwr_spec_binned, const Interp_obj &pwr_spec_input,
	FTYPE_t growth, string out_dir, string suffix)
{
    out_dir += "pwr_diff/";
    string file_name = out_dir + "pwr_spec_diff" + suffix + ".dat";
    Ofstream File(file_name);
    
	cout << "Writing power spectrum difference into file " << file_name << "\n";
    File << "# This file contains relative difference between power spectrum P(k)\n"
            "# and lineary extrapolated input power spectrum\n"
            "# depending on wavenumber k in units [h/Mpc].\n"
	        "# k [h/Mpc]\t(P(k, z)-P_lin(k, z))/P_lin(k, z)\n";

	FTYPE_t k, P_k, P_lin;
    const unsigned size = pwr_spec_binned.size();
	for (unsigned j = 0; j < size; j++){
        k = pwr_spec_binned[0][j];
        if (k < pwr_spec_input.x_min) continue;
        else if(k > pwr_spec_input.x_max) break;
        else
        {
            P_k = pwr_spec_binned[1][j];
            P_lin = pwr_spec_input(k) * pow2(growth);
            File << scientific << k << "\t" << fixed << (P_k-P_lin)/P_lin << "\n";
        }
	}
}

void print_pow_spec_diff(const Data_Vec<FTYPE_t, 2> &pwr_spec_binned, const Data_Vec<FTYPE_t, 2> &pwr_spec_binned_0,
    const Interp_obj &pwr_spec_input, FTYPE_t growth_now, FTYPE_t growth_init, string out_dir, string suffix)
{
    out_dir += "pwr_diff/";
    string file_name = out_dir + "pwr_spec_diff" + suffix + ".dat";
    Ofstream File(file_name);
    
	cout << "Writing power spectrum difference into file " << file_name << "\n";
    File << "# This file contains relative difference between power spectrum P(k)\n"
            "# and lineary extrapolated 'hybrid' power spectrum\n"
            "# depending on wavenumber k in units [h/Mpc].\n"
	        "# k [h/Mpc]\t(P(k, z)-P_lin(k, z))/P_lin(k, z)\n";

	FTYPE_t k, P_k, P_input, P_par;
    cout.precision(10);
    const unsigned size = pwr_spec_binned.size();
	for (unsigned j = 0; j < size; j++){
        if (is_err(pwr_spec_binned[0], pwr_spec_binned_0[0], j)) continue;
        k = pwr_spec_binned_0[0][j];
        if (k < pwr_spec_input.x_min) continue;
        else if(k > pwr_spec_input.x_max) break;
        else
        {
            P_k = pwr_spec_binned[1][j];
            P_input = pwr_spec_input(k) * pow2(growth_now);
            P_par = pwr_spec_binned_0[1][j] * pow2(growth_now / growth_init);
            File << scientific << k << "\t" << fixed << P_k/P_input - P_k/P_par << "\n";
        }
	}
}

void print_vel_pow_spec_diff(const Data_Vec<FTYPE_t, 2> &pwr_spec_binned, const Data_Vec<FTYPE_t, 2> &pwr_spec_binned_0,
	FTYPE_t dDda, string out_dir, string suffix)
{
    out_dir += "vel_pwr_diff/";
    string file_name = out_dir + "vel_pwr_spec_diff" + suffix + ".dat";
    Ofstream File(file_name);
    
	cout << "Writing power velocity divergence spectrum difference into file " << file_name << "\n";
    File << "# This file contains relative difference between velocity divergence power spectrum P(k)\n"
            "# and lineary extrapolated velocity divergence power spectrum depending on wavenumber k in units [h/Mpc].\n"
	        "# k [h/Mpc]\t(P(k, z)-P_lin(k, z))/P_lin(k, z)\n";
	
	FTYPE_t P_k, P_ZA;
    cout.precision(10);
    const unsigned size = pwr_spec_binned.size();
	for (unsigned j = 0; j < size; j++){
        if (is_err(pwr_spec_binned[0], pwr_spec_binned_0[0], j)) continue;
        P_k = pwr_spec_binned[1][j];
        P_ZA = pwr_spec_binned_0[1][j] * pow2(dDda);
        File << scientific << pwr_spec_binned_0[0][j] << "\t" << fixed << (P_k-P_ZA)/P_ZA << "\n";
	}
}

void print_rho_map(const Mesh& delta, const Sim_Param &sim, string out_dir, string suffix)
{
    out_dir += "rho_map/";
    const FTYPE_t x_0 = sim.x_0_pwr();
    string file_name = out_dir + "rho_map" + suffix + ".dat";
    Ofstream File(file_name);

	cout << "Writing density map into file " << file_name << "\n";
	File << "# This file contains density map delta(x).\n";
    File << "# x [Mpc/h]\tz [Mpc/h]\tdelta\n";
    const unsigned N = sim.box_opt.mesh_num_pwr;
	for (unsigned i = 0; i < N; i++){
		for (unsigned j = 0; j < N; j++){
            File << i*x_0 << "\t" << j*x_0 << "\t" << delta(i, sim.box_opt.mesh_num_pwr/2, j) << "\n";
		}
		File << "\n";
	}
}

void print_projected_rho(const Mesh& delta, const Sim_Param &sim, string out_dir, string suffix)
{
    out_dir += "rho_map/";
    const FTYPE_t x_0 = sim.x_0_pwr();
    string file_name = out_dir + "rho_map_projected" + suffix + ".dat";
    Ofstream File(file_name);
    
	cout << "Writing density map into file " << file_name << "\n";
	File << "# This file contains density map delta(x).\n"
	        "# x [Mpc/h]\tz [Mpc/h]\tdelta\n";
    FTYPE_t rho, rho_tmp;
    const unsigned N = sim.box_opt.mesh_num_pwr;
	for (unsigned i = 0; i < N; i++){
		for (unsigned j = 0; j < N; j++){
			rho = 0;
			for (unsigned k = 0; k < N; k++){
				rho_tmp = delta(i, k, j);
			//	if (rho_tmp != -1) printf("Density in (%i, %i, %i) = %f\n", i, j, k, rho_tmp);
				rho+=rho_tmp + 1;
            }
            File << i*x_0 << "\t" << j*x_0 << "\t" << rho -1 << "\n";
		}
		File << "\n";
	}
}

void print_dens_bin(const vector<int> &dens_binned, string out_dir, string suffix){
    out_dir += "rho_bin/";
    string file_name = out_dir + "rho_bin" + suffix + ".dat";
    Ofstream File(file_name);
    
	cout << "Writing binned density into file " << file_name << "\n";
	File << "# This file contains binned density field.\n"
	        "# dens\tbin_num\n";
    FTYPE_t dens;
    const unsigned N = dens_binned.size();
	for (unsigned j = 0; j < N; j++)
	{
        dens = j*0.1-0.9;
        File << dens << "\t" << dens_binned[j] << "\n";       
	}
}

template void print_par_pos_cut_small(const std::vector<Particle_x<FTYPE_t>>&, const Sim_Param&, string, string);
template void print_par_pos_cut_small(const std::vector<Particle_v<FTYPE_t>>&, const Sim_Param&, string, string);