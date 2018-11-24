/**
 * @brief functions handling output of the program
 * 
 * @file core_out.cpp
 * @author Michal Vrastil
 * @date 2018-07-11
 */

#include <boost/filesystem.hpp>
#include "core_out.h"
#include "class_particles.hpp"

namespace fs = boost::filesystem;

#define BUFF_SIZE 1024 * 1024 * 16 // 16 MB buffer

Ofstream::Ofstream(std::string file_name) : std::ofstream(file_name), buf(new char[BUFF_SIZE])
{
    if (!this->is_open()) throw std::runtime_error("Error while opening '" + file_name + "'");
    this->rdbuf()->pubsetbuf(buf, sizeof(buf));
}

Ofstream::~Ofstream()
{
    delete[] buf;
    if (this->is_open()) this->close();
}

Ifstream::Ifstream(std::string file_name) : std::ifstream(file_name)
{
    if (!this->is_open()) throw std::runtime_error("Error while opening '" + file_name + "'");
}

Ifstream::~Ifstream()
{
    if (this->is_open()) this->close();
}

// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
std::string currentDateTime()
{
	const time_t     now = time(0);
	struct tm  tstruct;
	char       buf[80];
	gmtime_r(&now, &tstruct);
	strftime(buf, sizeof(buf), "%y%m%d_%H%M%S", &tstruct);
	
	std::string returnval(buf);
    return returnval;
}

std::string std_out_dir(const std::string& pre_subdir, const Sim_Param &sim)
{
    /// directory name = YYMMDD_HHMMSS_m_p_M_b
    const std::string out_dir = sim.out_opt.out_dir + pre_subdir + currentDateTime() + "_" + std::to_string(sim.box_opt.mesh_num) +"m_" +
                      std::to_string(sim.box_opt.Ng) + "p_" + std::to_string(sim.box_opt.mesh_num_pwr) +"M_" + 
                      std::to_string((size_t)sim.box_opt.box_size) + "b";

    /// check if directory exists
    if (!fs::exists(fs::path(out_dir.c_str()))) return out_dir + "/";
    
    /// try appending numbers starting at 2
    else
    {
        for(size_t i = 2; ; ++i)
        {
            const std::string out_dir_new = out_dir  + "_" + std::to_string(i);
            if (!fs::exists(fs::path(out_dir_new.c_str()))) return out_dir_new + "/";
        }
    }
}

void create_dir(const std::string &out_dir)
{
	const fs::path dir(out_dir.c_str());
	if(fs::create_directories(dir))
    {
        std::cout << "Directory created: "<< out_dir << "\n";
    }
}

void remove_dir(const std::string &out_dir)
{
    const fs::path dir(out_dir.c_str());
    if (fs::remove_all(dir))
    {
        std::cout << "Directory removed: "<< out_dir << "\n";
    }
}

void remove_all_files(const std::string &out_dir)
{
    const fs::path dir(out_dir.c_str());
    size_t i = 0;

    for(auto & p : fs::directory_iterator(dir))
    {
        if (fs::is_regular_file(p))
        {
            fs::remove(p);
            ++i;
        }
    }
    std::cout << "Removed " << i << " file(s) in directory: "<< out_dir << "\n";
}

template <class T>
void print_par_pos_cut_small(const std::vector<T>& particles, const Sim_Param &sim, std::string out_dir, std::string suffix)
{
   out_dir += "par_cut/";
   std::string file_name = out_dir + "par_cut" + suffix + ".dat";
   Ofstream File(file_name);
   
   std::cout << "Writing small cut through the box of particles into file " << file_name << "\n";
   File << "# This file contains positions of particles in units [Mpc/h].\n";
   FTYPE_t x, y, z, dx;
   const FTYPE_t x_0 = sim.x_0();
   const size_t N = sim.box_opt.par_num;
   for(size_t i=0; i < N; i++)
   {
       x = particles[i].position[0];
       y = particles[i].position[1];
       z = particles[i].position[2];			
       dx = std::abs(y - sim.box_opt.mesh_num/2);
       if ((dx < 0.5) && (x < sim.box_opt.mesh_num/4.) && (z < sim.box_opt.mesh_num/4.))
       {
           // cut (L/4 x L/4 x 0.5)
           File << x*x_0 << "\t" << z*x_0 << "\t" << y*x_0 << "\n";
       }
   }
}

void print_pow_spec(const Data_Vec<FTYPE_t, 2> &pwr_spec_binned, std::string out_dir, std::string suffix)
{
	out_dir += "pwr_spec/";
	std::string file_name = out_dir + "pwr_spec" + suffix + ".dat";
	Ofstream File(file_name);
	
	std::cout << "Writing power spectrum into file " << file_name << "\n";
	File << "# This file contains power spectrum P(k) in units [(Mpc/h)^3] depending on wavenumber k in units [h/Mpc].\n"
	        "# k [h/Mpc]\tP(k) [(Mpc/h)^3]\n";

    File << std::scientific;
    const size_t size = pwr_spec_binned.size();
	for (size_t j = 0; j < size; j++){
        for (size_t i = 0; i < 2; i++){
            File << pwr_spec_binned[i][j] << "\t";
        }
        File << "\n";
	}
}

void print_vel_pow_spec(const Data_Vec<FTYPE_t, 2> &pwr_spec_binned, std::string out_dir, std::string suffix)
{
	out_dir += "vel_pwr_spec/";
	std::string file_name = out_dir + "vel_pwr_spec" + suffix + ".dat";
	Ofstream File(file_name);
	
	std::cout << "Writing velocity divergence power spectrum into file " << file_name << "\n";
	File << "# This file contains velocity divergence power spectrum P(k) in units [(Mpc/h)^3] depending on wavenumber k in units [h/Mpc].\n"
	        "# k [h/Mpc]\tP(k) [(Mpc/h)^3]\n";

    File << std::scientific;
    const size_t size = pwr_spec_binned.size();
	for (size_t j = 0; j < size; j++){
        for (size_t i = 0; i < 2; i++){
            File << pwr_spec_binned[i][j] << "\t";
        }
        File << "\n";
	}
}

void print_corr_func(const Data_Vec<FTYPE_t, 2> &pwr_spec_binned, std::string out_dir, std::string suffix)
{
	out_dir += "corr_func/";
	std::string file_name = out_dir + "corr_func" + suffix + ".dat";
	Ofstream File(file_name);
	
	std::cout << "Writing correlation function into file " << file_name << "\n";
	File << "# This file contains correlation function depending on distance r in units [Mpc/h].\n"
	        "# x [Mpc/h]\txsi(r)\n";
    
    const size_t N = pwr_spec_binned.size();
	for (size_t j = 0; j < N; j++){
		if (pwr_spec_binned[1][j]) File << pwr_spec_binned[0][j] << "\t" << pwr_spec_binned[1][j] << "\n";
	}
}

template<typename T>
T rel_error(const T& a, const T& b)
{
    return a ? fabs((a-b)/a) : fabs(a-b);
}

template<typename T>
bool is_err(const std::vector<T>& vec1, const std::vector<T>& vec2, size_t bin)
{
    const T err = rel_error( vec1[bin], vec2[bin]);
    constexpr T prec_err = std::is_same<T, float>::value ? 1e-3f : 1e-7;
    constexpr T prec_war = std::is_same<T, float>::value ? 1e-5f : 1e-12;

    if (err > prec_err){
        std::cout << "ERROR! Different values of k in bin " << bin << "! Relative error = " << err << "\n";
        return true;
    }
    else if (err > prec_war) std::cout << "WARNING! Different values of k in bin " << bin << "! Relative error = " << err << "\n";
    return false;
}

void print_pow_spec_diff(const Data_Vec<FTYPE_t, 2> &pwr_spec_binned, const Data_Vec<FTYPE_t, 2> &pwr_spec_binned_0,
	FTYPE_t growth, std::string out_dir, std::string suffix)
{
    out_dir += "pwr_diff/";
    std::string file_name = out_dir + "pwr_spec_diff" + suffix + ".dat";
    Ofstream File(file_name);
    
	std::cout << "Writing power spectrum difference into file " << file_name << "\n";
    File << "# This file contains relative difference between power spectrum P(k)\n"
            "# and lineary extrapolated power spectrum of initial position of particles\n"
            "# depending on wavenumber k in units [h/Mpc].\n"
	        "# k [h/Mpc]\t(P(k, z)-P_lin(k, z))/P_lin(k, z)\n";

	FTYPE_t P_k, P_lin;
    std::cout.precision(10);
    const size_t size = pwr_spec_binned.size();
	for (size_t j = 0; j < size; j++){
        if (is_err(pwr_spec_binned[0], pwr_spec_binned_0[0], j)) continue;
        P_k = pwr_spec_binned[1][j];
        P_lin = pwr_spec_binned_0[1][j] * pow2(growth);
        File << std::scientific << pwr_spec_binned_0[0][j] << "\t" << std::fixed << (P_k-P_lin)/P_lin << "\n";
	}
}

void print_pow_spec_diff(const Data_Vec<FTYPE_t, 2> &pwr_spec_binned, const Interp_obj &pwr_spec_input,
	FTYPE_t growth, std::string out_dir, std::string suffix)
{
    out_dir += "pwr_diff/";
    std::string file_name = out_dir + "pwr_spec_diff" + suffix + ".dat";
    Ofstream File(file_name);
    
	std::cout << "Writing power spectrum difference into file " << file_name << "\n";
    File << "# This file contains relative difference between power spectrum P(k)\n"
            "# and lineary extrapolated input power spectrum\n"
            "# depending on wavenumber k in units [h/Mpc].\n"
	        "# k [h/Mpc]\t(P(k, z)-P_lin(k, z))/P_lin(k, z)\n";

	FTYPE_t k, P_k, P_lin;
    const size_t size = pwr_spec_binned.size();
	for (size_t j = 0; j < size; j++){
        k = pwr_spec_binned[0][j];
        if (k < pwr_spec_input.x_min) continue;
        else if(k > pwr_spec_input.x_max) break;
        else
        {
            P_k = pwr_spec_binned[1][j];
            P_lin = pwr_spec_input(k) * pow2(growth);
            File << std::scientific << k << "\t" << std::fixed << (P_k-P_lin)/P_lin << "\n";
        }
	}
}

void print_pow_spec_diff(const Data_Vec<FTYPE_t, 2> &pwr_spec_binned, const Data_Vec<FTYPE_t, 2> &pwr_spec_binned_0,
    const Interp_obj &pwr_spec_input, FTYPE_t growth_now, FTYPE_t growth_init, std::string out_dir, std::string suffix)
{
    out_dir += "pwr_diff/";
    std::string file_name = out_dir + "pwr_spec_diff" + suffix + ".dat";
    Ofstream File(file_name);
    
	std::cout << "Writing power spectrum difference into file " << file_name << "\n";
    File << "# This file contains relative difference between power spectrum P(k)\n"
            "# and lineary extrapolated 'hybrid' power spectrum\n"
            "# depending on wavenumber k in units [h/Mpc].\n"
	        "# k [h/Mpc]\t(P(k, z)-P_lin(k, z))/P_lin(k, z)\n";

	FTYPE_t k, P_k, P_input, P_par;
    std::cout.precision(10);
    const size_t size = pwr_spec_binned.size();
	for (size_t j = 0; j < size; j++){
        if (is_err(pwr_spec_binned[0], pwr_spec_binned_0[0], j)) continue;
        k = pwr_spec_binned_0[0][j];
        if (k < pwr_spec_input.x_min) continue;
        else if(k > pwr_spec_input.x_max) break;
        else
        {
            P_k = pwr_spec_binned[1][j];
            P_input = pwr_spec_input(k) * pow2(growth_now);
            P_par = pwr_spec_binned_0[1][j] * pow2(growth_now / growth_init);
            File << std::scientific << k << "\t" << std::fixed << P_k/P_input - P_k/P_par << "\n";
        }
	}
}

void print_vel_pow_spec_diff(const Data_Vec<FTYPE_t, 2> &pwr_spec_binned, const Data_Vec<FTYPE_t, 2> &pwr_spec_binned_0,
	FTYPE_t dDda, std::string out_dir, std::string suffix)
{
    out_dir += "vel_pwr_diff/";
    std::string file_name = out_dir + "vel_pwr_spec_diff" + suffix + ".dat";
    Ofstream File(file_name);
    
	std::cout << "Writing power velocity divergence spectrum difference into file " << file_name << "\n";
    File << "# This file contains relative difference between velocity divergence power spectrum P(k)\n"
            "# and lineary extrapolated velocity divergence power spectrum depending on wavenumber k in units [h/Mpc].\n"
	        "# k [h/Mpc]\t(P(k, z)-P_lin(k, z))/P_lin(k, z)\n";
	
	FTYPE_t P_k, P_ZA;
    std::cout.precision(10);
    const size_t size = pwr_spec_binned.size();
	for (size_t j = 0; j < size; j++){
        if (is_err(pwr_spec_binned[0], pwr_spec_binned_0[0], j)) continue;
        P_k = pwr_spec_binned[1][j];
        P_ZA = pwr_spec_binned_0[1][j] * pow2(dDda);
        File << std::scientific << pwr_spec_binned_0[0][j] << "\t" << std::fixed << (P_k-P_ZA)/P_ZA << "\n";
	}
}

void print_rho_map(const Mesh& delta, const Sim_Param &sim, std::string out_dir, std::string suffix)
{
    out_dir += "rho_map/";
    const FTYPE_t x_0 = sim.x_0_pwr();
    std::string file_name = out_dir + "rho_map" + suffix + ".dat";
    Ofstream File(file_name);

	std::cout << "Writing density map into file " << file_name << "\n";
	File << "# This file contains density map delta(x).\n";
    File << "# x [Mpc/h]\tz [Mpc/h]\tdelta\n";
    const size_t N = sim.box_opt.mesh_num_pwr;
	for (size_t i = 0; i < N; i++){
		for (size_t j = 0; j < N; j++){
            File << i*x_0 << "\t" << j*x_0 << "\t" << delta(i, sim.box_opt.mesh_num_pwr/2, j) << "\n";
		}
		File << "\n";
	}
}

void print_projected_rho(const Mesh& delta, const Sim_Param &sim, std::string out_dir, std::string suffix)
{
    out_dir += "rho_map/";
    const FTYPE_t x_0 = sim.x_0_pwr();
    std::string file_name = out_dir + "rho_map_projected" + suffix + ".dat";
    Ofstream File(file_name);
    
	std::cout << "Writing density map into file " << file_name << "\n";
	File << "# This file contains density map delta(x).\n"
	        "# x [Mpc/h]\tz [Mpc/h]\tdelta\n";
    FTYPE_t rho, rho_tmp;
    const size_t N = sim.box_opt.mesh_num_pwr;
	for (size_t i = 0; i < N; i++){
		for (size_t j = 0; j < N; j++){
			rho = 0;
			for (size_t k = 0; k < N; k++){
				rho_tmp = delta(i, k, j);
			//	if (rho_tmp != -1) printf("Density in (%i, %i, %i) = %f\n", i, j, k, rho_tmp);
				rho+=rho_tmp + 1;
            }
            File << i*x_0 << "\t" << j*x_0 << "\t" << rho -1 << "\n";
		}
		File << "\n";
	}
}

void print_dens_bin(const std::vector<size_t> &dens_binned, std::string out_dir, std::string suffix){
    out_dir += "rho_bin/";
    std::string file_name = out_dir + "rho_bin" + suffix + ".dat";
    Ofstream File(file_name);
    
	std::cout << "Writing binned density into file " << file_name << "\n";
	File << "# This file contains binned density field.\n"
	        "# dens\tbin_num\n";
    FTYPE_t dens;
    const size_t N = dens_binned.size();
	for (size_t j = 0; j < N; j++)
	{
        dens = j*0.1-0.9;
        File << dens << "\t" << dens_binned[j] << "\n";       
	}
}

template void print_par_pos_cut_small(const std::vector<Particle_x<FTYPE_t>>&, const Sim_Param&, std::string, std::string);
template void print_par_pos_cut_small(const std::vector<Particle_v<FTYPE_t>>&, const Sim_Param&, std::string, std::string);