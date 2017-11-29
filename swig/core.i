// C++ code to wrap
%{
#include "core.h"
%}

%rename(__assign__) *::operator=;
%rename(__getitem__) *::operator[];

%ignore Mesh_base::Mesh_base(Mesh_base&& that);

template <typename T>
class Mesh_base
{
public:
	// CONSTRUCTORS & DESTRUCTOR
	Mesh_base(unsigned n1, unsigned n2, unsigned n3);
    Mesh_base(const Mesh_base& that);
    Mesh_base(Mesh_base<T>&& that) noexcept;
    Mesh_base& operator=(Mesh_base that);
    template <class U> friend void swap(Mesh_base<U>& first, Mesh_base<U>& second);
	~Mesh_base();
	
	// VARIABLES
	unsigned N1, N2, N3, length; // acces dimensions and length of mesh
	
	// METHODS
	T* real() const { return data;} // acces data
	void assign(T val);
	
	// OPERATORS
	T &operator[](int i){ return data[i]; }
	const T &operator[](int i) const{ return data[i]; }
	
	T& operator()(int i, int j, int k){ return data[i*N2*N3+j*N3+k]; }
	const T& operator()(int i, int j, int k) const{ return data[i*N2*N3+j*N3+k]; }
	
	T& operator()(int i, int j){ return data[i*N3+j]; }
	const T& operator()(int i, int j) const{ return data[i*N3+j]; }
	
	T& operator()(Vec_3D<int> pos);
	const T& operator()(Vec_3D<int> pos) const;
	
	Mesh_base& operator+=(const T& rhs);
	Mesh_base& operator-=(const T& rhs){ return *this+=-rhs; }
	Mesh_base& operator*=(const T& rhs);
	Mesh_base& operator/=(const T& rhs);

protected:
	// VARIABLES
	T* data;
};
%template(Mesh_base_d) Mesh_base<double>;

template <class T> 
class App_Var
{
public:
	// CONSTRUCTORS
    App_Var(const Sim_Param &sim, std::string app_str);

	// DESTRUCTOR
	~App_Var();
	
    // VARIABLES
    const Sim_Param &sim;

	int err, step, print_every;
    double b, b_out, db;
    double D_init, dDda_init;
    const std::string app_str, z_suffix_const, out_dir_app;
    
    // LARGE FIELDS
	std::vector<Mesh> app_field;
    std::vector<Mesh> power_aux;
    T* particles;

    // OTHER VARIABLES
    Data_Vec<double, 3> pwr_spec_binned, pwr_spec_binned_0, vel_pwr_spec_binned_0;
    Data_Vec<double, 2> corr_func_binned;
    Interp_obj pwr_spec_input;
	fftw_plan p_F, p_B, p_F_pwr, p_B_pwr;
	Tracking track;
	std::vector<int> dens_binned;
	
	// METHODS
	double z() const{ return 1./b - 1.;}
	double b_half() const { return b - db/2.; }
	bool integrate() const { return (b <= b_out) && (db > 0);}
	bool printing() const { return print_every ? ((step % print_every) == 0) or (b == b_out) : print_every ; }
    void print_output();
    void upd_time();
    void print_mem() const;
    void print_info() const;	
	std::string z_suffix();
	bool is_init_pwr_spec_0, is_init_vel_pwr_spec_0; //< be careful about setting these to true
protected:	
    std::stringstream z_suffix_num;
    uint64_t memory_alloc; // only the largest chunks
};
%template(App_Var_x) App_Var<Particle_x>;
%template(App_Var_v) App_Var<Particle_v>;

// Parse original header file
%ignore Mesh_base;
%ignore App_Var;
%include "core.h"

// Instantiate templates
%template(App_x) App_Var<Particle_x>;
%template(App_v) App_Var<Particle_v>;
%template(Vec_d) std::vector<double>;
%template(Data_d2) Data_Vec<double, 2>;