// SWIG interface
%module fastsim
%include "std_string.i"

%{
#include "core.h"
%}

%rename(__assign__) *::operator=;
%rename(__getitem__) *::operator[];

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

// Parse original header files
%ignore Mesh_base;
%include "core.h"

// Instantiate templates
 
%template(App_x) App_Var<Particle_x>;
%template(App_v) App_Var<Particle_v>;