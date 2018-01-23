// SWIG interface
// module name
%module fastsim

// SWIG STL exception
%include exception.i

%exception {
  try {
    $action
  } catch (const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}

// C++ STL
%include "std_string.i"
%include "std_vector.i"
%include "std_map.i"

%template(vectorf) std::vector<float>;
%template(vectord) std::vector<double>;
%template(vectorl) std::vector<long double>;

// templates generator

%define TEMP_CLASS_GEN(temp_name)
#if PRECISION == 1
%template(temp_name ## _2) temp_name<float, 2>;
%template(temp_name ## _3) temp_name<float, 3>;
#elif PRECISION == 2
%template(temp_name ## _2) temp_name<double, 2>;
%template(temp_name ## _3) temp_name<double, 3>;
#elif PRECISION == 3
%template(temp_name ## _2) temp_name<long double, 2>;
%template(temp_name ## _3) temp_name<long double, 3>;
#endif
%enddef

%define TEMP_FUNC_GEN(func_name, class_name)
#if PRECISION == 1
%template(func_name) func_name<class_name<float, 2> >;
%template(func_name) func_name<class_name<float, 3> >;
#elif PRECISION == 2
%template(func_name) func_name<class_name<double, 2> >;
%template(func_name) func_name<class_name<double, 3> >;
#elif PRECISION == 3
%template(func_name) func_name<class_name<long double, 2> >;
%template(func_name) func_name<class_name<long double, 3> >;
#endif
%enddef

// Individual headers
%include "core.i"
%include "power.i"
%include "ccl.i"
