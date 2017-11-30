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

// Individual headers
%include "core.i"
%include "emu.i"
%include "power.i"
%include "ccl.i"
