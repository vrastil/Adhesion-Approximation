// SWIG interface
// module name
%module fastsim

// C++ STL
%include "std_string.i"
%include "std_vector.i"
%include "std_map.i"

// Individual headers
%include "core.i"
%include "emu.i"
%include "power.i"
%include "ccl.i"
