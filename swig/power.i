// C++ code to wrap
%{
#include "core_power.h"
%}

// Parse original header file
%include "core_power.h"

// Instantiate templates
%extend Extrap_Pk {
    %template(Extrap_Pk) Extrap_Pk<2>;
}