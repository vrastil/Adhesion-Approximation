// C++ code to wrap
%{
#include "core_power.h"
%}

// make Extrap_Pk callable from python
%rename(__call__) Extrap_Pk::eval;

// Parse original header file
%include "core_power.h"

// Instantiate templates
%extend Extrap_Pk {
    %template(Extrap_Pk) Extrap_Pk<2>;
    %template(Extrap_Pk) Extrap_Pk<3>;
}
%template(gen_corr_func_binned_gsl_qawf) gen_corr_func_binned_gsl_qawf<Extrap_Pk>;
%template(gen_corr_func_binned_gsl_qawf) gen_corr_func_binned_gsl_qawf<Extrap_Pk_Nl>;