// C++ code to wrap
%{
#include "core_power.h"
%}

// Parse original header file
%include "core_power.h"

// Instantiate templates
TEMP_CLASS_GEN(Extrap_Pk)
TEMP_CLASS_GEN(Extrap_Pk_Nl)
TEMP_FUNC_GEN(gen_corr_func_binned_gsl_qawf, Extrap_Pk)
TEMP_FUNC_GEN(gen_corr_func_binned_gsl_qawf, Extrap_Pk_Nl)
TEMP_FUNC_GEN(gen_sigma_binned_gsl_qawf, Extrap_Pk)
TEMP_FUNC_GEN(gen_sigma_binned_gsl_qawf, Extrap_Pk_Nl)
