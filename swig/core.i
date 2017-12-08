// C++ code to wrap
%{
#include "core.h"
%}

%rename(__assign__) *::operator=;
%rename(__getitem__) *::operator[];

// Parse original header file
%ignore Cosmo_Param::operator void*;
%include "core.h"

// Instantiate templates
%template(Data_d2) Data_Vec<double, 2>;
%template(Data_d3) Data_Vec<double, 3>;