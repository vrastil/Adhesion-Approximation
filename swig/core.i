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
TEMP_CLASS_GEN(Data_Vec)