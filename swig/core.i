// C++ code to wrap
%{
#include "precision.hpp"
#include "class_data_vec.hpp"
#include "params.hpp"
%}

%rename(__assign__) *::operator=;
%rename(__getitem__) *::operator[];

// Parse original header file
%ignore Cosmo_Param::operator void*;
%include "precision.hpp"
%include "class_data_vec.hpp"
%include "params.hpp"

// Instantiate templates
TEMP_CLASS_GEN(Data_Vec)