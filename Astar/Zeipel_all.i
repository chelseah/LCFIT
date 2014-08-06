%module Zeipel_all
%{
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"
%init %{
import_array();
%}

%apply (double * INPLACE_ARRAY1, int DIM1){(double *phase, int np), (double *F, int nf)};
%apply (double * INPLACE_ARRAY1, int DIM1){(double *F0, int np)};

%{
#include "Zeipel_all.h"
%}

%include "Zeipel_all.h"
