%module Zeipel_all
%{
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"
%init %{
import_array();
%}

%apply (double * INPLACE_ARRAY1, int DIM1){(double *phase, int np), (double *F, int nf)};

%{
#include "Zeipel_all.h"
%}

%include "Zeipel_all.h"
