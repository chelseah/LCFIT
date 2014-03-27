%module Larafilter
%{
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"
%init %{
import_array();
%}

%apply (double * INPLACE_ARRAY1, int DIM1){(double *filter, int nf)};

%{
#include "model.h"
#include "integrator.h"
#include "quadpolar.h"
#include "quad2d.h"
#include "Laramodelfilter.h"
#include "Larafilter.h"
%}

%include "model.h"
%include "integrator.h"
%include "quadpolar.h"
%include "quad2d.h"
%include "Laramodelfilter.h"
%include "Larafilter.h"

%clear (double *filter, int nf);
