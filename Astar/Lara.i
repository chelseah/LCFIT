%module Lara
%{
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"
%init %{
import_array();
%}

%apply (double * INPLACE_ARRAY1, int DIM1){(double *x, int nx), (double *y, int ny),(double *Txy, int nT)};
%apply (double * INPLACE_ARRAY1, int DIM1){(double *phase, int np), (double *F, int nf)};
%apply (double * INPLACE_ARRAY1, int DIM1){(double *F0, int np)};
%apply (double * INPLACE_ARRAY1, int DIM1){(double *theta, int nth), (double *T, int nT)};

%{
#include "model.h"
#include "integrator.h"
#include "quadpolar.h"
#include "quadsphere.h"
#include "Laramodel.h"
#include "Lara.h"
%}

%include "model.h"
%include "integrator.h"
%include "quadpolar.h"
%include "quadsphere.h"
#include "Laramodel.h"
%include "Lara.h"
