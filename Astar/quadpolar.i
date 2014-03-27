%module quadpolar
%{
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"
%init %{
import_array();
%}


%{
#include "model.h"
#include "integrator.h"
#include "quadpolar.h"
%}

%include "model.h"
%include "integrator.h"
%include "quadpolar.h"
