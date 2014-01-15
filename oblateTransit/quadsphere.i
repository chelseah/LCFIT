%module quadsphere
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
#include "quadsphere.h"
%}

%include "model.h"
%include "integrator.h"
%include "quadsphere.h"
