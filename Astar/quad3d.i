%module quad3d
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
#include "quad3d.h"
%}

%include "model.h"
%include "integrator.h"
%include "quad3d.h"
