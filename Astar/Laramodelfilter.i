%module LaraModelfilter
%{
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"
%init %{
import_array();
%}


%{
#include "model.h"
#include "LaraModelfilter.h"
%}

%include "model.h"
%include "LaraModelfilter.h"
