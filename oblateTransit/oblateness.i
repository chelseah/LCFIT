%module oblateness
%{
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"
%init %{
import_array();
%}

%{
#include "oblateness.h"
%}

%include "oblateness.h"
