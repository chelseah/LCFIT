%module ZeipelModel
%{
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"
%init %{
import_array();
%}


%{
#include "model.h"
#include "ZeipelModel.h"
%}

%include "model.h"
%include "ZeipelModel.h"
