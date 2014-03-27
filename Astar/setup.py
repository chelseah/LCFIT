#! /usr/bin/env python

from distutils.core import *
from distutils import sysconfig
from distutils.extension import Extension

import numpy

numpy_include = numpy.get_include()
_elliptic = Extension("_elliptic",
        ["elliptic_wrap.cxx",
            "elliptic.cc"],
        include_dirs = [numpy_include],
        )
_quadpolar = Extension("_quadpolar",
        ["quadpolar_wrap.cxx",
            "quadpolar.cc"],
        include_dirs = [numpy_include],
        )
_quad2d = Extension("_quad2d",
        ["quad2d_wrap.cxx",
            "quad2d.cc"],
        include_dirs = [numpy_include],
        )


_quadsphere = Extension("_quadsphere",
        ["quadsphere_wrap.cxx",
            "quadsphere.cc"],
        include_dirs = [numpy_include],
        )

_ZeipelModel = Extension("_ZeipelModel",
        ["ZeipelModel_wrap.cxx",
            "ZeipelModel.cc"],
        include_dirs = [numpy_include],
        )

_Zeipel = Extension("_Zeipel",
        ["Zeipel_wrap.cxx",
            "Zeipel.cc"],
        include_dirs = [numpy_include],
        )

_LaraModel = Extension("_LaraModel",
        ["LaraModel_wrap.cxx",
            "LaraModel.cc"],
        include_dirs = [numpy_include],
        )

_Lara = Extension("_Lara",
        ["Lara_wrap.cxx",
            "Lara.cc"],
        include_dirs = [numpy_include],
        )

_LaraModelfilter = Extension("_LaraModelfilter",
        ["LaraModelfilter_wrap.cxx",
            "LaraModelfilter.cc"],
        include_dirs = [numpy_include],
        )

_Larafilter = Extension("_Larafilter",
        ["Larafilter_wrap.cxx",
            "Larafilter.cc"],
        include_dirs = [numpy_include],
        )


#_Zeipel_all = Extension("_Zeipel_all",
#        ["Zeipel_all_wrap.cxx",
#            "Zeipel_all.cc"],
#        include_dirs = [numpy_include],
#        )



setup(name="Oblate",
        description = "model for an oblate star",
        author = "Xu Huang",
        author_email = "",
        url = "",
        version = "0.0.0",
        py_modules = ["elliptic","quadpolar","quad2d","quadsphere","ZeipelModel","Zeipel","LaraModel","Lara","LaraModelfilter","Larafilter"],
        ext_modules = [_elliptic,_quadpolar,_quad2d,_quadsphere,_ZeipelModel,_Zeipel,_LaraModel,_Lara,_LaraModelfilter,_Larafilter])
