#! /usr/bin/env python

from distutils.core import *
from distutils import sysconfig
from distutils.extension import Extension

import numpy

numpy_include = numpy.get_include()
_oblateness = Extension("_oblateness",
        ["oblateness_wrap.cxx",
            "oblateness.cc"],
        include_dirs = [numpy_include],
        )

_elliptic = Extension("_elliptic",
        ["elliptic_wrap.cxx",
            "elliptic.cc"],
        include_dirs = [numpy_include],
        )
_oblatenessfast = Extension("_oblatenessfast",
        ["oblatenessfast_wrap.cxx",
            "oblatenessfast.cc"],
        include_dirs = [numpy_include],
        )

_quadpolar = Extension("_quadpolar",
        ["quadpolar_wrap.cxx",
            "quadpolar.cc"],
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


setup(name="Oblate",
        description = "model for an oblate planet and a oblate star",
        author = "Wei Zhu, Xu Huang",
        author_email = "",
        url = "",
        version = "0.0.0",
        py_modules = ["oblateness","elliptic","oblatenessfast","quadpolar","quadsphere","ZeipelModel","Zeipel"],
        ext_modules = [_oblateness,_elliptic,_oblatenessfast,_quadpolar,_quadsphere,_ZeipelModel,_Zeipel])
