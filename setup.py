#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension

setup(name="heat",
      ext_modules=[
          Extension("heat",
                    sources=["heatmodule.cpp"],
                    libraries = ["boost_python"])
    ])
