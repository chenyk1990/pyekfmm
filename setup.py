#!/usr/bin/env python
# -*- encoding: utf8 -*-
import glob
import inspect
import io
import os

from setuptools import find_packages
from setuptools import setup


long_description = """
Source code: https://github.com/chenyk1990/pyekfmm""".strip() 

def read(*names, **kwargs):
    return io.open(
        os.path.join(os.path.dirname(__file__), *names),
        encoding=kwargs.get("encoding", "utf8")).read()


from distutils.core import Extension
import numpy
eikonalc_module = Extension('eikonalc', sources=['pyekfmm/src/eikonal.c'], 
										include_dirs=[numpy.get_include()])

eikonalvtic_module = Extension('eikonalvtic', sources=['pyekfmm/src/eikonalvti.c'], 
										include_dirs=[numpy.get_include()])
# from numpy.distutils.core import setup 
# from distutils.core import setup
setup(
    name="pyekfmm",
    version="0.0.8.3",
    license='GNU General Public License, Version 3 (GPLv3)',
    description="Fast Marching Method for Traveltime Calculation",
    long_description=long_description,
    author="pyekfmm developing team",
    author_email="chenyk2016@gmail.com",
    url="https://github.com/chenyk1990/pyekfmm",
    ext_modules=[eikonalc_module,eikonalvtic_module],
    packages=['pyekfmm'],
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Operating System :: Unix",
        "Operating System :: POSIX",
        "Operating System :: Microsoft :: Windows",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: Implementation :: CPython",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)"
    ],
    keywords=[
        "seismology", "exploration seismology", "array seismology", "denoising", "science", "signal-to-noise ratio", "damped rank reduction method"
    ],
    install_requires=[
        "numpy", "scipy", "matplotlib"
    ],
    extras_require={
        "docs": ["sphinx", "ipython", "runipy"]
    }
    
)
