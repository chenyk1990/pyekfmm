#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
:copyright:
    Yangkang Chen (chenyk2016@gmail.com), 2022-Present 
:license:
    GNU General Public License, Version 3
    (http://www.gnu.org/copyleft/gpl.html)
"""

__version__ = "0.0.0"

from .fmm import eikonal


# This is my first C-extension (more will be coming), the first step to make the codes faster
from calculator import *

from eikonal import *

# from spam import *
import npufunc

# import npufunc
# from spam2 import *


# def configuration(parent_package='', top_path=None):
#     from numpy.distutils.misc_util import Configuration
# 
#     config = Configuration('src',
#                            parent_package,
#                            top_path)
#     config.add_extension('npufunc', ['spam2.c'])
# 
#     return config
#     
# if __name__ == '__main__':
# 	print('Compile this one')
# 	from numpy.distutils.core import setup
# 	setup(configuration=configuration)
    









