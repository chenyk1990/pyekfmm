#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
:copyright:
    Yangkang Chen (chenyk2016@gmail.com), 2022-Present 
:license:
    GNU General Public License, Version 3
    (http://www.gnu.org/copyleft/gpl.html)
"""

__version__ = "0.0.3"

from .fmm import eikonal

from eikonal import *

# This is my first C-extension (more will be coming), the first step to make the codes faster
from calculator import *
from eikonalc import *












