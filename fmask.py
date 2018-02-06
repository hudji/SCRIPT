"""
Implement the cloud and shadow algorithms known collectively as Fmask,
as published in

Zhu, Z. and Woodcock, C.E. (2012).
Object-based cloud and cloud shadow detection in Landsat imagery
Remote Sensing of Environment 118 (2012) 83-94.

and

Zhu, Z., Wang, S. and Woodcock, C.E. (2015).
Improvement and expansion of the Fmask algorithm: cloud, cloud
shadow, and snow detection for Landsats 4-7, 8, and Sentinel 2 images
Remote Sensing of Environment 159 (2015) 269-277.

Taken from Neil Flood's implementation by permission.

The notation and variable names are largely taken from the paper. Equation
numbers are also from the paper.

Input is a top of atmosphere (TOA) reflectance file.

The output file is a single thematic raster layer with codes representing
null, clear, cloud, shadow, snow and water. These are the values 0-5
respectively, but there are constants defined for the different codes,
as fmask.fmask.OUTCODE_*

"""
# This file is part of 'python-fmask' - a cloud masking module
# Copyright (C) 2015  Neil Flood
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
from __future__ import print_function, division

import sys
import os
import subprocess
import tempfile

import numpy

numpy.seterr(all='raise')
from osgeo import gdal

gdal.UseExceptions()
from scipy.ndimage import uniform_filter, maximum_filter, label
import scipy.stats

# We use RIOS intensively here
from rios import applier
from rios import pixelgrid
from rios import rat
from rios import imageio
from rios import fileinfo

# our wrappers for bits of C that are installed with this package
from . import fillminima
from . import valueindexes
# configuration classes
from . import config
# exceptions
from . import fmaskerrors
# so we can check if thermal all zeroes
from . import zerocheck

# Bands in the saturation mask, if supplied
SATURATION_BLUE = 0
SATURATION_GREEN = 1
SATURATION_RED = 2
