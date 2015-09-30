'''
compile_test.py - check pyximport
=================================
test script for checking if compilation against
pysam and tabix works.
'''
# clean up previous compilation
import os
try:
    os.unlink('pysam_ext.c')
    os.unlink('pysam_ext.pyxbldc')
except OSError:
    pass

#import pyximport
#pyximport.install(build_in_temp=False)
#import pysam_ext

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import pysam

extensions = Extension(name = "pysam_ext",
                       sources = ["pysam_ext.pyx"],
                       include_dirs = pysam.get_include(),
                       define_macros = pysam.get_defines())

setup(
    name="pysam_ext",
    ext_modules = cythonize(extensions),
)

