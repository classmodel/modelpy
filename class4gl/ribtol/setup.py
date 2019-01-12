# build with "python setup.py build_ext --inplace"
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy as np
import os

os.environ["CC"] = "g++-7"

setup(
    ext_modules = cythonize((Extension("ribtol", sources=["ribtol.pyx"], include_dirs=[np.get_include()], ), ))
)
