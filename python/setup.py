from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy as np
import os

os.environ["CC"] = "mpicxx"
setup(
        ext_modules = cythonize([
            Extension(
                "python.generic_io",
                sources=["./python/wrapper.pyx"],
                language="c++",
                extra_link_args=["-fopenmp"],
                extra_objects=[
                    "./mpi/GenericIO.o",
                    "./mpi/thirdparty/blosc/combined_blosc.o",
                ],
            ),
        ]),
        include_dirs = [
                np.get_include(),
        ]
)
