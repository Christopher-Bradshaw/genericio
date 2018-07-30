from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy as np

# setup(
#         ext_modules = cythonize([
#             Extension(
#                 "wrapper",
#                 sources=["wrapper.pyx"],
#                 language="c++",
#                 extra_link_args=["-fopenmp"],
#                 extra_compile_args=["-DGENERICIO_NO_MPI"],
#                 extra_objects=[
#                     "./frontend/GenericIO.o",
#                     "./frontend/thirdparty/blosc/combined_blosc.o",
#                 ],
#             ),
#         ]),
#         include_dirs = [
#                 np.get_include(),
#         ]
# )

import os
os.environ["CXX"] = "mpicxx"
setup(
        ext_modules = cythonize([
            Extension(
                "mpi_wrapper",
                sources=["wrapper.pyx"],
                language="c++",
                extra_link_args=["-fopenmp", "-L/usr/lib64/openmpi/lib/" "-lmpi"],
                extra_objects=[
                    "./mpi/GenericIO.o",
                    "./mpi/thirdparty/blosc/combined_blosc.o",
                ],
            ),
        ]),
        include_dirs = [
                np.get_include(),
                "/usr/include/openmpi-x86_64/"
        ]
)
