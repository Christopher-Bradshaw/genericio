# GenericIO

Forked from https://trac.alcf.anl.gov/projects/genericio

This fork includes a python (N.B. only tested with python3) wrapper to the C++ api. You should probably check that the C++ api hasn't changed. If this is useful it might be useful to merge back.

## Getting started

To install, make sure that you:
1. Have pip installed (we need install some dependencies with pip, not conda)
2. Are in a virtualenv (so pip doesn't try and fail to install in the root python dir)
3. Have an mpi compiler and mpi headers available (e.g. via `module load ...`)

Then run:

```
make py_deps
make py_build
make py_test
```

This installs the python dependencies (cython, numpy, pandas and mpi4py), builds the object files that the python wrapper depends on (Blosc and GenericIO), compiles the cython wrapper into C and then into a shared object and then links that all into a python module that can be imported.

If this all worked fine, you should see the GenericIO logs from the tests.


## API

The tests are probably the easiest way to learn the API. To see all features that are supported you probably just need to read the [wrapper code](./python/wrapper.pyx). A tiny example is:

```
import generic_io

gio = generic_io.Generic_IO("some/file, MPI_Comm)
data = numpy_structured_array or pandas_dataframe

gio.write(data) # Write this data to the file

col_headers = gio.read_column_headers() # Get info about the column names and types

read_data = gio.read_columns() # Read all the data
read_data = gio.read_columns(["x", "y"]) # Read only column x and y
read_data = gio.read_column("x") # Read only column x

metadata = gio.read_metadata() # Dict with General info about the gio file
# E.g. (number of ranks, location of origin, scale) etc
```
