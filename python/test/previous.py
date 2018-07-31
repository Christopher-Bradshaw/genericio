# Test reading a previously written file with a different number of ranks
import sys
import os
sys.path.insert(1, os.path.join(sys.path[0], '..'))

from mpi4py import MPI
import mpi_wrapper as wrapper
import numpy as np

# The previous one was written with 8
assert MPI.COMM_WORLD.Get_size() == 2
comm = MPI.COMM_WORLD
# comm = MPI.COMM_WORLD.Create_cart(
#         [2,2,1],
#         periods=[True, True, True],
# )

f = os.path.dirname(os.path.abspath(__file__)) + "/_data/split_file"
f = "/home/christopher/research/argonne/um_anl_catalogs/_data/core_cat/07_13_17.AlphaQ.499.coreproperties"

gio = wrapper.GenericIO_(comm, f)
        # should_compress=True, partition=comm.Get_rank())

out_data = gio.readColumns(["x", "y"])
print(len(out_data))
print(out_data)
