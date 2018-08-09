# Test reading a previously written file with a different number of ranks
import sys
import os
sys.path.insert(1, os.path.join(sys.path[0], '..'))

from mpi4py import MPI
import generic_io
import numpy as np

# The previous one was written with 8
assert MPI.COMM_WORLD.Get_size() == 4
comm = MPI.COMM_WORLD.Create_cart(
        [2,2,1],
        periods=[True, True, True],
)

f = os.path.dirname(os.path.abspath(__file__)) + "/_data/split_file"

gio = generic_io.Generic_IO(f, comm)

out_data = gio.read_columns(["x", "y"], ranks=[comm.Get_rank(), comm.Get_rank() + 4])
data_len = 1000000
assert np.all(
        out_data["x"] == np.concatenate((
            np.arange(data_len) * (comm.Get_rank() + 1),
            np.arange(data_len) * (comm.Get_rank() + 5)))
)
