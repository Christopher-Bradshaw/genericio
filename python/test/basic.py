# Test basic reading and writing. Ensure that we get out what we put in
import sys
import os
sys.path.insert(1, os.path.join(sys.path[0], '..'))

from mpi4py import MPI
import generic_io
import numpy as np
import pandas as pd

assert MPI.COMM_WORLD.Get_size() == 8

comm = MPI.COMM_WORLD.Create_cart(
        [2,2,2],
        periods=[True, True, True],
)

f = os.path.dirname(os.path.abspath(__file__)) + "/_data/basic"

gio = generic_io.Generic_IO(f, comm)
in_data = pd.DataFrame({
    "x": 1,
    "y": np.arange(comm.Get_rank() + 1),
})
gio.write(in_data)

out_headers = gio.read_header()
assert np.all(out_headers["name"] == np.array(["x", "y"]))
MPI.COMM_WORLD.barrier()

out_data = gio.read_columns(["x", "y"])
assert out_data.equals(in_data)

out_data = gio.read_columns(["x"])
assert out_data["x"].equals(in_data["x"])

out_data = gio.read_column("x")
assert out_data.equals(in_data["x"])
