# Test basic reading and writing. Ensure that we get out what we put in
# Do this with the split file format.
import sys
import os
sys.path.insert(1, os.path.join(sys.path[0], '..'))

from mpi4py import MPI
import mpi_wrapper as wrapper
import numpy as np
import pandas as pd


assert MPI.COMM_WORLD.Get_size() == 8

comm = MPI.COMM_WORLD.Create_cart(
        [2,2,2],
        periods=[True, True, True],
)

f = os.path.dirname(os.path.abspath(__file__)) + "/_data/split_file"

gio = wrapper.GenericIO_(comm, f,
        should_compress=True,
        partition=comm.Get_rank())


data_len = 1000000

in_data = pd.DataFrame({
        "x": np.arange(data_len) * (comm.Get_rank() + 1),
})
in_data["y"] = np.sqrt(in_data["x"])


gio.write(in_data)

out_headers = gio.readHeader()
assert np.all(out_headers["name"] == np.array(["x", "y"]))

# Reading before we were done writing would be bad...
MPI.COMM_WORLD.barrier()

out_data = gio.readColumns(["x", "y"])
assert out_data.equals(in_data)

out_data = gio.readColumns(["x"])
assert out_data["x"].equals(in_data["x"])

out_data = gio.readColumn("x")
assert out_data.equals(in_data["x"])
