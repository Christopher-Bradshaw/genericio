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
in_data = np.zeros(1000000, dtype=[("x", "i8"), ("y", "f8")])
in_data["x"] = np.arange(len(in_data)) * (comm.Get_rank() + 1)
in_data["y"] = np.sqrt(in_data["x"])


toWrite = pd.DataFrame(in_data)

gio.write(toWrite)

# out_headers = gio.readHeader()
# assert np.all(out_headers["name"] == np.array(["x", "y"]))

# # Reading before we were done writing would be bad...
# MPI.COMM_WORLD.barrier()

# out_data = gio.readColumns(["x", "y"])
# assert np.all(out_data == in_data)

# out_data = gio.readColumns(["x"])
# assert np.all(out_data["x"] == in_data["x"])

# out_data = gio.readColumn("x")
# assert np.all(out_data == in_data["x"])
