# Test basic reading and writing. Ensure that we get out what we put in
import sys
import os
sys.path.insert(1, os.path.join(sys.path[0], '..'))

from mpi4py import MPI
import mpi_wrapper as wrapper
import numpy as np

assert MPI.COMM_WORLD.Get_size() == 8

comm = MPI.COMM_WORLD.Create_cart(
        [2,2,2],
        periods=[True, True, True],
)

f = os.path.dirname(os.path.abspath(__file__)) + "/_data/new"

gio = wrapper.GenericIO_(comm, f, 1)
in_data = np.zeros(comm.Get_rank() + 1, dtype=[("x", "f8"), ("y", "i8")])
in_data["x"] = 1
in_data["y"] = np.arange(comm.Get_rank() + 1)
gio.write(in_data)

out_headers = gio.readHeader()
assert np.all(out_headers["name"] == np.array(["x", "y"]))
MPI.COMM_WORLD.barrier()

out_data = gio.readColumns(["x", "y"])
assert np.all(out_data == in_data)

out_data = gio.readColumns(["x"])
assert np.all(out_data["x"] == in_data["x"])

out_data = gio.readColumn("x")
assert np.all(out_data == in_data["x"])
