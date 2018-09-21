# Test basic reading and writing. Ensure that we get out what we put in
import sys
import os
sys.path.insert(1, os.path.join(sys.path[0], '..'))
import generic_io

import numpy as np
import pandas as pd


f = os.path.dirname(os.path.abspath(__file__)) + "/_data/no_mpi"

gio = generic_io.Generic_IO(f, None)
in_data = pd.DataFrame({
    "x": 1,
    "y": np.arange(10),
})
gio.write(in_data)

out_headers = gio.read_column_headers()
assert np.all(out_headers["name"] == np.array(["x", "y"]))

out_data = gio.read_columns(["x", "y"])
assert out_data.equals(in_data)

out_data = gio.read_columns(["x"])
assert out_data["x"].equals(in_data["x"])

out_data = gio.read_column("x")
assert out_data.equals(in_data["x"])
