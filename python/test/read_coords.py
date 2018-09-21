# Test basic reading and writing. Ensure that we get out what we put in
import sys
import os
sys.path.insert(1, os.path.join(sys.path[0], '..'))

from mpi4py import MPI
import generic_io
import numpy as np
import pandas as pd


f = os.path.dirname(os.path.abspath(__file__)) + "/_data/input/07_13_17.AlphaQ.42.coreproperties"

gio = generic_io.Generic_IO(f, None)

metadata = gio.read_metadata(1)
expected = {'num_ranks': 256, 'coords': [0, 0, 1], 'dims': [8, 8, 4], 'origin': [0.0, 0.0, 0.0], 'scale': [256.0, 256.0, 256.0], 'global_rank': 1}
assert metadata.keys() == expected.keys()
for k in metadata.keys():
    assert metadata[k] == expected[k], "Error for key: {}".format(k)
