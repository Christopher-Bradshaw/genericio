# distutils: language = c++
# cython: profile=True
# Not sure if this is even worth it... But I don't like seeing yellow!
# Undo cython: cdivision=True, boundscheck=False

# https://stackoverflow.com/questions/40845304
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

import numpy as np
import pandas as pd

cimport numpy as cnp
cimport cython

# Very useful: https://github.com/mpi4py/mpi4py/tree/master/demo/cython
from mpi4py import MPI
from mpi4py cimport MPI
from mpi4py.MPI cimport Intracomm as IntracommType
from mpi4py cimport libmpi as mpi

from libcpp.string cimport string
from libcpp.vector cimport vector

ctypedef fused gio_numeric:
    cnp.int64_t
    cnp.int32_t
    cnp.uint64_t
    cnp.uint32_t
    cnp.float64_t
    cnp.float32_t

cdef extern from "../GenericIO.h" namespace "gio":
    cdef cppclass GenericIO:

        # TODO: Work out why we need these string tags to be able to access these types
        # I currently have *no idea*...
        # See https://github.com/cython/cython/issues/1603
        enum MismatchBehavior:
            MismatchAllowed "gio::GenericIO::MismatchBehavior::MismatchAllowed",
            MismatchDisallowed "gio::GenericIO::MismatchBehavior::MismatchDisallowed",
            MismatchRedistribute "gio::GenericIO::MismatchBehavior::MismatchRedistribute",

        enum VariableFlags:
            VarHasExtraSpace "gio::GenericIO::VariableFlags::VarHasExtraSpace" = (1 << 0),
            VarIsPhysCoordX "gio::GenericIO::VariableFlags::VarIsPhysCoordX" = (1 << 1),
            VarIsPhysCoordY "gio::GenericIO::VariableFlags::VarIsPhysCoordY" = (1 << 2),
            VarIsPhysCoordZ "gio::GenericIO::VariableFlags::VarIsPhysCoordZ" = (1 << 3),
            VarMaybePhysGhost "gio::GenericIO::VariableFlags::VarMaybePhysGhost" = (1 << 4),

        struct VariableInfo:
            string Name
            size_t Size
            bint IsFloat
            bint IsSigned
            bint IsPhysCoordX, IsPhysCoordY, IsPhysCoordZ
            bint MaybePhysGhost
            size_t ElementSize

        # Initialization
        GenericIO(void *c, string filename) except +# Don't support the optional FIOT arg
        void setDefaultShouldCompress(bint shouldCompress)
        void setPartition(int partition)

        # Variables
        void setNumElems(size_t num_elems)
        void addScalarizedVariable[T](string varname, T *data, size_t num_elems,
                unsigned int flags)
        void clearVariables()
        void getVariableInfo(vector[VariableInfo] &VI)

        # Writing
        void write()

        # Reading
        void openAndReadHeader(MismatchBehavior MB, int EffRank, bint CheckPartMap)
        long readNRanks() # The number of ranks that wrote the file
        long readNumElems(int rank)
        void readData(int EffRank, bint PrintStats, bint)

        # Arb
        size_t requestedExtraSpace()

        # Unclear/testing
        void readCoords(int coords[3], int EffRank)
        void readDims(int dims[3])
        void readPhysOrigin(double orig[3])
        void readPhysScale(double scale[3])
        int readGlobalRankNumber(int EffRank)

        # Unused
        # long readTotalNumElems()
        # void getSourceRanks(vector[int] &SR);

        # void readDims(int dims[3])

        # void addVariable(string varname, long *data, unsigned int flags)




cdef class Generic_IO:
    cdef GenericIO *_thisptr
    cdef bint verbose
    cdef bint header_is_read

    def __cinit__(self, str filename, MPI.Comm world,
            bint should_compress = False, int partition = 0, bint verbose = False):
        if world is None:
            world = MPI.COMM_WORLD
        self._thisptr = new GenericIO(
                world.ob_mpi,
                bytes(filename, "ascii"),
        )

        self._thisptr.setDefaultShouldCompress(should_compress)
        self._thisptr.setPartition(partition)

        self.verbose = verbose
        self.header_is_read = False

    def __dealloc__(self):
        del self._thisptr

    def read_meta(self, rank = None):
        if not self.header_is_read:
            self._thisptr.openAndReadHeader(GenericIO.MismatchBehavior.MismatchAllowed, -1, True)
            self.header_is_read = True

        cdef int world_rank
        if rank is None:
            mpi.MPI_Comm_rank(mpi.MPI_COMM_WORLD, &world_rank)
            ranks = world_rank

        # Number of ranks
        NR = self._thisptr.readNRanks()

        # Coords
        cdef int[3] coords = np.zeros(3, np.int32)
        self._thisptr.readCoords(coords, rank)

        # Dims
        cdef int[3] dims = np.zeros(3, np.int32)
        self._thisptr.readDims(dims)

        # Orig
        cdef double[3] orig = np.zeros(3, np.float64)
        self._thisptr.readPhysOrigin(orig)

        # Scale
        cdef double[3] scale = np.zeros(3, np.float64)
        self._thisptr.readPhysScale(scale)

        # Global rank number
        gr = self._thisptr.readGlobalRankNumber(rank)

        # These aren't that well named
        return {
                "num_ranks": NR,
                "coords": coords,
                "dims": dims,
                "origin": orig,
                "scale": scale,
                "global_rank": gr,
        }

    def write(self, to_write):
        if type(to_write) is np.ndarray:
            to_write = pd.DataFrame(to_write)

        assert type(to_write) is pd.core.frame.DataFrame
        cdef int i
        cdef str colname
        cdef type typ

        cdef list tmp = []

        for i in range(len(to_write.columns)):
            colname = to_write.columns[i]
            data = to_write[colname].values
            typ = data.dtype.type
            assert data.flags["C_CONTIGUOUS"]

            if typ is np.int32:
                self._add_variable[cnp.int32_t](data, colname)
            elif typ is np.int64:
                self._add_variable[cnp.int64_t](data, colname)
            elif typ is np.uint32:
                self._add_variable[cnp.uint32_t](data, colname)
            elif typ is np.uint64:
                self._add_variable[cnp.uint64_t](data, colname)
            elif typ is np.float32:
                self._add_variable[cnp.float32_t](data, colname)
            elif typ is np.float64:
                self._add_variable[cnp.float64_t](data, colname)
            else:
                raise Exception("Type not allowed")

        self._thisptr.setNumElems(len(to_write))
        self._thisptr.write()

    # Should maybe rename to something like read_variable_info
    def read_header(self):
        if not self.header_is_read:
            self._thisptr.openAndReadHeader(GenericIO.MismatchBehavior.MismatchAllowed, -1, True)
            self.header_is_read = True

        cdef vector[GenericIO.VariableInfo] vi
        self._thisptr.getVariableInfo(vi)

        cols = np.zeros(vi.size(), dtype=[
            ("name", str, 100),
            ("size", np.int64),
            ("elemsize", np.int64),
            ("type", str, 2),
        ])
        for i in range(vi.size()):
            cols[i] = (
                    vi[i].Name,
                    vi[i].Size,
                    vi[i].ElementSize,
                    self._type_from_variable_info(vi[i]))
        return cols

    def read_columns(self, colnames = None, ranks = None, bint as_numpy_array = False):
        # This always needs to go first
        header_cols = self.read_header()

        cdef int world_rank, world_size
        mpi.MPI_Comm_size(mpi.MPI_COMM_WORLD, &world_size)
        assert (world_size >= self._thisptr.readNRanks(),
            "Can't read with more ranks than were used to write")

        # If not specified, read your own rank
        mpi.MPI_Comm_rank(mpi.MPI_COMM_WORLD, &world_rank)
        if ranks is None:
            ranks = [world_rank]

        # Find the cols we are looking for. Error if they don't exist
        colnames = colnames or header_cols["name"]
        col_index = np.where(np.isin(header_cols["name"], colnames))[0]
        if len(col_index) != len(colnames):
            raise Exception("One or more cols not found: got {}, found {}".format(
                colnames, header_cols["name"][col_index]))

        # Get info about the rows
        cdef long [:] elems_in_rank = np.zeros(len(ranks), np.int64)
        cdef long rank
        for i in range(len(ranks)):
            elems_in_rank[i] = self._thisptr.readNumElems(ranks[i])

        cdef long tot_rows = sum(elems_in_rank)
        # TODO why??? Also this is very important. Was having intermittant crashes with this == 5
        # Why was it set to five? Wild guess
        cdef long extra_space = self._thisptr.requestedExtraSpace()
        cdef long max_rows = max(elems_in_rank) + extra_space

        cdef long idx
        results = pd.DataFrame()

        cdef int field_count
        for idx in col_index:
            field_count = header_cols[idx]["size"] / header_cols[idx]["elemsize"]
            colname_str = str(header_cols[idx]["name"]) # Previously was numpy.str_
            results[colname_str] = np.zeros(tot_rows, dtype=header_cols[idx]["type"])
            data = results[colname_str].values

            if header_cols[idx]["type"] == "f4":
                self._load_data[cnp.float32_t](
                        np.zeros(max_rows, "f4"), data,
                        colname_str, field_count,
                        ranks, elems_in_rank)
            elif header_cols[idx]["type"] == "f8":
                self._load_data[cnp.float64_t](
                        np.zeros(max_rows, "f8"), data,
                        colname_str, field_count,
                        ranks, elems_in_rank)
            elif header_cols[idx]["type"] == "i4":
                self._load_data[cnp.int32_t](
                        np.zeros(max_rows, "i4"), data,
                        colname_str, field_count,
                        ranks, elems_in_rank)
            elif header_cols[idx]["type"] == "i8":
                self._load_data[cnp.int64_t](
                        np.zeros(max_rows, "i8"), data,
                        colname_str, field_count,
                        ranks, elems_in_rank)
            elif header_cols[idx]["type"] == "u4":
                self._load_data[cnp.uint32_t](
                        np.zeros(max_rows, "u4"), data,
                        colname_str, field_count,
                        ranks, elems_in_rank)
            elif header_cols[idx]["type"] == "u8":
                self._load_data[cnp.uint64_t](
                        np.zeros(max_rows, "u8"), data,
                        colname_str, field_count,
                        ranks, elems_in_rank)
            else:
                raise Exception("Unknown type")

        if as_numpy_array:
            recarr = results.to_records(index=False)
            structarr = recarr.view(recarr.dtype.fields, np.ndarray)
            return structarr

        return results

    def read_column(self, str colname, ranks = None, bint as_numpy_array = False):
        return self.read_columns([colname], ranks, as_numpy_array)[colname]

    # return the number of writers a file had
    def num_writers(self):
        if not self.header_is_read:
            self._thisptr.openAndReadHeader(GenericIO.MismatchBehavior.MismatchAllowed, -1, True)
            self.header_is_read = True
        return self._thisptr.readNRanks()

    # Private
    cdef _add_variable(self, gio_numeric [:] data, str colname):
        self._thisptr.addScalarizedVariable(
                bytes(colname, "ascii"),
                &data[0] if len(data) else NULL, # handle empty input
                1,
                0,
        )
                # (GenericIO.VariableFlags.VarHasExtraSpace & # No clue what this does
                # GenericIO.VariableFlags.VarIsPhysCoordX)) # or this...

    cdef _load_data(self, gio_numeric [:] rank_data, gio_numeric [:] results,
            str colname, int field_count,
            list ranks, long [:] elems_in_rank):

        self._thisptr.clearVariables()
        self._thisptr.addScalarizedVariable(
                bytes(colname, "ascii"),
                &rank_data[0],
                field_count,
                GenericIO.VariableFlags.VarHasExtraSpace)
                # We don't need this here as we have extraSpace. But we might one day when
                # I understand that...
                # &rank_data[0] if len(rank_data) else NULL,

        cdef long loc = 0
        cdef long rank
        for i in range(len(ranks)):
            self._thisptr.readData(ranks[i], self.verbose, self.verbose)
            results[loc:loc + elems_in_rank[i]] = rank_data[:elems_in_rank[i]]
            loc += elems_in_rank[i]

        return results

    cdef str _type_from_variable_info(self, GenericIO.VariableInfo vi):
        if vi.IsFloat and vi.Size == 4:
            return "f4"
        elif vi.IsFloat and vi.Size == 8:
            return "f8"
        elif vi.IsSigned and vi.Size == 4:
            return "i4"
        elif vi.IsSigned and vi.Size == 8:
            return "i8"
        elif not vi.IsSigned and vi.Size == 4:
            return "u4"
        elif not vi.IsSigned and vi.Size == 8:
            return "u8"
