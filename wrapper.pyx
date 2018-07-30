# distutils: language = c++
import numpy as np
cimport numpy as cnp
cimport cython

from libcpp.string cimport string
from libcpp.vector cimport vector

ctypedef fused gio_numeric:
    cython.int
    cython.long
    cython.float
    cython.double

cdef extern from "GenericIO.h" namespace "gio":
    int x
    cdef cppclass GenericIO:
        GenericIO(string filename, unsigned int FIOT)

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


        void openAndReadHeader(
                MismatchBehavior MB,
                int EffRank, bint CheckPartMap)

        int getNumberOfVariables()
        long readTotalNumElems()
        long readNRanks()
        long readNumElems(int rank)

        void getVariableInfo(vector[VariableInfo] &VI)
        void readDims(int dims[3])
        void readData(int EffRank, bint PrintStats, bint)

        void addVariable(string varname, long *data, unsigned int flags)
        void addScalarizedVariable[T](
                string varname, T *data, size_t numelems, unsigned int flags)

cdef class GenericIO_:
    cdef GenericIO *_thisptr

    def __cinit__(self, bytes filename, unsigned int FIOT):
        self._thisptr = new GenericIO(filename, FIOT)

    def readHeader(self):
        self._thisptr.openAndReadHeader(GenericIO.MismatchBehavior.MismatchAllowed, -1, True)
        # Get info about the cols
        cdef vector[GenericIO.VariableInfo] VI
        self._thisptr.getVariableInfo(VI)

        cols = np.zeros(VI.size(), dtype=[
            ("name", bytes, 100), ("size", np.int64), ("elemsize", np.int64)
        ])
        for i in range(VI.size()):
            cols[i] = (VI[i].Name, VI[i].Size, VI[i].ElementSize)
        return cols

    def readColumn(self, bytes colname):
        self._thisptr.openAndReadHeader(GenericIO.MismatchBehavior.MismatchAllowed, -1, True)

        # Get info about the rows
        cdef long num_ranks = self.readNRanks()
        cdef long [:] elems_in_rank = np.zeros(num_ranks, np.int64)
        for rank in range(num_ranks):
            elems_in_rank[rank] = self._thisptr.readNumElems(rank)

        cdef long tot_rows = sum(elems_in_rank)
        cdef long extra_space = 5 # TODO why???
        cdef long max_rows = max(elems_in_rank) + extra_space

        # Info about the cols
        header_cols = self.readHeader()
        col_index = np.where(header_cols["name"] == colname)[0]
        try:
            col_index = col_index[0]
        except IndexError:
            raise Exception("Colname not found in data")

        # TODO: Work out what this is for...
        # It is used in the num elems call. I think this is the number of elements per row.
        # e.g. a position three tuple would have size = 3 elemsize
        field_count = header_cols[col_index]["size"] / header_cols[col_index]["elemsize"]

        cdef long [:] rank_data = np.zeros(max_rows, np.int64)
        self._thisptr.addScalarizedVariable(colname, &rank_data[0], field_count,
                GenericIO.VariableFlags.VarHasExtraSpace)

        cdef long [:] results = np.zeros(tot_rows, np.int64)
        cdef long loc = 0
        for rank in range(num_ranks):
            self._thisptr.readData(rank, False, False)
            results[loc:loc + elems_in_rank[rank]] = rank_data[:elems_in_rank[rank]]
            loc += elems_in_rank[rank]

        return np.array(results)


    ############ Probably shouldn't be part of the public interface

    def openAndReadHeader(self, GenericIO.MismatchBehavior MB, int EffRank, bint CheckPartMap):
        # pass
        return self._thisptr.openAndReadHeader(MB, EffRank, CheckPartMap)

    def getNumberOfVariables(self):
        return self._thisptr.getNumberOfVariables()

    def readTotalNumElems(self):
        return self._thisptr.readTotalNumElems()

    def readNumElems(self, rank):
        return self._thisptr.readNumElems(rank)

    def readNRanks(self):
        return self._thisptr.readNRanks()


    def readDims(self):
        cdef int [:] dims = np.zeros(3, np.int32)
        self._thisptr.readDims(&dims[0])
        return np.array(dims)


    # The process of reading is something like:
    # 1) Add the variable you want to read
    # 2) read it
    # 3) get var info and in the .data field there is what you read
    def addVariable(self, name):
        cdef long [:] data = np.zeros(10, np.int64)
        self._thisptr.addVariable(name, &data[0], 1)

    def readData(self):
        self._thisptr.readData(0, True, True)

    def getVariableInfo(self):
        cdef vector[GenericIO.VariableInfo] VI
        self._thisptr.getVariableInfo(VI)
        print(VI[0].Name)
        print(VI.size())
