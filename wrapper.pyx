# distutils: language = c++
import numpy as np
cimport numpy as cnp
cimport cython
from cython.view cimport array as cvarray


from libcpp.string cimport string
from libcpp.vector cimport vector

ctypedef fused gio_numeric:
    cnp.int64_t
    cnp.int32_t
    cnp.uint64_t
    cnp.uint32_t
    cnp.float64_t
    cnp.float32_t

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
        void clearVariables()

cdef class GenericIO_:
    cdef GenericIO *_thisptr

    def __cinit__(self, bytes filename, unsigned int FIOT):
        self._thisptr = new GenericIO(filename, FIOT)

    def readHeader(self):
        self._thisptr.openAndReadHeader(GenericIO.MismatchBehavior.MismatchAllowed, -1, True)
        # Get info about the cols
        cdef vector[GenericIO.VariableInfo] vi
        self._thisptr.getVariableInfo(vi)

        cols = np.zeros(vi.size(), dtype=[
            ("name", bytes, 100),
            ("size", np.int64),
            ("elemsize", np.int64),
            ("type", bytes, 2),
        ])
        for i in range(vi.size()):
            cols[i] = (vi[i].Name, vi[i].Size, vi[i].ElementSize, self._type_from_variable_info(vi[i]))
        return cols

    def readColumns(self, colnames):
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
        col_index = np.where(np.isin(header_cols["name"], colnames))[0]
        if len(col_index) != len(colnames):
            raise Exception("One or more cols not found: got {}, found {}".format(
                colnames, header_cols["name"][col_index]))

        results = np.zeros(tot_rows, dtype=[
            (
                header_cols[idx]["name"].decode("utf-8"),
                header_cols[idx]["type"].decode("utf-8"),
            ) for idx in col_index])

        cdef int field_count
        for idx in col_index:
            field_count = header_cols[idx]["size"] / header_cols[idx]["elemsize"]
            colname_str = header_cols[idx]["name"].decode("utf-8")
            colname_byt = bytes(header_cols[idx]["name"])
            print(colname_str)

            if header_cols[idx]["type"] == b"f4":
                results[colname_str] = self._loadData(<cnp.float32_t>1, "f4",
                    max_rows, tot_rows, colname_byt, field_count, num_ranks, elems_in_rank)
            elif header_cols[idx]["type"] == b"f8":
                results[colname_str] = self._loadData(<cnp.float64_t>1, "f8",
                    max_rows, tot_rows, colname_byt, field_count, num_ranks, elems_in_rank)
            elif header_cols[idx]["type"] == b"i4":
                results[colname_str] = self._loadData(<cnp.int32_t>1, "i4",
                    max_rows, tot_rows, colname_byt, field_count, num_ranks, elems_in_rank)
            elif header_cols[idx]["type"] == b"i8":
                print("here-1")
                results[colname_str] = self._loadData(<cnp.int64_t>1, "i8",
                    max_rows, tot_rows, colname_byt, field_count, num_ranks, elems_in_rank)
            elif header_cols[idx]["type"] == b"u4":
                results[colname_str] = self._loadData(<cnp.uint32_t>1, "u4",
                    max_rows, tot_rows, colname_byt, field_count, num_ranks, elems_in_rank)
            elif header_cols[idx]["type"] == b"u8":
                results[colname_str] = self._loadData(<cnp.uint64_t>1, "u8",
                    max_rows, tot_rows, colname_byt, field_count, num_ranks, elems_in_rank)
            else:
                raise Exception("Unknown type")
        return results

    def readColumn(self, bytes colname):
        return self.readColumns([colname])[colname.decode("utf-8")]

    cdef _loadData(self, gio_numeric a, str python_type,
            long max_rows, long tot_rows, bytes colname,
            int field_count, long num_ranks, elems_in_rank):

        cdef gio_numeric [:] rank_data = np.zeros(max_rows, dtype=python_type)
        cdef gio_numeric [:] results = np.zeros(tot_rows, dtype=python_type)
        # It seems like we can read more than one col at a time!
        print(self._thisptr.clearVariables())
        self._thisptr.addScalarizedVariable(colname, &rank_data[0], field_count,
                GenericIO.VariableFlags.VarHasExtraSpace)

        cdef long loc = 0
        for rank in range(num_ranks):
            self._thisptr.readData(rank, False, False)
            results[loc:loc + elems_in_rank[rank]] = rank_data[:elems_in_rank[rank]]
            loc += elems_in_rank[rank]

        return results

    cdef bytes _type_from_variable_info(self, GenericIO.VariableInfo vi):
        if vi.IsFloat and vi.Size == 4:
            return b"f4"
        elif vi.IsFloat and vi.Size == 8:
            return b"f8"
        elif vi.IsSigned and vi.Size == 4:
            return b"i4"
        elif vi.IsSigned and vi.Size == 8:
            return b"i8"
        elif not vi.IsSigned and vi.Size == 4:
            return b"u4"
        elif not vi.IsSigned and vi.Size == 8:
            return b"u8"


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
