# distutils: language = c++
# Not sure if this is even worth it... But I don't like seeing yellow!
# Undo cython: cdivision=True, boundscheck=False
import numpy as np
cimport numpy as cnp
cimport cython

# Very useful: https://github.com/mpi4py/mpi4py/tree/master/demo/cython
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
        GenericIO(void *c, string filename) # Don't support the optional FIOT arg
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

        # Unused
        # long readTotalNumElems()
        # void getSourceRanks(vector[int] &SR);

        # void readDims(int dims[3])

        # void addVariable(string varname, long *data, unsigned int flags)




cdef class GenericIO_:
    cdef GenericIO *_thisptr

    # TODO: what is FIOT?
    def __cinit__(self, MPI.Comm world, str filename,
            bint should_compress = False, int partition = 0):
        self._thisptr = new GenericIO(
                world.ob_mpi,
                bytes(filename, "ascii"),
        )

        self._thisptr.setDefaultShouldCompress(should_compress)
        self._thisptr.setPartition(partition)

    def __dealloc__(self):
        del self._thisptr

    def write(self, cnp.ndarray toWrite):
        cdef int i
        cdef str colname
        cdef bytes colname_byt
        cdef type typ

        for i in range(len(toWrite.dtype)):
            colname = toWrite.dtype.names[i]
            colname_byt = bytes(colname, "ascii")
            contig = np.ascontiguousarray(toWrite[colname])
            typ = toWrite.dtype[i].type

            if typ is np.int32:
                self._addVariable[cnp.int32_t](contig, colname_byt)
            elif typ is np.int64:
                self._addVariable[cnp.int64_t](contig, colname_byt)
            elif typ is np.uint32:
                self._addVariable[cnp.uint32_t](contig, colname_byt)
            elif typ is np.uint64:
                self._addVariable[cnp.uint64_t](contig, colname_byt)
            elif typ is np.float32:
                self._addVariable[cnp.float32_t](contig, colname_byt)
            elif typ is np.float64:
                self._addVariable[cnp.float64_t](contig, colname_byt)
            else:
                raise Exception("Type not allowed")
            self._thisptr.setNumElems(len(toWrite))

        self._thisptr.write()

    def readHeader(self):
        self._thisptr.openAndReadHeader(GenericIO.MismatchBehavior.MismatchAllowed, -1, True)
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

    def readColumns(self, list colnames):
        cdef int world_rank, world_size
        mpi.MPI_Comm_rank(mpi.MPI_COMM_WORLD, &world_rank)
        mpi.MPI_Comm_size(mpi.MPI_COMM_WORLD, &world_size)

        cdef list colnames_byt = [bytes(cn, "ascii") for cn in colnames]

        # Find the cols we are looking. Error if they don't exist
        header_cols = self.readHeader()
        col_index = np.where(np.isin(header_cols["name"], colnames))[0]
        if len(col_index) != len(colnames):
            raise Exception("One or more cols not found: got {}, found {}".format(
                colnames, header_cols["name"][col_index]))

        # Get info about the rows
        cdef long num_writer_ranks = self._thisptr.readNRanks()
        assert (world_size >= num_writer_ranks,
            "Can't read with more ranks than were used to write")
        # Which of those ranks this reader should read
        # This could probably be ordered better if we know the topology is cartesian
        cdef long my_start_rank = world_rank * (num_writer_ranks / world_size)
        cdef long my_end_rank = (world_rank+1) * (num_writer_ranks / world_size)
        cdef long my_num_ranks = my_end_rank - my_start_rank

        cdef long [:] elems_in_rank = np.zeros(my_num_ranks, np.int64)
        cdef long rank
        for rank in range(my_num_ranks):
            elems_in_rank[rank] = self._thisptr.readNumElems(my_start_rank + rank)

        cdef long tot_rows = sum(elems_in_rank)
        cdef long extra_space = 5 # TODO why???
        cdef long max_rows = max(elems_in_rank) + extra_space

        cdef long idx
        results = np.zeros(tot_rows, dtype=
                [(header_cols[idx]["name"], header_cols[idx]["type"]) for idx in col_index])

        cdef int field_count
        for idx in col_index:
            field_count = header_cols[idx]["size"] / header_cols[idx]["elemsize"]
            colname_str = str(header_cols[idx]["name"]) # Previously was numpy.str_

            if header_cols[idx]["type"] == "f4":
                self._loadData[cnp.float32_t](
                        np.zeros(max_rows, "f4"), results[colname_str],
                        colname_str, field_count,
                        my_start_rank, my_num_ranks, elems_in_rank)
            elif header_cols[idx]["type"] == "f8":
                self._loadData[cnp.float64_t](
                        np.zeros(max_rows, "f8"), results[colname_str],
                        colname_str, field_count,
                        my_start_rank, my_num_ranks, elems_in_rank)
            elif header_cols[idx]["type"] == "i4":
                self._loadData[cnp.int32_t](
                        np.zeros(max_rows, "i4"), results[colname_str],
                        colname_str, field_count,
                        my_start_rank, my_num_ranks, elems_in_rank)
            elif header_cols[idx]["type"] == "i8":
                self._loadData[cnp.int64_t](
                        np.zeros(max_rows, "i8"), results[colname_str],
                        colname_str, field_count,
                        my_start_rank, my_num_ranks, elems_in_rank)
            elif header_cols[idx]["type"] == "u4":
                self._loadData[cnp.uint32_t](
                        np.zeros(max_rows, "u4"), results[colname_str],
                        colname_str, field_count,
                        my_start_rank, my_num_ranks, elems_in_rank)
            elif header_cols[idx]["type"] == "u8":
                self._loadData[cnp.uint64_t](
                        np.zeros(max_rows, "u8"), results[colname_str],
                        colname_str, field_count,
                        my_start_rank, my_num_ranks, elems_in_rank)
            else:
                raise Exception("Unknown type")
        return results

    def readColumn(self, str colname):
        return self.readColumns([colname])[colname]

    # Private
    cdef _addVariable(self, gio_numeric [:] data, bytes colname):
        self._thisptr.addScalarizedVariable(colname, &data[0], 1,
                (GenericIO.VariableFlags.VarHasExtraSpace & # No clue what this does
                GenericIO.VariableFlags.VarIsPhysCoordX)) # or this...

    cdef _loadData(self, gio_numeric [:] rank_data, gio_numeric [:] results,
            str colname, int field_count,
            long my_start_rank, long my_num_ranks, long [:] elems_in_rank):

        self._thisptr.clearVariables()
        self._thisptr.addScalarizedVariable(
                bytes(colname, "ascii"), &rank_data[0], field_count,
                GenericIO.VariableFlags.VarHasExtraSpace)

        cdef long loc = 0
        cdef long rank
        for rank in range(my_num_ranks):
            self._thisptr.readData(my_start_rank + rank, False, False)
            results[loc:loc + elems_in_rank[rank]] = rank_data[:elems_in_rank[rank]]
            loc += elems_in_rank[rank]

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
