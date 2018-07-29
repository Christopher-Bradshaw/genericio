# distutils: language = c++
import numpy as np
cimport numpy as cnp

from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "GenericIO.h" namespace "gio":
    int x
    cdef cppclass GenericIO:
        GenericIO(string filename, unsigned int FIOT)

        enum MismatchBehavior:
            MismatchAllowed
            MismatchDisallowed
            MismatchRedistribute

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


cdef class GenericIO_:
    cdef GenericIO *_thisptr

    def __cinit__(self, bytes filename, unsigned int FIOT):
        self._thisptr = new GenericIO(filename, FIOT)

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



def test():
    return x
