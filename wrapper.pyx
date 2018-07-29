cdef extern from "GenericIO.h" namespace "gio":
    int x
    cdef cppclass GenericIO:
        GenericIO(char * filename, unsigned int FIOT)
        enum MismatchBehavior:
            MismatchAllowed
            MismatchDisallowed
            MismatchRedistribute

        void openAndReadHeader(
                MismatchBehavior MB,
                int EffRank, bint CheckPartMap)

        int getNumberOfVariables()
        long readTotalNumElems()
        long readNRanks()
        long readNumElems(int rank)


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

def test():
    return x
