#define _XOPEN_SOURCE 600
#include "CRC64.h"
#include "GenericIO.h"

extern "C" {
#include "blosc.h"
}

#include <sstream>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstring>

#ifndef GENERICIO_NO_MPI
#include <ctime>
#endif

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>

#ifdef __bgq__
#include <spi/include/kernel/location.h>
#include <spi/include/kernel/process.h>
#include <firmware/include/personality.h>
#endif

#ifndef MPI_UINT64_T
#define MPI_UINT64_T (sizeof(long) == 8 ? MPI_LONG : MPI_LONG_LONG)
#endif

using namespace std;

namespace gio {


#ifndef GENERICIO_NO_MPI
GenericFileIO_MPI::~GenericFileIO_MPI() {
  (void) MPI_File_close(&FH);
}

void GenericFileIO_MPI::open(const std::string &FN, bool ForReading) {
  FileName = FN;

  int amode = ForReading ? MPI_MODE_RDONLY : (MPI_MODE_WRONLY | MPI_MODE_CREATE);
  if (MPI_File_open(Comm, const_cast<char *>(FileName.c_str()), amode,
                    MPI_INFO_NULL, &FH) != MPI_SUCCESS)
    throw runtime_error((!ForReading ? "Unable to create the file: " :
                                       "Unable to open the file: ") +
                        FileName);
}

void GenericFileIO_MPI::setSize(size_t sz) {
  if (MPI_File_set_size(FH, sz) != MPI_SUCCESS)
    throw runtime_error("Unable to set size for file: " + FileName);
}

void GenericFileIO_MPI::read(void *buf, size_t count, off_t offset,
                             const std::string &D) {
  while (count > 0) {
    MPI_Status status;
    if (MPI_File_read_at(FH, offset, buf, count, MPI_BYTE, &status) != MPI_SUCCESS)
      throw runtime_error("Unable to read " + D + " from file: " + FileName);

    int scount;
    (void) MPI_Get_count(&status, MPI_BYTE, &scount);

    count -= scount;
    buf = ((char *) buf) + scount;
    offset += scount;
  }
}

void GenericFileIO_MPI::write(const void *buf, size_t count, off_t offset,
                              const std::string &D) {
  while (count > 0) {
    MPI_Status status;
    if (MPI_File_write_at(FH, offset, (void *) buf, count, MPI_BYTE, &status) != MPI_SUCCESS)
      throw runtime_error("Unable to write " + D + " to file: " + FileName);

    int scount;
    (void) MPI_Get_count(&status, MPI_BYTE, &scount);

    count -= scount;
    buf = ((char *) buf) + scount;
    offset += scount;
  }
}

void GenericFileIO_MPICollective::read(void *buf, size_t count, off_t offset,
                             const std::string &D) {
  int Continue = 0;

  do {
    MPI_Status status;
    if (MPI_File_read_at_all(FH, offset, buf, count, MPI_BYTE, &status) != MPI_SUCCESS)
      throw runtime_error("Unable to read " + D + " from file: " + FileName);

    int scount;
    (void) MPI_Get_count(&status, MPI_BYTE, &scount);

    count -= scount;
    buf = ((char *) buf) + scount;
    offset += scount;

    int NeedContinue = (count > 0);
    MPI_Allreduce(&NeedContinue, &Continue, 1, MPI_INT, MPI_SUM, Comm);
  } while (Continue);
}

void GenericFileIO_MPICollective::write(const void *buf, size_t count, off_t offset,
                              const std::string &D) {
  int Continue = 0;

  do {
    MPI_Status status;
    if (MPI_File_write_at_all(FH, offset, (void *) buf, count, MPI_BYTE, &status) != MPI_SUCCESS)
      throw runtime_error("Unable to write " + D + " to file: " + FileName);

    int scount;
    (void) MPI_Get_count(&status, MPI_BYTE, &scount);

    count -= scount;
    buf = ((char *) buf) + scount;
    offset += scount;

    int NeedContinue = (count > 0);
    MPI_Allreduce(&NeedContinue, &Continue, 1, MPI_INT, MPI_SUM, Comm);
  } while (Continue);
}
#endif

GenericFileIO_POSIX::~GenericFileIO_POSIX() {
  if (FH != -1) close(FH);
}

void GenericFileIO_POSIX::open(const std::string &FN, bool ForReading) {
  FileName = FN;

  int flags = ForReading ? O_RDONLY : (O_WRONLY | O_CREAT);
  int mode = S_IRUSR | S_IWUSR | S_IRGRP;
  errno = 0;
  if ((FH = ::open(FileName.c_str(), flags, mode)) == -1)
    throw runtime_error((!ForReading ? "Unable to create the file: " :
                                       "Unable to open the file: ") +
                        FileName + ": " + strerror(errno));
}

void GenericFileIO_POSIX::setSize(size_t sz) {
  if (ftruncate(FH, sz) == -1)
    throw runtime_error("Unable to set size for file: " + FileName);
}

void GenericFileIO_POSIX::read(void *buf, size_t count, off_t offset,
                               const std::string &D) {
  while (count > 0) {
    ssize_t scount;
    errno = 0;
    if ((scount = pread(FH, buf, count, offset)) == -1) {
      if (errno == EINTR)
        continue;

      throw runtime_error("Unable to read " + D + " from file: " + FileName);
    }

    count -= scount;
    buf = ((char *) buf) + scount;
    offset += scount;
  }
}

void GenericFileIO_POSIX::write(const void *buf, size_t count, off_t offset,
                                const std::string &D) {
  while (count > 0) {
    ssize_t scount;
    errno = 0;
    if ((scount = pwrite(FH, buf, count, offset)) == -1) {
      if (errno == EINTR)
        continue;

      throw runtime_error("Unable to write " + D + " to file: " + FileName);
    }

    count -= scount;
    buf = ((char *) buf) + scount;
    offset += scount;
  }
}

static bool isBigEndian() {
  const uint32_t one = 1;
  return !(*((char *)(&one)));
}

static const size_t CRCSize = 8;

static const size_t MagicSize = 8;
static const char *MagicBE = "HACC01B";
static const char *MagicLE = "HACC01L";

struct GlobalHeader {
  char Magic[MagicSize];
  uint64_t HeaderSize;
  uint64_t NElems; // The global total
  uint64_t Dims[3];
  uint64_t NVars;
  uint64_t VarsSize;
  uint64_t VarsStart;
  uint64_t NRanks;
  uint64_t RanksSize;
  uint64_t RanksStart;
  uint64_t GlobalHeaderSize;
  double   PhysOrigin[3];
  double   PhysScale[3];
  uint64_t BlocksSize;
  uint64_t BlocksStart;
} __attribute__((packed));

enum {
  FloatValue          = (1 << 0),
  SignedValue         = (1 << 1),
  ValueIsPhysCoordX   = (1 << 2),
  ValueIsPhysCoordY   = (1 << 3),
  ValueIsPhysCoordZ   = (1 << 4),
  ValueMaybePhysGhost = (1 << 5)
};

static const size_t NameSize = 256;
struct VariableHeader {
  char Name[NameSize];
  uint64_t Flags;
  uint64_t Size;
} __attribute__((packed));

struct RankHeader {
  uint64_t Coords[3];
  uint64_t NElems;
  uint64_t Start;
  uint64_t GlobalRank;
} __attribute__((packed));

static const size_t FilterNameSize = 8;
static const size_t MaxFilters = 4;
struct BlockHeader {
  char Filters[MaxFilters][FilterNameSize];
  uint64_t Start;
  uint64_t Size;
} __attribute__((packed));

struct CompressHeader {
  uint64_t OrigCRC;
} __attribute__((packed));
const char *CompressName = "BLOSC";

unsigned GenericIO::DefaultFileIOType = FileIOPOSIX;
int GenericIO::DefaultPartition = 0;
bool GenericIO::DefaultShouldCompress = false;

#ifndef GENERICIO_NO_MPI
std::size_t GenericIO::CollectiveMPIIOThreshold = 0;
#endif

static bool blosc_initialized = false;

#ifndef GENERICIO_NO_MPI
// Note: writing errors are not currently recoverable (one rank may fail
// while the others don't).
void GenericIO::write() {
  const char *Magic = isBigEndian() ? MagicBE : MagicLE;

  uint64_t FileSize = 0;

  int NRanks, Rank;
  MPI_Comm_rank(Comm, &Rank);
  MPI_Comm_size(Comm, &NRanks);

#ifdef __bgq__
  MPI_Barrier(Comm);
#endif
  MPI_Comm_split(Comm, Partition, Rank, &SplitComm);

  int SplitNRanks, SplitRank;
  MPI_Comm_rank(SplitComm, &SplitRank);
  MPI_Comm_size(SplitComm, &SplitNRanks);

  string LocalFileName;
  if (SplitNRanks != NRanks) {
    if (Rank == 0) {
      // In split mode, the specified file becomes the rank map, and the real
      // data is partitioned.

      vector<int> MapRank, MapPartition;
      MapRank.resize(NRanks);
      for (int i = 0; i < NRanks; ++i) MapRank[i] = i;

      MapPartition.resize(NRanks);
      MPI_Gather(&Partition, 1, MPI_INT, &MapPartition[0], 1, MPI_INT, 0, Comm);

      GenericIO GIO(MPI_COMM_SELF, FileName, FileIOType);
      GIO.setNumElems(NRanks);
      GIO.addVariable("$rank", MapRank); /* this is for use by humans; the reading
                                            code assumes that the partitions are in
                                            rank order */
      GIO.addVariable("$partition", MapPartition);

      vector<int> CX, CY, CZ;
      int TopoStatus;
      MPI_Topo_test(Comm, &TopoStatus);
      if (TopoStatus == MPI_CART) {
        CX.resize(NRanks);
        CY.resize(NRanks);
        CZ.resize(NRanks);

        for (int i = 0; i < NRanks; ++i) {
          int C[3];
          MPI_Cart_coords(Comm, i, 3, C);

          CX[i] = C[0];
          CY[i] = C[1];
          CZ[i] = C[2];
        }

        GIO.addVariable("$x", CX);
        GIO.addVariable("$y", CY);
        GIO.addVariable("$z", CZ);
      }

      GIO.write();
    } else {
      MPI_Gather(&Partition, 1, MPI_INT, 0, 0, MPI_INT, 0, Comm);
    }

    stringstream ss;
    ss << FileName << "#" << Partition;
    LocalFileName = ss.str();
  } else {
    LocalFileName = FileName;
  }

  RankHeader RHLocal;
  int Dims[3], Periods[3], Coords[3];

  int TopoStatus;
  MPI_Topo_test(Comm, &TopoStatus);
  if (TopoStatus == MPI_CART) {
    MPI_Cart_get(Comm, 3, Dims, Periods, Coords);
  } else {
    Dims[0] = NRanks;
    std::fill(Dims + 1, Dims + 3, 1);
    std::fill(Periods, Periods + 3, 0);
    Coords[0] = Rank;
    std::fill(Coords + 1, Coords + 3, 0);
  }

  std::copy(Coords, Coords + 3, RHLocal.Coords);
  RHLocal.NElems = NElems;
  RHLocal.Start = 0;
  RHLocal.GlobalRank = Rank;

  bool ShouldCompress = DefaultShouldCompress;
  const char *EnvStr = getenv("GENERICIO_COMPRESS");
  if (EnvStr) {
    int Mod = atoi(EnvStr);
    ShouldCompress = (Mod > 0);
  }

  bool NeedsBlockHeaders = ShouldCompress;
  EnvStr = getenv("GENERICIO_FORCE_BLOCKS");
  if (!NeedsBlockHeaders && EnvStr) {
    int Mod = atoi(EnvStr);
    NeedsBlockHeaders = (Mod > 0);
  }

  vector<BlockHeader> LocalBlockHeaders;
  vector<void *> LocalData;
  vector<bool> LocalHasExtraSpace;
  vector<vector<unsigned char> > LocalCData;
  if (NeedsBlockHeaders) {
    LocalBlockHeaders.resize(Vars.size());
    LocalData.resize(Vars.size());
    LocalHasExtraSpace.resize(Vars.size());
    if (ShouldCompress)
      LocalCData.resize(Vars.size());

    for (size_t i = 0; i < Vars.size(); ++i) {
      // Filters null by default, leave null starting address (needs to be
      // calculated by the header-writing rank).
      memset(&LocalBlockHeaders[i], 0, sizeof(BlockHeader));
      if (ShouldCompress) {
        LocalCData[i].resize(sizeof(CompressHeader));

        CompressHeader *CH = (CompressHeader*) &LocalCData[i][0];
        CH->OrigCRC = crc64_omp(Vars[i].Data, Vars[i].Size*NElems);

#ifdef _OPENMP
#pragma omp master
  {
#endif

       if (!blosc_initialized) {
         blosc_init();
         blosc_initialized = true;
       }

#ifdef _OPENMP
       blosc_set_nthreads(omp_get_max_threads());
  }
#endif

        LocalCData[i].resize(LocalCData[i].size() + NElems*Vars[i].Size);
        if (blosc_compress(9, 1, Vars[i].Size, NElems*Vars[i].Size, Vars[i].Data,
                           &LocalCData[i][0] + sizeof(CompressHeader),
                           NElems*Vars[i].Size) <= 0)
          goto nocomp;

        strncpy(LocalBlockHeaders[i].Filters[0], CompressName, FilterNameSize);
        size_t CNBytes, CCBytes, CBlockSize;
        blosc_cbuffer_sizes(&LocalCData[i][0] + sizeof(CompressHeader),
                            &CNBytes, &CCBytes, &CBlockSize);
        LocalCData[i].resize(CCBytes + sizeof(CompressHeader));

        LocalBlockHeaders[i].Size = LocalCData[i].size();
        LocalCData[i].resize(LocalCData[i].size() + CRCSize);
        LocalData[i] = &LocalCData[i][0];
        LocalHasExtraSpace[i] = true;
      } else {
nocomp:
        LocalBlockHeaders[i].Size = NElems*Vars[i].Size;
        LocalData[i] = Vars[i].Data;
        LocalHasExtraSpace[i] = Vars[i].HasExtraSpace;
      }
    }
  }

  double StartTime = MPI_Wtime();

  if (SplitRank == 0) {
    uint64_t HeaderSize = sizeof(GlobalHeader) + Vars.size()*sizeof(VariableHeader) +
                          SplitNRanks*sizeof(RankHeader) + CRCSize;
    if (NeedsBlockHeaders)
      HeaderSize += SplitNRanks*Vars.size()*sizeof(BlockHeader);

    vector<char> Header(HeaderSize, 0);
    GlobalHeader *GH = (GlobalHeader *) &Header[0];
    std::copy(Magic, Magic + MagicSize, GH->Magic);
    GH->HeaderSize = HeaderSize - CRCSize;
    GH->NElems = NElems; // This will be updated later
    std::copy(Dims, Dims + 3, GH->Dims);
    GH->NVars = Vars.size();
    GH->VarsSize = sizeof(VariableHeader);
    GH->VarsStart = sizeof(GlobalHeader);
    GH->NRanks = SplitNRanks;
    GH->RanksSize = sizeof(RankHeader);
    GH->RanksStart = GH->VarsStart + Vars.size()*sizeof(VariableHeader);
    GH->GlobalHeaderSize = sizeof(GlobalHeader);
    std::copy(PhysOrigin, PhysOrigin + 3, GH->PhysOrigin);
    std::copy(PhysScale,  PhysScale  + 3, GH->PhysScale);
    if (!NeedsBlockHeaders) {
      GH->BlocksSize = GH->BlocksStart = 0;
    } else {
      GH->BlocksSize = sizeof(BlockHeader);
      GH->BlocksStart = GH->RanksStart + SplitNRanks*sizeof(RankHeader);
    }

    uint64_t RecordSize = 0;
    VariableHeader *VH = (VariableHeader *) &Header[GH->VarsStart];
    for (size_t i = 0; i < Vars.size(); ++i, ++VH) {
      string VName(Vars[i].Name);
      VName.resize(NameSize);

      std::copy(VName.begin(), VName.end(), VH->Name);
      if (Vars[i].IsFloat)  VH->Flags |= FloatValue;
      if (Vars[i].IsSigned) VH->Flags |= SignedValue;
      if (Vars[i].IsPhysCoordX) VH->Flags |= ValueIsPhysCoordX;
      if (Vars[i].IsPhysCoordY) VH->Flags |= ValueIsPhysCoordY;
      if (Vars[i].IsPhysCoordZ) VH->Flags |= ValueIsPhysCoordZ;
      if (Vars[i].MaybePhysGhost) VH->Flags |= ValueMaybePhysGhost;
      RecordSize += VH->Size = Vars[i].Size;
    }

    MPI_Gather(&RHLocal, sizeof(RHLocal), MPI_BYTE,
               &Header[GH->RanksStart], sizeof(RHLocal),
               MPI_BYTE, 0, SplitComm);

    if (NeedsBlockHeaders) {
      MPI_Gather(&LocalBlockHeaders[0],
                 Vars.size()*sizeof(BlockHeader), MPI_BYTE,
                 &Header[GH->BlocksStart],
                 Vars.size()*sizeof(BlockHeader), MPI_BYTE,
                 0, SplitComm);

      BlockHeader *BH = (BlockHeader *) &Header[GH->BlocksStart];
      for (int i = 0; i < SplitNRanks; ++i)
      for (size_t j = 0; j < Vars.size(); ++j, ++BH) {
        if (i == 0 && j == 0)
          BH->Start = HeaderSize;
        else
          BH->Start = BH[-1].Start + BH[-1].Size + CRCSize;
      }

      RankHeader *RH = (RankHeader *) &Header[GH->RanksStart];
      RH->Start = HeaderSize; ++RH;
      for (int i = 1; i < SplitNRanks; ++i, ++RH) {
        RH->Start =
          ((BlockHeader *) &Header[GH->BlocksStart])[i*Vars.size()].Start;
        GH->NElems += RH->NElems;
      }

      // Compute the total file size.
      uint64_t LastData = BH[-1].Size + CRCSize;
      FileSize = BH[-1].Start + LastData;
    } else {
      RankHeader *RH = (RankHeader *) &Header[GH->RanksStart];
      RH->Start = HeaderSize; ++RH;
      for (int i = 1; i < SplitNRanks; ++i, ++RH) {
        uint64_t PrevNElems = RH[-1].NElems;
        uint64_t PrevData = PrevNElems*RecordSize + CRCSize*Vars.size();
        RH->Start = RH[-1].Start + PrevData;
        GH->NElems += RH->NElems;
      }

      // Compute the total file size.
      uint64_t LastNElems = RH[-1].NElems;
      uint64_t LastData = LastNElems*RecordSize + CRCSize*Vars.size();
      FileSize = RH[-1].Start + LastData;
    }

    // Now that the starting offset has been computed, send it back to each rank.
    MPI_Scatter(&Header[GH->RanksStart], sizeof(RHLocal),
                MPI_BYTE, &RHLocal, sizeof(RHLocal),
                MPI_BYTE, 0, SplitComm);

    if (NeedsBlockHeaders)
      MPI_Scatter(&Header[GH->BlocksStart],
                  sizeof(BlockHeader)*Vars.size(), MPI_BYTE,
                  &LocalBlockHeaders[0],
                  sizeof(BlockHeader)*Vars.size(), MPI_BYTE,
                  0, SplitComm);

    uint64_t HeaderCRC = crc64_omp(&Header[0], HeaderSize - CRCSize);
    crc64_invert(HeaderCRC, &Header[HeaderSize - CRCSize]);

    if (FileIOType == FileIOMPI)
      FH.get() = new GenericFileIO_MPI(MPI_COMM_SELF);
    else if (FileIOType == FileIOMPICollective)
      FH.get() = new GenericFileIO_MPICollective(MPI_COMM_SELF);
    else
      FH.get() = new GenericFileIO_POSIX();

    FH.get()->open(LocalFileName);
    FH.get()->setSize(FileSize);
    FH.get()->write(&Header[0], HeaderSize, 0, "header");

    close();
  } else {
    MPI_Gather(&RHLocal, sizeof(RHLocal), MPI_BYTE, 0, 0, MPI_BYTE, 0, SplitComm);
    if (NeedsBlockHeaders)
      MPI_Gather(&LocalBlockHeaders[0], Vars.size()*sizeof(BlockHeader),
                 MPI_BYTE, 0, 0, MPI_BYTE, 0, SplitComm);
    MPI_Scatter(0, 0, MPI_BYTE, &RHLocal, sizeof(RHLocal), MPI_BYTE, 0, SplitComm);
    if (NeedsBlockHeaders)
      MPI_Scatter(0, 0, MPI_BYTE, &LocalBlockHeaders[0], sizeof(BlockHeader)*Vars.size(),
                  MPI_BYTE, 0, SplitComm);
  }

  MPI_Barrier(SplitComm);

  if (FileIOType == FileIOMPI)
    FH.get() = new GenericFileIO_MPI(SplitComm);
  else if (FileIOType == FileIOMPICollective)
    FH.get() = new GenericFileIO_MPICollective(SplitComm);
  else
    FH.get() = new GenericFileIO_POSIX();

  FH.get()->open(LocalFileName);

  uint64_t Offset = RHLocal.Start;
  for (size_t i = 0; i < Vars.size(); ++i) {
    uint64_t WriteSize = NeedsBlockHeaders ?
                         LocalBlockHeaders[i].Size : NElems*Vars[i].Size;
    void *Data = NeedsBlockHeaders ? LocalData[i] : Vars[i].Data;
    uint64_t CRC = crc64_omp(Data, WriteSize);
    bool HasExtraSpace = NeedsBlockHeaders ?
                         LocalHasExtraSpace[i] : Vars[i].HasExtraSpace;
    char *CRCLoc = HasExtraSpace ?  ((char *) Data) + WriteSize : (char *) &CRC;

    if (NeedsBlockHeaders)
      Offset = LocalBlockHeaders[i].Start;

    // When using extra space for the CRC write, preserve the original contents.
    char CRCSave[CRCSize];
    if (HasExtraSpace)
      std::copy(CRCLoc, CRCLoc + CRCSize, CRCSave);

    crc64_invert(CRC, CRCLoc);

    if (HasExtraSpace) {
      FH.get()->write(Data, WriteSize + CRCSize, Offset, Vars[i].Name + " with CRC");
    } else {
      FH.get()->write(Data, WriteSize, Offset, Vars[i].Name);
      FH.get()->write(CRCLoc, CRCSize, Offset + WriteSize, Vars[i].Name + " CRC");
    }

    if (HasExtraSpace)
       std::copy(CRCSave, CRCSave + CRCSize, CRCLoc);

    Offset += WriteSize + CRCSize;
  }

  close();
  MPI_Barrier(Comm);

  double EndTime = MPI_Wtime();
  double TotalTime = EndTime - StartTime;
  double MaxTotalTime;
  MPI_Reduce(&TotalTime, &MaxTotalTime, 1, MPI_DOUBLE, MPI_MAX, 0, Comm);

  if (SplitNRanks != NRanks) {
    uint64_t ContribFileSize = (SplitRank == 0) ? FileSize : 0;
    MPI_Reduce(&ContribFileSize, &FileSize, 1, MPI_UINT64_T, MPI_SUM, 0, Comm);
  }

  if (Rank == 0) {
    double Rate = ((double) FileSize) / MaxTotalTime / (1024.*1024.);
    cout << "Wrote " << Vars.size() << " variables to " << FileName <<
            " (" << FileSize << " bytes) in " << MaxTotalTime << "s: " <<
            Rate << " MB/s" << endl;
  }

  MPI_Comm_free(&SplitComm);
  SplitComm = MPI_COMM_NULL;
}
#endif // GENERICIO_NO_MPI

// Note: Errors from this function should be recoverable. This means that if
// one rank throws an exception, then all ranks should.
void GenericIO::openAndReadHeader(bool MustMatch, int EffRank, bool CheckPartMap) {
  const char *Magic = isBigEndian() ? MagicBE : MagicLE;
  const char *MagicInv = isBigEndian() ? MagicLE : MagicBE;

  int NRanks, Rank;
#ifndef GENERICIO_NO_MPI
  MPI_Comm_rank(Comm, &Rank);
  MPI_Comm_size(Comm, &NRanks);
#else
  Rank = 0;
  NRanks = 1;
#endif

  if (EffRank == -1)
    EffRank = Rank;

  if (RankMap.empty() && CheckPartMap) {
    // First, check to see if the file is a rank map.
    unsigned long RanksInMap = 0;
    if (Rank == 0) {
      try {
#ifndef GENERICIO_NO_MPI
        GenericIO GIO(MPI_COMM_SELF, FileName, FileIOType);
#else
        GenericIO GIO(FileName, FileIOType);
#endif
        GIO.openAndReadHeader(true, 0, false);
        RanksInMap = GIO.readNumElems();

        RankMap.resize(RanksInMap + GIO.requestedExtraSpace()/sizeof(int));
        GIO.addVariable("$partition", RankMap, true);

        GIO.readData(0, false);
        RankMap.resize(RanksInMap);
      } catch (...) {
        RankMap.clear();
        RanksInMap = 0;
      }
    }

#ifndef GENERICIO_NO_MPI
    MPI_Bcast(&RanksInMap, 1, MPI_UNSIGNED_LONG, 0, Comm);
    if (RanksInMap > 0) {
      RankMap.resize(RanksInMap);
      MPI_Bcast(&RankMap[0], RanksInMap, MPI_INT, 0, Comm);
    }
#endif
  }

#ifndef GENERICIO_NO_MPI
  if (SplitComm != MPI_COMM_NULL)
    MPI_Comm_free(&SplitComm);
#endif

  string LocalFileName;
  if (RankMap.empty()) {
    LocalFileName = FileName;
#ifndef GENERICIO_NO_MPI
    MPI_Comm_dup(Comm, &SplitComm);
#endif
  } else {
    stringstream ss;
    ss << FileName << "#" << RankMap[EffRank];
    LocalFileName = ss.str();
#ifndef GENERICIO_NO_MPI
#ifdef __bgq__
    MPI_Barrier(Comm);
#endif
    MPI_Comm_split(Comm, RankMap[EffRank], Rank, &SplitComm);
#endif
  }

  if (LocalFileName == OpenFileName)
    return;
  FH.close();

  int SplitNRanks, SplitRank;
#ifndef GENERICIO_NO_MPI
  MPI_Comm_rank(SplitComm, &SplitRank);
  MPI_Comm_size(SplitComm, &SplitNRanks);
#else
  SplitRank = 0;
  SplitNRanks = 1;
#endif

  uint64_t HeaderSize;
  vector<char> Header;

  if (SplitRank == 0) {
#ifndef GENERICIO_NO_MPI
    if (FileIOType == FileIOMPI)
      FH.get() = new GenericFileIO_MPI(MPI_COMM_SELF);
    else if (FileIOType == FileIOMPICollective)
      FH.get() = new GenericFileIO_MPICollective(MPI_COMM_SELF);
    else
#endif
      FH.get() = new GenericFileIO_POSIX();

#ifndef GENERICIO_NO_MPI
    char True = 1, False = 0;
#endif

    try {
      FH.get()->open(LocalFileName, true);

      GlobalHeader GH;
      FH.get()->read(&GH, sizeof(GlobalHeader), 0, "global header");

      if (string(GH.Magic, GH.Magic + MagicSize - 1) != Magic) {
        string Error;
        if (string(GH.Magic, GH.Magic + MagicSize - 1) == MagicInv) {
          Error = "wrong endianness";
        } else {
          Error = "invalid file-type identifier";
        }
        throw runtime_error("Won't read " + LocalFileName + ": " + Error);
      }

      if (MustMatch) {
        if (SplitNRanks != (int) GH.NRanks) {
          stringstream ss;
          ss << "Won't read " << LocalFileName << ": communicator-size mismatch: " <<
                "current: " << SplitNRanks << ", file: " << GH.NRanks;
          throw runtime_error(ss.str());
        }

#ifndef GENERICIO_NO_MPI
        int TopoStatus;
        MPI_Topo_test(Comm, &TopoStatus);
        if (TopoStatus == MPI_CART) {
          int Dims[3], Periods[3], Coords[3];
          MPI_Cart_get(Comm, 3, Dims, Periods, Coords);

          bool DimsMatch = true;
          for (int i = 0; i < 3; ++i) {
            if ((uint64_t) Dims[i] != GH.Dims[i]) {
              DimsMatch = false;
              break;
            }
          }

          if (!DimsMatch) {
            stringstream ss;
            ss << "Won't read " << LocalFileName <<
                  ": communicator-decomposition mismatch: " <<
                  "current: " << Dims[0] << "x" << Dims[1] << "x" << Dims[2] <<
                  ", file: " << GH.Dims[0] << "x" << GH.Dims[1] << "x" <<
                  GH.Dims[2];
            throw runtime_error(ss.str());
          }
        }
#endif
      }

      HeaderSize = GH.HeaderSize;
      Header.resize(HeaderSize + CRCSize, 0xFE /* poison */);
      FH.get()->read(&Header[0], HeaderSize + CRCSize, 0, "header");

      uint64_t CRC = crc64_omp(&Header[0], HeaderSize + CRCSize);
      if (CRC != (uint64_t) -1) {
        throw runtime_error("Header CRC check failed: " + LocalFileName);
      }

#ifndef GENERICIO_NO_MPI
      close();
      MPI_Bcast(&True, 1, MPI_BYTE, 0, SplitComm);
#endif
    } catch (...) {
#ifndef GENERICIO_NO_MPI
      MPI_Bcast(&False, 1, MPI_BYTE, 0, SplitComm);
#endif
      close();
      throw;
    }
  } else {
#ifndef GENERICIO_NO_MPI
    char Okay;
    MPI_Bcast(&Okay, 1, MPI_BYTE, 0, SplitComm);
    if (!Okay)
      throw runtime_error("Failure broadcast from rank 0");
#endif
  }

#ifndef GENERICIO_NO_MPI
  MPI_Bcast(&HeaderSize, 1, MPI_UINT64_T, 0, SplitComm);
#endif

  Header.resize(HeaderSize, 0xFD /* poison */);
#ifndef GENERICIO_NO_MPI
  MPI_Bcast(&Header[0], HeaderSize, MPI_BYTE, 0, SplitComm);
#endif

  FH.getHeaderCache().clear();
  FH.getHeaderCache().swap(Header);
  OpenFileName = LocalFileName;

#ifndef GENERICIO_NO_MPI
  MPI_Barrier(Comm);

  if (FileIOType == FileIOMPI)
    FH.get() = new GenericFileIO_MPI(SplitComm);
  else if (FileIOType == FileIOMPICollective)
    FH.get() = new GenericFileIO_MPICollective(SplitComm);
  else
    FH.get() = new GenericFileIO_POSIX();

  int OpenErr = 0, TotOpenErr;
  try {
    FH.get()->open(LocalFileName, true);
    MPI_Allreduce(&OpenErr, &TotOpenErr, 1, MPI_INT, MPI_SUM, Comm);
  } catch (...) {
    OpenErr = 1;
    MPI_Allreduce(&OpenErr, &TotOpenErr, 1, MPI_INT, MPI_SUM, Comm);
    throw;
  }

  if (TotOpenErr > 0) {
    stringstream ss;
    ss << TotOpenErr << " ranks failed to open file: " << LocalFileName;
    throw runtime_error(ss.str());
  }
#endif
}

int GenericIO::readNRanks() {
  if (RankMap.size())
    return RankMap.size();

  assert(FH.getHeaderCache().size() && "HeaderCache must not be empty");
  GlobalHeader *GH = (GlobalHeader *) &FH.getHeaderCache()[0];
  return (int) GH->NRanks;
}

void GenericIO::readDims(int Dims[3]) {
  assert(FH.getHeaderCache().size() && "HeaderCache must not be empty");
  GlobalHeader *GH = (GlobalHeader *) &FH.getHeaderCache()[0];
  std::copy(GH->Dims, GH->Dims + 3, Dims);
}

uint64_t GenericIO::readTotalNumElems() {
  if (RankMap.size())
    return (uint64_t) -1;

  assert(FH.getHeaderCache().size() && "HeaderCache must not be empty");
  GlobalHeader *GH = (GlobalHeader *) &FH.getHeaderCache()[0];
  return GH->NElems;
}

void GenericIO::readPhysOrigin(double Origin[3]) {
  assert(FH.getHeaderCache().size() && "HeaderCache must not be empty");
  GlobalHeader *GH = (GlobalHeader *) &FH.getHeaderCache()[0];
  if (offsetof(GlobalHeader, PhysOrigin) >= GH->GlobalHeaderSize) {
    std::fill(Origin, Origin + 3, 0.0);
    return;
  }

  std::copy(GH->PhysOrigin, GH->PhysOrigin + 3, Origin);
}

void GenericIO::readPhysScale(double Scale[3]) {
  assert(FH.getHeaderCache().size() && "HeaderCache must not be empty");
  GlobalHeader *GH = (GlobalHeader *) &FH.getHeaderCache()[0];
  if (offsetof(GlobalHeader, PhysScale) >= GH->GlobalHeaderSize) {
    std::fill(Scale, Scale + 3, 0.0);
    return;
  }

  std::copy(GH->PhysScale, GH->PhysScale + 3, Scale);
}

static size_t getRankIndex(int EffRank, GlobalHeader *GH,
                           vector<int> &RankMap, vector<char> &HeaderCache) {
  if (offsetof(RankHeader, GlobalRank) >= GH->RanksSize || RankMap.empty())
    return EffRank;

  for (size_t i = 0; i < GH->NRanks; ++i) {
    RankHeader *RH = (RankHeader *) &HeaderCache[GH->RanksStart +
                                                 i*GH->RanksSize];
    if ((int) RH->GlobalRank == EffRank)
      return i;
  }

  assert(false && "Index requested of an invalid rank");
  return (size_t) -1;
}

int GenericIO::readGlobalRankNumber(int EffRank) {
  if (EffRank == -1) {
#ifndef GENERICIO_NO_MPI
    MPI_Comm_rank(Comm, &EffRank);
#else
    EffRank = 0;
#endif
  }

  openAndReadHeader(false, EffRank, false);

  assert(FH.getHeaderCache().size() && "HeaderCache must not be empty");

  GlobalHeader *GH = (GlobalHeader *) &FH.getHeaderCache()[0];
  size_t RankIndex = getRankIndex(EffRank, GH, RankMap, FH.getHeaderCache());

  assert(RankIndex < GH->NRanks && "Invalid rank specified");

  if (offsetof(RankHeader, GlobalRank) >= GH->RanksSize)
    return EffRank;

  RankHeader *RH = (RankHeader *) &FH.getHeaderCache()[GH->RanksStart +
                                               RankIndex*GH->RanksSize];

  return (int) RH->GlobalRank;
}

size_t GenericIO::readNumElems(int EffRank) {
  if (EffRank == -1) {
#ifndef GENERICIO_NO_MPI
    MPI_Comm_rank(Comm, &EffRank);
#else
    EffRank = 0;
#endif
  }

  openAndReadHeader(false, EffRank, false);

  assert(FH.getHeaderCache().size() && "HeaderCache must not be empty");

  GlobalHeader *GH = (GlobalHeader *) &FH.getHeaderCache()[0];
  size_t RankIndex = getRankIndex(EffRank, GH, RankMap, FH.getHeaderCache());

  assert(RankIndex < GH->NRanks && "Invalid rank specified");

  RankHeader *RH = (RankHeader *) &FH.getHeaderCache()[GH->RanksStart +
                                               RankIndex*GH->RanksSize];
  return (size_t) RH->NElems;
}

void GenericIO::readCoords(int Coords[3], int EffRank) {
  if (EffRank == -1) {
#ifndef GENERICIO_NO_MPI
    MPI_Comm_rank(Comm, &EffRank);
#else
    EffRank = 0;
#endif
  }

  openAndReadHeader(false, EffRank, false);

  assert(FH.getHeaderCache().size() && "HeaderCache must not be empty");

  GlobalHeader *GH = (GlobalHeader *) &FH.getHeaderCache()[0];
  size_t RankIndex = getRankIndex(EffRank, GH, RankMap, FH.getHeaderCache());

  assert(RankIndex < GH->NRanks && "Invalid rank specified");

  RankHeader *RH = (RankHeader *) &FH.getHeaderCache()[GH->RanksStart +
                                               RankIndex*GH->RanksSize];

  std::copy(RH->Coords, RH->Coords + 3, Coords);
}

// Note: Errors from this function should be recoverable. This means that if
// one rank throws an exception, then all ranks should.
void GenericIO::readData(int EffRank, bool PrintStats, bool CollStats) {
  int Rank;
#ifndef GENERICIO_NO_MPI
  MPI_Comm_rank(Comm, &Rank);
#else
  Rank = 0;
#endif

  openAndReadHeader(false, EffRank, false);

  assert(FH.getHeaderCache().size() && "HeaderCache must not be empty");

  if (EffRank == -1)
    EffRank = Rank;

  GlobalHeader *GH = (GlobalHeader *) &FH.getHeaderCache()[0];
  size_t RankIndex = getRankIndex(EffRank, GH, RankMap, FH.getHeaderCache());

  assert(RankIndex < GH->NRanks && "Invalid rank specified");

  RankHeader *RH = (RankHeader *) &FH.getHeaderCache()[GH->RanksStart +
                                               RankIndex*GH->RanksSize];

  uint64_t TotalReadSize = 0;
#ifndef GENERICIO_NO_MPI
  double StartTime = MPI_Wtime();
#else
  double StartTime = double(clock())/CLOCKS_PER_SEC;
#endif

  int NErrs[3] = { 0, 0, 0 };
  for (size_t i = 0; i < Vars.size(); ++i) {
    uint64_t Offset = RH->Start;
    bool VarFound = false;
    for (uint64_t j = 0; j < GH->NVars; ++j) {
      VariableHeader *VH = (VariableHeader *) &FH.getHeaderCache()[GH->VarsStart +
                                                           j*GH->VarsSize];

      string VName(VH->Name, VH->Name + NameSize);
      size_t VNameNull = VName.find('\0');
      if (VNameNull < NameSize)
        VName.resize(VNameNull);

      uint64_t ReadSize = RH->NElems*VH->Size + CRCSize;
      if (VName != Vars[i].Name) {
        Offset += ReadSize;
        continue;
      }

      VarFound = true;
      bool IsFloat = (bool) (VH->Flags & FloatValue),
           IsSigned = (bool) (VH->Flags & SignedValue);
      if (VH->Size != Vars[i].Size) {
        stringstream ss;
        ss << "Size mismatch for variable " << Vars[i].Name <<
              " in: " << OpenFileName << ": current: " << Vars[i].Size <<
              ", file: " << VH->Size;
        throw runtime_error(ss.str());
      } else if (IsFloat != Vars[i].IsFloat) {
        string Float("float"), Int("integer");
        stringstream ss;
        ss << "Type mismatch for variable " << Vars[i].Name <<
              " in: " << OpenFileName << ": current: " <<
              (Vars[i].IsFloat ? Float : Int) <<
              ", file: " << (IsFloat ? Float : Int);
        throw runtime_error(ss.str());
      } else if (IsSigned != Vars[i].IsSigned) {
        string Signed("signed"), Uns("unsigned");
        stringstream ss;
        ss << "Type mismatch for variable " << Vars[i].Name <<
              " in: " << OpenFileName << ": current: " <<
              (Vars[i].IsSigned ? Signed : Uns) <<
              ", file: " << (IsSigned ? Signed : Uns);
        throw runtime_error(ss.str());
      }

      vector<unsigned char> LData;
      void *Data = Vars[i].Data;
      bool HasExtraSpace = Vars[i].HasExtraSpace;
      if (offsetof(GlobalHeader, BlocksStart) < GH->GlobalHeaderSize &&
          GH->BlocksSize > 0) {
        BlockHeader *BH = (BlockHeader *)
          &FH.getHeaderCache()[GH->BlocksStart +
                               (RankIndex*GH->NVars + j)*GH->BlocksSize];
        ReadSize = BH->Size + CRCSize;
        Offset = BH->Start;

        if (strncmp(BH->Filters[0], CompressName, FilterNameSize) == 0) {
          LData.resize(ReadSize);
          Data = &LData[0];
          HasExtraSpace = true;
        } else if (BH->Filters[0][0] != '\0') {
          stringstream ss;
          ss << "Unknown filter \"" << BH->Filters[0] << "\" on variable " << Vars[i].Name;
          throw runtime_error(ss.str());
        }
      }

      assert(HasExtraSpace && "Extra space required for reading");

      char CRCSave[CRCSize];
      char *CRCLoc = ((char *) Data) + ReadSize - CRCSize;
      if (HasExtraSpace)
        std::copy(CRCLoc, CRCLoc + CRCSize, CRCSave);

      int Retry = 0;
      {
        int RetryCount = 300;
        const char *EnvStr = getenv("GENERICIO_RETRY_COUNT");
        if (EnvStr)
          RetryCount = atoi(EnvStr);

        int RetrySleep = 100; // ms
        EnvStr = getenv("GENERICIO_RETRY_SLEEP");
        if (EnvStr)
          RetrySleep = atoi(EnvStr);

        for (; Retry < RetryCount; ++Retry) {
          try {
            FH.get()->read(Data, ReadSize, Offset, Vars[i].Name);
            break;
          } catch (...) { }

          usleep(1000*RetrySleep);
        }

        if (Retry == RetryCount) {
          ++NErrs[0];
          break;
        } else if (Retry > 0) {
          EnvStr = getenv("GENERICIO_VERBOSE");
          if (EnvStr) {
            int Mod = atoi(EnvStr);
            if (Mod > 0) {
              int Rank;
#ifndef GENERICIO_NO_MPI
              MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
#else
              Rank = 0;
#endif

              std::cerr << "Rank " << Rank << ": " << Retry <<
                           " I/O retries were necessary for reading " <<
                           Vars[i].Name << " from: " << OpenFileName << "\n";

              std::cerr.flush();
            }
          }
        }
      }

      TotalReadSize += ReadSize;

      uint64_t CRC = crc64_omp(Data, ReadSize);
      if (CRC != (uint64_t) -1) {
        ++NErrs[1];

        int Rank;
#ifndef GENERICIO_NO_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
#else
        Rank = 0;
#endif

        // All ranks will do this and have a good time!
        string dn = "gio_crc_errors";
        mkdir(dn.c_str(), 0777);

        srand(time(0));
        int DumpNum = rand();
        stringstream ssd;
        ssd << dn << "/gio_crc_error_dump." << Rank << "." << DumpNum << ".bin";

        stringstream ss;
        ss << dn << "/gio_crc_error_log." << Rank << ".txt";

        ofstream ofs(ss.str().c_str(), ofstream::out | ofstream::app);
        ofs << "On-Disk CRC Error Report:\n";
        ofs << "Variable: " << Vars[i].Name << "\n";
        ofs << "File: " << OpenFileName << "\n";
        ofs << "I/O Retries: " << Retry << "\n";
        ofs << "Size: " << ReadSize << " bytes\n";
        ofs << "Offset: " << Offset << " bytes\n";
        ofs << "CRC: " << CRC << " (expected is -1)\n";
        ofs << "Dump file: " << ssd.str() << "\n";
        ofs << "\n";
        ofs.close();

        ofstream dofs(ssd.str().c_str(), ofstream::out);
        dofs.write((const char *) Data, ReadSize);
        dofs.close();
        break;
      }

      if (HasExtraSpace)
        std::copy(CRCSave, CRCSave + CRCSize, CRCLoc);

      if (LData.size()) {
        CompressHeader *CH = (CompressHeader*) &LData[0];

#ifdef _OPENMP
#pragma omp master
  {
#endif

       if (!blosc_initialized) {
         blosc_init();
         blosc_initialized = true;
       }

#ifdef _OPENMP
       blosc_set_nthreads(omp_get_max_threads());
  }
#endif

        blosc_decompress(&LData[0] + sizeof(CompressHeader),
                         Vars[i].Data, Vars[i].Size*RH->NElems);

        if (CH->OrigCRC != crc64_omp(Vars[i].Data, Vars[i].Size*RH->NElems)) {
          ++NErrs[2];
          break;
        }
      }

      break;
    }

    if (!VarFound)
      throw runtime_error("Variable " + Vars[i].Name +
                          " not found in: " + OpenFileName);

    // This is for debugging.
    if (NErrs[0] || NErrs[1] || NErrs[2]) {
      const char *EnvStr = getenv("GENERICIO_VERBOSE");
      if (EnvStr) {
        int Mod = atoi(EnvStr);
        if (Mod > 0) {
          int Rank;
#ifndef GENERICIO_NO_MPI
          MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
#else
          Rank = 0;
#endif

          std::cerr << "Rank " << Rank << ": " << NErrs[0] << " I/O error(s), " <<
          NErrs[1] << " CRC error(s) and " << NErrs[2] <<
          " decompression CRC error(s) reading: " << Vars[i].Name <<
          " from: " << OpenFileName << "\n";

          std::cerr.flush();
        }
      }
    }

    if (NErrs[0] || NErrs[1] || NErrs[2])
      break;
  }

  int AllNErrs[3];
#ifndef GENERICIO_NO_MPI
  MPI_Allreduce(NErrs, AllNErrs, 3, MPI_INT, MPI_SUM, Comm);
#else
  AllNErrs[0] = NErrs[0]; AllNErrs[1] = NErrs[1]; AllNErrs[2] = NErrs[2];
#endif

  if (AllNErrs[0] > 0 || AllNErrs[1] > 0 || AllNErrs[2] > 0) {
    stringstream ss;
    ss << "Experienced " << AllNErrs[0] << " I/O error(s), " <<
          AllNErrs[1] << " CRC error(s) and " << AllNErrs[2] <<
          " decompression CRC error(s) reading: " << OpenFileName;
    throw runtime_error(ss.str());
  }

#ifndef GENERICIO_NO_MPI
  MPI_Barrier(Comm);
#endif

#ifndef GENERICIO_NO_MPI
  double EndTime = MPI_Wtime();
#else
  double EndTime = double(clock())/CLOCKS_PER_SEC;
#endif

  double TotalTime = EndTime - StartTime;
  double MaxTotalTime;
#ifndef GENERICIO_NO_MPI
  if (CollStats)
    MPI_Reduce(&TotalTime, &MaxTotalTime, 1, MPI_DOUBLE, MPI_MAX, 0, Comm);
  else
#endif
  MaxTotalTime = TotalTime;

  uint64_t AllTotalReadSize;
#ifndef GENERICIO_NO_MPI
  if (CollStats)
    MPI_Reduce(&TotalReadSize, &AllTotalReadSize, 1, MPI_UINT64_T, MPI_SUM, 0, Comm);
  else
#endif
  AllTotalReadSize = TotalReadSize;

  if (Rank == 0 && PrintStats) {
    double Rate = ((double) AllTotalReadSize) / MaxTotalTime / (1024.*1024.);
    cout << "Read " << Vars.size() << " variables from " << FileName <<
            " (" << AllTotalReadSize << " bytes) in " << MaxTotalTime << "s: " <<
            Rate << " MB/s [excluding header read]" << endl;
  }

}

void GenericIO::getVariableInfo(vector<VariableInfo> &VI) {
  assert(FH.getHeaderCache().size() && "HeaderCache must not be empty");

  GlobalHeader *GH = (GlobalHeader *) &FH.getHeaderCache()[0];
  for (uint64_t j = 0; j < GH->NVars; ++j) {
    VariableHeader *VH = (VariableHeader *) &FH.getHeaderCache()[GH->VarsStart +
                                                         j*GH->VarsSize];

    string VName(VH->Name, VH->Name + NameSize);
    size_t VNameNull = VName.find('\0');
    if (VNameNull < NameSize)
      VName.resize(VNameNull);

    bool IsFloat = (bool) (VH->Flags & FloatValue),
         IsSigned = (bool) (VH->Flags & SignedValue),
         IsPhysCoordX = (bool) (VH->Flags & ValueIsPhysCoordX),
         IsPhysCoordY = (bool) (VH->Flags & ValueIsPhysCoordY),
         IsPhysCoordZ = (bool) (VH->Flags & ValueIsPhysCoordZ),
         MaybePhysGhost = (bool) (VH->Flags & ValueMaybePhysGhost);
    VI.push_back(VariableInfo(VName, (size_t) VH->Size, IsFloat, IsSigned,
                              IsPhysCoordX, IsPhysCoordY, IsPhysCoordZ,
                              MaybePhysGhost));
  }
}

void GenericIO::setNaturalDefaultPartition() {
#ifdef __bgq__
  Personality_t pers;
  Kernel_GetPersonality(&pers, sizeof(pers));

  // Nodes in an ION Partition
  int SPLIT_A = 2;
  int SPLIT_B = 2;
  int SPLIT_C = 4;
  int SPLIT_D = 4;
  int SPLIT_E = 2;

  int Anodes, Bnodes, Cnodes, Dnodes, Enodes;
  int Acoord, Bcoord, Ccoord, Dcoord, Ecoord;
  int A_color, B_color, C_color, D_color, E_color;
  int A_blocks, B_blocks, C_blocks, D_blocks, E_blocks;
  uint32_t id_on_node;
  int ranks_per_node, color;

  Anodes = pers.Network_Config.Anodes;
  Acoord = pers.Network_Config.Acoord;

  Bnodes = pers.Network_Config.Bnodes;
  Bcoord = pers.Network_Config.Bcoord;

  Cnodes = pers.Network_Config.Cnodes;
  Ccoord = pers.Network_Config.Ccoord;

  Dnodes = pers.Network_Config.Dnodes;
  Dcoord = pers.Network_Config.Dcoord;

  Enodes = pers.Network_Config.Enodes;
  Ecoord = pers.Network_Config.Ecoord;

  A_color  = Acoord /  SPLIT_A;
  B_color  = Bcoord /  SPLIT_B;
  C_color  = Ccoord /  SPLIT_C;
  D_color  = Dcoord /  SPLIT_D;
  E_color  = Ecoord /  SPLIT_E;

  // Number of blocks
  A_blocks = Anodes / SPLIT_A;
  B_blocks = Bnodes / SPLIT_B;
  C_blocks = Cnodes / SPLIT_C;
  D_blocks = Dnodes / SPLIT_D;
  E_blocks = Enodes / SPLIT_E;

  color = (A_color * (B_blocks * C_blocks * D_blocks * E_blocks))
    + (B_color * (C_blocks * D_blocks * E_blocks))
    + (C_color * ( D_blocks * E_blocks))
    + (D_color * ( E_blocks))
    + E_color;

  DefaultPartition = color;
#else
#ifndef GENERICIO_NO_MPI
  bool UseName = true;
  const char *EnvStr = getenv("GENERICIO_PARTITIONS_USE_NAME");
  if (EnvStr) {
    int Mod = atoi(EnvStr);
    UseName = (Mod != 0);
  }

  if (UseName) {
    // This is a heuristic to generate ~256 partitions based on the
    // names of the nodes.
    char Name[MPI_MAX_PROCESSOR_NAME];
    int Len = 0;

    MPI_Get_processor_name(Name, &Len);
    unsigned char color = 0;
    for (int i = 0; i < Len; ++i)
      color += (unsigned char) Name[i];

    DefaultPartition = color;
  }

  // This is for debugging.
  EnvStr = getenv("GENERICIO_RANK_PARTITIONS");
  if (EnvStr) {
    int Mod = atoi(EnvStr);
    if (Mod > 0) {
      int Rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
      DefaultPartition += Rank % Mod;
    }
  }
#endif
#endif
}

} /* END namespace cosmotk */
