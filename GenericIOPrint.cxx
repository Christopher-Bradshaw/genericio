#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <stdint.h>

#include "GenericIO.h"

using namespace std;
using namespace gio;

class PrinterBase {
public:
  virtual ~PrinterBase() {}
  virtual void print(ostream &os, size_t i) = 0;
};

template <class T>
class Printer : public PrinterBase {
public:
  Printer(GenericIO &G, size_t MNE, const string &N)
    : Data(MNE + G.requestedExtraSpace()/sizeof(T)) {
    G.addVariable(N, Data, true);
  }

  virtual void print(ostream &os, size_t i) {
    os << scientific << setprecision(numeric_limits<T>::digits10) << Data[i];
  }

protected:
  vector<T> Data;
};

template <typename T>
PrinterBase *addPrinter(GenericIO::VariableInfo &V,
                GenericIO &GIO, size_t MNE) {
  if (sizeof(T) != V.Size)
    return 0;

  if (V.IsFloat == numeric_limits<T>::is_integer)
    return 0;
  if (V.IsSigned != numeric_limits<T>::is_signed)
    return 0;

  return new Printer<T>(GIO, MNE, V.Name);
}

int main(int argc, char *argv[]) {
#ifndef GENERICIO_NO_MPI
  MPI_Init(&argc, &argv);
#endif

  bool ShowMap = false;
  bool NoData = false;
  bool PrintRankInfo = true;
  int FileNameIdx = 1;
  if (argc > 2) {
    if (string(argv[1]) == "--no-rank-info") {
      PrintRankInfo = false;
      ++FileNameIdx;
      --argc;
    } else if (string(argv[1]) == "--no-data") {
      NoData = true;
      ++FileNameIdx;
      --argc;
    } else if (string(argv[1]) == "--show-map") {
      ShowMap = true;
      ++FileNameIdx;
      --argc;
    }
  }

  if (argc != 2) {
    cerr << "Usage: " << argv[0] << " [--no-rank-info|--no-data|--show-map] <mpiioName>" << endl;
    exit(-1);
  }

  string FileName(argv[FileNameIdx]);

  int Rank, NRanks;
#ifndef GENERICIO_NO_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
  MPI_Comm_size(MPI_COMM_WORLD, &NRanks);
#else
  Rank = 0;
  NRanks = 1;
#endif

  if (Rank == 0) {
    unsigned Method = GenericIO::FileIOPOSIX;
#ifndef GENERICIO_NO_MPI
    const char *EnvStr = getenv("GENERICIO_USE_MPIIO");
    if (EnvStr && string(EnvStr) == "1")
      Method = GenericIO::FileIOMPI;
#endif

#ifndef GENERICIO_NO_MPI
    GenericIO GIO(MPI_COMM_SELF, FileName, Method);
#else
    GenericIO GIO(FileName, Method);
#endif
    GIO.openAndReadHeader(false, -1, !ShowMap);

    int NR = GIO.readNRanks();

    vector<GenericIO::VariableInfo> VI;
    GIO.getVariableInfo(VI);

    size_t MaxNElem = 0;
    for (int i = 0; i < NR; ++i) {
      size_t NElem = GIO.readNumElems(i);
      MaxNElem = max(MaxNElem, NElem);
    }

    vector<PrinterBase *> Printers;
    for (size_t i = 0; i < VI.size(); ++i) {
      PrinterBase *P = 0;

#define ADD_PRINTER(T) \
      if (!P) P = addPrinter<T>(VI[i], GIO, MaxNElem)
      ADD_PRINTER(float);
      ADD_PRINTER(double);
      ADD_PRINTER(unsigned char);
      ADD_PRINTER(signed char);
      ADD_PRINTER(int16_t);
      ADD_PRINTER(uint16_t);
      ADD_PRINTER(int32_t);
      ADD_PRINTER(uint32_t);
      ADD_PRINTER(int64_t);
      ADD_PRINTER(uint64_t);
#undef ADD_PRINTER 

      if (!P) throw runtime_error("Don't know how to print variable: " + VI[i].Name);
      Printers.push_back(P);
    }

    int Dims[3];
    GIO.readDims(Dims);

    cout << "# " << FileName << ": " << NR << " rank(s): " <<
            Dims[0] << "x" << Dims[1] << "x" << Dims[2];

    uint64_t TotalNumElems = GIO.readTotalNumElems();
    if (TotalNumElems != (uint64_t) -1)
      cout  << ": " << GIO.readTotalNumElems() << " row(s)";
    cout << endl;

    double PhysOrigin[3], PhysScale[3];
    GIO.readPhysOrigin(PhysOrigin);
    GIO.readPhysScale(PhysScale);
    if (PhysScale[0] != 0.0 || PhysScale[1] != 0.0 || PhysScale[2] != 0.0) {
      cout << "# physical coordinates: (" << PhysOrigin[0] << "," <<
              PhysOrigin[1] << "," << PhysOrigin[2] << ") -> (" <<
              PhysScale[0] << "," << PhysScale[1] << "," <<
              PhysScale[2] << ")" << endl;


      vector<GenericIO::VariableInfo> VIX, VIY, VIZ;
      for (size_t i = 0; i < VI.size(); ++i) {
        if (VI[i].IsPhysCoordX) VIX.push_back(VI[i]);
        if (VI[i].IsPhysCoordY) VIY.push_back(VI[i]);
        if (VI[i].IsPhysCoordZ) VIZ.push_back(VI[i]);
      }

      if (!VIX.empty()) {
        cout << "# x variables: ";
        for (size_t i = 0; i < VIX.size(); ++i) {
          if (i != 0) cout << ", ";
          cout << VIX[i].Name;
          if (VIX[i].MaybePhysGhost)
            cout << " [maybe ghost]";
        }
        cout << endl;
      }
      if (!VIY.empty()) {
        cout << "# y variables: ";
        for (size_t i = 0; i < VIY.size(); ++i) {
          if (i != 0) cout << ", ";
          cout << VIY[i].Name;
          if (VIY[i].MaybePhysGhost)
            cout << " [maybe ghost]";
        }
        cout << endl;
      }
      if (!VIZ.empty()) {
        cout << "# z variables: ";
        for (size_t i = 0; i < VIZ.size(); ++i) {
          if (i != 0) cout << ", ";
          cout << VIZ[i].Name;
          if (VIZ[i].MaybePhysGhost)
            cout << " [maybe ghost]";
        }
        cout << endl;
      }
    }

    cout << "# ";
    for (size_t i = 0; i < VI.size(); ++i) {
      cout << VI[i].Name;
      if (i != VI.size() - 1)
        cout << "\t";
    }
    cout << endl;

    for (int i = 0; i < NR; ++i) {
      size_t NElem = GIO.readNumElems(i);
      int Coords[3];
      GIO.readCoords(Coords, i);
      if (PrintRankInfo)
        cout << "# rank " << GIO.readGlobalRankNumber(i) << ": " <<
                Coords[0] << "," << Coords[1] << "," <<
                Coords[2] << ": " << NElem << " row(s)" << endl;
      if (NoData)
        continue;

      GIO.readData(i, false);
      for (size_t j = 0; j < NElem; ++j) {
        for (size_t k = 0; k < Printers.size(); ++k) {
          Printers[k]->print(cout, j);
          if (k != Printers.size() - 1)
            cout << "\t";
        }
        cout << endl;
      }
    }

    for (size_t i = 0; i < Printers.size(); ++i) {
      delete Printers[i];
    }
  }

#ifndef GENERICIO_NO_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
#endif

  return 0;
}

