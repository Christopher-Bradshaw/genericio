/* With the current master version of GIO I think that this hangs.
 * With the code in this branch I think it works. */
#include "../GenericIO.h"
#include <vector>

using namespace std;
using namespace gio;

int main() {
    MPI_Init(NULL, NULL);

    int rank, nRanks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nRanks);

    unsigned method = GenericIO::FileIOMPI;
    string filename = "test_file";

    GenericIO* GIO = new GenericIO(MPI_COMM_WORLD, filename, method);

    vector<int> data;
    for (int i = 0; i < 10; i++) {
        data.push_back(i * (rank+1));
    }

    GIO->addScalarizedVariable("test", data, 1, 0);
    GIO->setNumElems(data.size());
    GIO->write();


    // READ
    GIO->clearVariables();

    GIO->openAndReadHeader(GenericIO::MismatchAllowed, -1, true);
    vector<GenericIO::VariableInfo> vi;
    GIO->getVariableInfo(vi);

    for (int i = 0; i < int(vi.size()); i++) {
        cout << vi[i].Name << " " << vi[i].Size << " " << vi[i].ElementSize << endl;
    }

    vector<int> readData;
    for (int i = 0; i < int(data.size() + GIO->requestedExtraSpace()); i++) {
        readData.push_back(0);
    }
    cout << readData.size() << endl;
    GIO->addScalarizedVariable("test", readData, 1, GenericIO::VarHasExtraSpace);
    GIO->readData(rank, false, false);

    for (int i = 0; i < int(data.size()); i++) {
        cout << readData[i] << endl;
    }


    delete GIO;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
