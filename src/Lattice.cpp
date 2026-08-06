#include "Lattice.h"
#include "Instrumentation.h"

// constructor
Lattice::Lattice(Parameters *param, int N, int length) {
    IPG_PROFILE_SCOPE("lattice.allocate");
    Nc = N;
    size = length * length;
    double a = param->getL() / static_cast<double>(length);

    std::cout << "Allocating square lattice of size " << length << "x" << length
              << " with a=" << a << " fm ...";

    // initialize the array of cells in one contiguous allocation
    cellStorage.reserve(size);
    cells.reserve(size);
    for (int i = 0; i < size; i++) cellStorage.emplace_back(Nc);
    for (int i = 0; i < size; i++) cells.push_back(&cellStorage[i]);

    for (int i = 0; i < length; i++) {
        for (int j = 0; j < length; j++) {
            // pos = i*length+j;
            pospX.push_back((std::min(length - 1, i + 1)) * length + j);
            pospY.push_back(i * length + std::min(length - 1, j + 1));

            posmX.push_back((std::max(0, i - 1)) * length + j);
            posmY.push_back(i * length + std::max(0, j - 1));
        }
    }
    std::cout << " done on rank " << param->getMPIRank() << "." << std::endl;
}

Lattice::~Lattice() { cells.clear(); }

// constructor
BufferLattice::BufferLattice(int N, int length) {
    Nc = N;
    size = length * length;

    cellStorage.reserve(size);
    cells.reserve(size);
    for (int i = 0; i < size; i++) cellStorage.emplace_back(Nc);
    for (int i = 0; i < size; i++) cells.push_back(&cellStorage[i]);
}

BufferLattice::~BufferLattice() { cells.clear(); }
