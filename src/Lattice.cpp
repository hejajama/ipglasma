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

    posmX.reserve(size);
    pospX.reserve(size);
    posmY.reserve(size);
    pospY.reserve(size);
    posmXpY.reserve(size);
    pospXmY.reserve(size);

    for (int i = 0; i < length; i++) {
        const int im = std::max(0, i - 1);
        const int ip = std::min(length - 1, i + 1);
        for (int j = 0; j < length; j++) {
            const int jm = std::max(0, j - 1);
            const int jp = std::min(length - 1, j + 1);

            pospX.push_back(ip * length + j);
            pospY.push_back(i * length + jp);
            posmX.push_back(im * length + j);
            posmY.push_back(i * length + jm);

            // Diagonal neighbors used in evolveE() and Tmunu().  These use
            // the same clamped boundary convention as the original inline
            // index calculations.
            posmXpY.push_back(im * length + jp);
            pospXmY.push_back(ip * length + jm);
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
