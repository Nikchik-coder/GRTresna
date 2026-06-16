/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#include "mpi.h"
#include <iostream>

#include "CTTK.hpp"
#include "CTTKHybrid.hpp"
#include "ComplexScalarField.hpp"
#include "GRParmParse.hpp"
#include "GRSolver.hpp"

using namespace std;

int main(int argc, char *argv[])
{
    int status = 0;
#ifdef CH_MPI
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
        cout << "Running with MPI" << endl;
#endif

    if (argc < 1)
    {
        cerr << " usage " << argv[0] << " <input_file_name> " << endl;
        exit(0);
    }

    char *inFile = argv[1];
    GRParmParse pp(argc - 2, argv + 2, NULL, inFile);

#if defined(USE_CTTKHybrid)
    GRSolver<CTTKHybrid<ComplexScalarField>, ComplexScalarField> solver(pp);
    pout() << "Using CTTK Hybrid with ComplexScalarField (boson star)" << endl;
#else
#error "BosonStarBH requires USE_CTTKHybrid"
#endif

    solver.setup();
    status = solver.run();

#ifdef CH_MPI
    MPI_Finalize();
#endif
    return status;
}
