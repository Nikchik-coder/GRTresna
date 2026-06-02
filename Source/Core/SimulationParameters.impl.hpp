/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#error "This file should only be included through SimulationParameters.hpp"
#endif

#include "CoarseAverage.H"
#include "FilesystemTools.hpp"
#include "GRParmParse.hpp"
#include "ProblemDomain.H"
#include "REAL.H"
#include "RealVect.H"

#include <fstream>
#include <iostream>

template <class method_t, class matter_t>
struct SimulationParameters<method_t, matter_t>::BaseParams
{

    int max_NL_iter;
    bool write_diagnostics;
    int diagnostic_interval;
    Real iter_tolerance;
    int max_iter;
    int numMGIter;
    int numMGSmooth;
    int preCondSolverDepth;
    Real alpha;
    Real beta;
    bool readin_matter_data;
    std::string input_filename;
    std::string output_filename;
    std::string output_path;
    std::string pout_path;
    std::string error_filename;
    Real G_Newton;
    int verbosity;
    // Outer (nonlinear) Newton-iteration robustness, used mainly for exotic
    // (negative-energy) matter where the linearised psi-operator can become
    // indefinite and a full Newton step overshoots psi into the unphysical
    // psi <= 0 region (pow(psi, -7) -> NaN). psi_relaxation in (0, 1] damps the
    // applied correction; psi_floor (> 0) clamps psi_reg from below so psi
    // stays positive. Defaults reproduce the original undamped full-step
    // behaviour.
    Real psi_relaxation;
    Real psi_floor;
};

template <class method_t, class matter_t>
void SimulationParameters<method_t, matter_t>::read_base_params(GRParmParse &pp)
{
#ifdef CH_MPI
    // setPoutBaseName must be called before pout is used
    if (pp.contains("pout_path"))
    {
        pp.get("pout_path", base_params.pout_path);
        // Create pout directory
        if (!FilesystemTools::directory_exists(base_params.pout_path))
            FilesystemTools::mkdir_recursive(base_params.pout_path);
        setPoutBaseName(base_params.pout_path + "pout");
    }
    else
    {
        base_params.pout_path = "pout/";
        if (!FilesystemTools::directory_exists(base_params.pout_path))
            FilesystemTools::mkdir_recursive(base_params.pout_path);
        setPoutBaseName(base_params.pout_path + "pout");
    }
#endif

    pp.load("max_NL_iterations", base_params.max_NL_iter, 100);
    pp.load("write_diagnostics", base_params.write_diagnostics, true);
    pp.load("diagnostic_interval", base_params.diagnostic_interval, 10);

    // Setup multigrid params, most of them defaulted

    pp.load("iter_tolerance", base_params.iter_tolerance, 1e-10);
    pp.load("max_iterations", base_params.max_iter, 100);
    pp.load("numMGIterations", base_params.numMGIter, 4);
    pp.load("numMGsmooth", base_params.numMGSmooth, 4);
    pp.load("preCondSolverDepth", base_params.preCondSolverDepth, -1);

    // Params for variable coefficient multigrid solver, solving the eqn
    // alpha*aCoef(x)*I - beta*bCoef(x) * laplacian = rhs
    // spatially-varying aCoef and bCoef are set in Methods
    // (for pure laplacian, alpha = 0, beta=-1)
    pp.load("alpha", base_params.alpha, 1.0);
    pp.load("beta", base_params.beta, -1.0);

    // Read from hdf5 file
    if (pp.contains("input_filename"))
    {
        pp.get("input_filename", base_params.input_filename);
        base_params.readin_matter_data = true;
    }
    else
    {
        base_params.input_filename = "";
        base_params.readin_matter_data = false;
    }

    // Error outputs
    if (pp.contains("error_filename"))
    {
        pp.get("error_filename", base_params.error_filename);
    }
    else
    {
        base_params.error_filename = "Ham_and_Mom_errors";
    }

    // Output for final file
    if (pp.contains("output_path"))
    {
        pp.get("output_path", base_params.output_path);
        if (!FilesystemTools::directory_exists(base_params.output_path))
            FilesystemTools::mkdir_recursive(base_params.output_path);
    }
    else
    {
        base_params.output_path = "";
    }

    if (pp.contains("output_filename"))
    {
        string filename;
        pp.get("output_filename", filename);
        base_params.output_filename = base_params.output_path + filename;
    }
    else
    {
        base_params.output_filename =
            base_params.output_path + "InitialDataFinal.3d.hdf5";
    }

    pp.load("G_Newton", base_params.G_Newton, 1.0);
    pp.load("verbosity", base_params.verbosity, 1);

    // Nonlinear-iteration robustness (default = original undamped behaviour).
    pp.load("psi_relaxation", base_params.psi_relaxation, 1.0);
    pp.load("psi_floor", base_params.psi_floor, -1.0);
}
