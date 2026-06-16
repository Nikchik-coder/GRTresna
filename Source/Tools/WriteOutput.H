/* GRTresna                                                                          
* Copyright 2024 The GRTL collaboration.
* Please refer to LICENSE in GRTresna's root directory.
*/

#ifndef _WRITEOUTPUT_H_
#define _WRITEOUTPUT_H_

#include "AMRIO.H"
#include "BRMeshRefine.H"
#include "BiCGStabSolver.H"
#include "CH_HDF5.H"
#include "DebugDump.H"
#include "DiagnosticVariables.hpp"
#include "FABView.H"
#include "FArrayBox.H"
#include "FourthOrderCFInterp.H"
#include "GRChomboVariables.hpp"
#include "Grids.hpp"
#include "LevelData.H"
#include "LoadBalance.H"
#include "MultigridVariables.hpp"
#include "MultilevelLinearOp.H"
#include "ParmParse.H"
#include "SimulationParameters.hpp"
#include <iostream>

#ifdef CH_Linux
// Should be undefined by default
// #define TRAP_FPE
#undef TRAP_FPE
#endif

#ifdef TRAP_FPE
static void enableFpExceptions();
#endif

using std::cerr;

/// This function outputs the variable data during the solver updates
/// which helps to check convergence and debug
template <typename method_t, typename matter_t>
void output_solver_data(
    const Vector<LevelData<FArrayBox> *> &a_constraint_vars,
    const Vector<LevelData<FArrayBox> *> &a_multigrid_vars,
    const Vector<LevelData<FArrayBox> *> &a_diagnostic_vars,
    const Vector<DisjointBoxLayout> &a_grids,
    const SimulationParameters<method_t, matter_t> &a_params, const int iter)
{
#ifdef CH_USE_HDF5

#if CH_SPACEDIM == 2
    string filename(a_params.base_params.output_path);
    filename += "NL_iteration";
#else
    string filename(a_params.base_params.output_path);
    filename += "NL_iteration";
#endif

    char suffix[30];
    sprintf(suffix, "_0%04i.hdf5", iter);
    filename += suffix;

    int nConComp = a_constraint_vars[0]->nComp();
    int nMGComp = a_multigrid_vars[0]->nComp();
    int nDiagComp = a_diagnostic_vars[0]->nComp();
    int totalComp = nConComp + nMGComp + nDiagComp;
    CH_assert(totalComp ==
              NUM_CONSTRAINT_VARS + NUM_MULTIGRID_VARS + NUM_DIAGNOSTIC_VARS);
    // set names
    Vector<string> variable_names(totalComp);

    for (int icomp = 0; icomp < NUM_CONSTRAINT_VARS; icomp++)
    {
        variable_names[icomp] =
            ConstraintVariables::variable_names[icomp];
    }
    for (int icomp = 0; icomp < NUM_METRIC_VARS; icomp++)
    {
        variable_names[NUM_CONSTRAINT_VARS + icomp] =
            MetricVariables::variable_names[icomp];
    }
    for (int icomp = 0; icomp < NUM_MULTIGRID_VARS - NUM_METRIC_VARS; icomp++)
    {
        variable_names[NUM_CONSTRAINT_VARS + NUM_METRIC_VARS + icomp] =
            MatterVariables::variable_names[icomp];
    }
    for (int icomp = 0; icomp < NUM_DIAGNOSTIC_VARS; icomp++)
    {
        variable_names[NUM_CONSTRAINT_VARS + NUM_MULTIGRID_VARS + icomp] =
            DiagnosticVariables::variable_names[icomp];
    }

    // check the domain sizes are the same
    CH_assert(a_constraint_vars.size() == a_multigrid_vars.size());
    CH_assert(a_constraint_vars.size() == a_diagnostic_vars.size());
    // Solver only has second order stencils so 1 ghost
    IntVect ghosts = 1 * IntVect::Unit;
    Vector<LevelData<FArrayBox> *> tempData(a_constraint_vars.size(), NULL);
    for (int level = 0; level < a_constraint_vars.size(); level++)
    {
        tempData[level] =
            new LevelData<FArrayBox>(a_grids[level], totalComp, ghosts);
        Interval conComps(0, NUM_CONSTRAINT_VARS - 1);
        Interval mgComps(NUM_CONSTRAINT_VARS,
                         NUM_CONSTRAINT_VARS + NUM_MULTIGRID_VARS - 1);
        Interval diagComps(NUM_CONSTRAINT_VARS + NUM_MULTIGRID_VARS,
                           NUM_CONSTRAINT_VARS + NUM_MULTIGRID_VARS +
                               NUM_DIAGNOSTIC_VARS - 1);
        ProblemDomain grown_domain_box = a_grids[level].physDomain();
        grown_domain_box.grow(ghosts);
        Copier boundary_copier;
        boundary_copier.ghostDefine(a_constraint_vars[level]->getBoxes(),
                                    a_constraint_vars[level]->getBoxes(), grown_domain_box,
                                    ghosts, ghosts);
        a_constraint_vars[level]->copyTo(a_constraint_vars[level]->interval(), *tempData[level],
                              conComps, boundary_copier);
        a_multigrid_vars[level]->copyTo(a_multigrid_vars[level]->interval(),
                                        *tempData[level], mgComps,
                                        boundary_copier);
        a_diagnostic_vars[level]->copyTo(a_diagnostic_vars[level]->interval(),
                                         *tempData[level], diagComps);
    }
    Real fakeTime = iter * 1.0;
    Real fakeDt = 1.0;
    WriteAMRHierarchyHDF5(filename, a_grids, tempData, variable_names,
                          a_params.grid_params.coarsestDomain.domainBox(),
                          a_params.grid_params.coarsestDx, fakeDt, fakeTime,
                          a_params.grid_params.refRatio,
                          a_params.grid_params.numLevels);
    // clean up temporary storage
    for (int level = 0; level < a_constraint_vars.size(); level++)
    {
        delete tempData[level];
        tempData[level] = NULL;
    }
#endif
}

/// This function outputs the final data after the solver is finished
/// in a format that GRChombo can read in for an AMRTimeDependent run

template <typename method_t, typename matter_t>
void output_final_data(const Vector<LevelData<FArrayBox> *> &a_multigrid_vars,
                       const Vector<DisjointBoxLayout> &a_grids,
                       const Vector<RealVect> &a_vectDx,
                       const Vector<ProblemDomain> &a_vectDomains,
                       const SimulationParameters<method_t, matter_t> &a_params,
                       const string filename = "InitialConditionsFinal.3d.hdf5")
{
#ifdef CH_USE_HDF5

    // set up full data and write to hdf5 file ready for GRChombo
    HDF5Handle handle(filename, HDF5Handle::CREATE);

    HDF5HeaderData header;
    header.m_int["max_level"] = a_params.grid_params.maxLevel;
    header.m_int["num_levels"] = a_params.grid_params.maxLevel + 1;
    header.m_int["iteration"] = 0;
    header.m_real["time"] = 0.0;

    for (int level = 0; level < a_multigrid_vars.size(); ++level)
    {
        char headername1[100];
        sprintf(headername1, "regrid_interval_%d", level);
        header.m_int[headername1] = 1;

        char headername2[100];
        sprintf(headername2, "steps_since_regrid_%d", level);
        // setting this to 1 forces a regrid on loading, which can be useful
        header.m_int[headername2] = 1;
    }

    header.m_int["num_components"] = NUM_GRCHOMBO_VARS;

    char comp_str[30];
    for (int comp = 0; comp < NUM_GRCHOMBO_VARS; ++comp)
    {
        sprintf(comp_str, "component_%d", comp);
        header.m_string[comp_str] = GRChomboVariables::variable_names[comp];
    }

    // write header data to file
    header.writeToFile(handle);

    // want output data to have 3 ghost cells to match GRChombo
    int num_ghosts = 3;
    IntVect ghost_vector = num_ghosts * IntVect::Unit;

    // loop over levels to prepare the output data
    Vector<LevelData<FArrayBox> *> grchombo_vars(a_multigrid_vars.size(), NULL);
    for (int level = 0; level < a_multigrid_vars.size(); level++)
    {
        grchombo_vars[level] = new LevelData<FArrayBox>(
            a_grids[level], NUM_GRCHOMBO_VARS, ghost_vector);

        // Set the values of the grchombo vars from the multigrid data
        // within the domain
        set_output_data(*grchombo_vars[level], *a_multigrid_vars[level],
                        a_params, a_vectDx[level]);

        // fill the boundary cells in all directions (may have more ghosts
        // than in solver so will need to fill them appropriately)
        BoundaryConditions solver_boundaries;
        solver_boundaries.define(a_vectDx[level][0],
                                 a_params.grid_params.boundary_params,
                                 a_vectDomains[level], num_ghosts);

        // this will populate the grchombo boundaries according to the BCs
        solver_boundaries.fill_grchombo_boundaries(Side::Lo,
                                                   *grchombo_vars[level]);
        solver_boundaries.fill_grchombo_boundaries(Side::Hi,
                                                   *grchombo_vars[level]);

        DisjointBoxLayout grown_grids;
        if (a_params.grid_params.boundary_params.nonperiodic_boundaries_exist)
        {
            solver_boundaries.expand_grids_to_boundaries(grown_grids,
                                                         a_grids[level]);
        }
        else
        { // nothing to do if periodic BCs
            grown_grids = a_grids[level];
        }
        Copier exchange_copier;
        exchange_copier.exchangeDefine(grown_grids, ghost_vector);
        grchombo_vars[level]->exchange(grchombo_vars[level]->interval(),
                                       exchange_copier);

        // fill interlevel ghosts from coarser level
        if (level > 0)
        {
            FourthOrderCFInterp m_patcher;
            m_patcher.define(a_grids[level], a_grids[level - 1],
                             NUM_GRCHOMBO_VARS, a_vectDomains[level - 1],
                             a_params.grid_params.refRatio[level], num_ghosts);
            m_patcher.coarseFineInterp(*grchombo_vars[level],
                                       *grchombo_vars[level - 1], 0, 0,
                                       NUM_GRCHOMBO_VARS);
        }
    }

    // Now write out the data
    for (int level = 0; level < a_multigrid_vars.size(); level++)
    {
        HDF5HeaderData level_header;
        // Now continue the write out
        char level_str[20];
        sprintf(level_str, "%d", level);
        const std::string label = std::string("level_") + level_str;
        handle.setGroup(label);

        // assign arbitrary values to the params
        // most can be overwritten in GRChombo setup
        level_header.m_int["ref_ratio"] = a_params.grid_params.refRatio[level];
        level_header.m_int["tag_buffer_size"] = 3;
        level_header.m_real["dx"] = a_vectDx[level][0];
        level_header.m_real["dt"] = 0.25 * a_vectDx[level][0];
        level_header.m_real["time"] = 0.0;
        level_header.m_box["prob_domain"] = a_vectDomains[level].domainBox();

        // Setup the periodicity info
        for (int dir = 0; dir < SpaceDim; ++dir)
        {
            bool periodic = a_params.grid_params.periodic[dir];
            char dir_str[20];
            sprintf(dir_str, "%d", dir);
            const std::string periodic_label =
                std::string("is_periodic_") + dir_str;
            level_header.m_int[periodic_label] = periodic;
        }

        level_header.writeToFile(handle);
        write(handle, a_grids[level]);
        write(handle, *grchombo_vars[level], "data", ghost_vector);
    }

    // shut the file
    handle.close();

    // clean up temporary storage
    for (int level = 0; level < a_multigrid_vars.size(); level++)
    {
        delete grchombo_vars[level];
        grchombo_vars[level] = NULL;
    }

#endif
}

template <typename method_t, typename matter_t>
void set_output_data(LevelData<FArrayBox> &a_grchombo_vars,
                     LevelData<FArrayBox> &a_multigrid_vars,
                     const SimulationParameters<method_t, matter_t> &a_params,
                     const RealVect &a_dx)
{

    CH_assert(a_grchombo_vars.nComp() == NUM_GRCHOMBO_VARS);
    CH_assert(a_multigrid_vars.nComp() == NUM_MULTIGRID_VARS);

    PsiAndAijFunctions psi_and_Aij_functions(
        a_params.psi_and_Aij_params);

    DataIterator dit = a_grchombo_vars.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &grchombo_vars_box = a_grchombo_vars[dit()];
        FArrayBox &multigrid_vars_box = a_multigrid_vars[dit()];

        // first set everything to zero
        for (int comp = 0; comp < NUM_GRCHOMBO_VARS; comp++)
        {
            grchombo_vars_box.setVal(0.0, comp);
        }

        // now set non zero terms - const across whole box
        // Conformally flat, and lapse = 1
        grchombo_vars_box.setVal(1.0, c_h11);
        grchombo_vars_box.setVal(1.0, c_h22);
        grchombo_vars_box.setVal(1.0, c_h33);
        grchombo_vars_box.setVal(1.0, c_lapse);

        // now non constant terms by location
        Box solver_ghosted_box = multigrid_vars_box.box();
        BoxIterator bit(solver_ghosted_box);

        for (bit.begin(); bit.ok(); ++bit)
        {
            // work out location on the grid
            IntVect iv = bit();
            RealVect loc;
            Grids::get_loc(loc, iv, a_dx, a_params.grid_params.center);

            // GRChombo conformal factor chi = psi^-4
            Real psi_bh = psi_and_Aij_functions.compute_bowenyork_psi(loc);
            Real chi = pow(multigrid_vars_box(iv, c_psi_reg) + psi_bh, -4.0);
            grchombo_vars_box(iv, c_chi) = chi;
            Real factor = pow(chi, 1.5);

            // Copy phi and Aij across - note this is now \tilde Aij not
            // \bar Aij
#ifdef USE_COMPLEX_SCALAR_MATTER
            grchombo_vars_box(iv, c_phi) = multigrid_vars_box(iv, c_phi_re);
            grchombo_vars_box(iv, c_Pi) = multigrid_vars_box(iv, c_Pi_re);
            grchombo_vars_box(iv, c_phi2) = multigrid_vars_box(iv, c_phi_im);
            grchombo_vars_box(iv, c_Pi2) = multigrid_vars_box(iv, c_Pi_im);
#else
            grchombo_vars_box(iv, c_phi) = multigrid_vars_box(iv, c_phi_0);
            grchombo_vars_box(iv, c_Pi) = multigrid_vars_box(iv, c_Pi_0);
#endif

            grchombo_vars_box(iv, c_K) = multigrid_vars_box(iv, c_K_0);
            grchombo_vars_box(iv, c_A11) =
                multigrid_vars_box(iv, c_A11_0) * factor;
            grchombo_vars_box(iv, c_A12) =
                multigrid_vars_box(iv, c_A12_0) * factor;
            grchombo_vars_box(iv, c_A13) =
                multigrid_vars_box(iv, c_A13_0) * factor;
            grchombo_vars_box(iv, c_A22) =
                multigrid_vars_box(iv, c_A22_0) * factor;
            grchombo_vars_box(iv, c_A23) =
                multigrid_vars_box(iv, c_A23_0) * factor;
            grchombo_vars_box(iv, c_A33) =
                multigrid_vars_box(iv, c_A33_0) * factor;
        }
    }
}

#endif
