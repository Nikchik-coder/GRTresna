/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#ifndef GRIDS_HPP_
#define GRIDS_HPP_

#include "BoundaryConditions.hpp"
#include "CoarseAverage.H"
#include "FilesystemTools.hpp"
#include "GRParmParse.hpp"
#include "IntVect.H"
#include "IntVectSet.H"
#include "MultilevelLinearOp.H"
#include "ProblemDomain.H"
#include "REAL.H"
#include "RealVect.H"
#include "TaggingCriterion.hpp"
#include "computeNorm.H"
#include "computeSum.H"

class Grids
{
  public:
    struct params_t
    {
        IntVect nCells; //
        int maxGridSize;
        int blockFactor;
        int bufferSize;
        Real fillRatio;
        Real refineThresh;
        Real regrid_radius;
        int coefficient_average_type;

        Vector<int> periodic;
        bool periodic_directions_exist;
        int domBcType;
        int maxLevel;  //
        int numLevels; //
        Vector<int> refRatio;
        ProblemDomain coarsestDomain;
        Real coarsestDx;
        std::array<double, SpaceDim> center; // grid center
        BoundaryConditions::params_t
            boundary_params; // set boundaries in each dir
        RealVect domainLength;
        RealVect probLo;
        RealVect probHi;

        int num_ghosts;
    };

    Grids(params_t a_grid_params, TaggingCriterion *a_tagging_criterion,
          bool a_readin_matter_data)
        : m_grid_params(a_grid_params), tagging_criterion(a_tagging_criterion),
          readin_matter_data(a_readin_matter_data){};

    // get location:
    // This takes an IntVect and writes the physical coordinates to a RealVect
    static void get_loc(RealVect &a_out_loc, const IntVect &a_iv,
                        const RealVect &a_dx,
                        const std::array<double, SpaceDim> center);

    void read_grids(std::string input_filename);

    void set_grids();

    static void read_params(GRParmParse &pp, params_t &grid_params);

    void define_operator(MultilevelLinearOp<FArrayBox> &mlOp,
                         Vector<RefCountedPtr<LevelData<FArrayBox>>> &aCoef,
                         Vector<RefCountedPtr<LevelData<FArrayBox>>> &bCoef,
                         const Real &a_alpha, const Real &a_beta);

    void
    fill_ghosts_correct_coarse(Vector<LevelData<FArrayBox> *> multigrid_vars,
                               bool filling_solver_vars);

    void update_psi0(Vector<LevelData<FArrayBox> *> multigrid_vars,
                     Vector<LevelData<FArrayBox> *> constraint_vars,
                     bool deactivate_zero_mode, Real relaxation = 1.0,
                     Real psi_floor = -1.0);

    params_t m_grid_params;

    TaggingCriterion *tagging_criterion;

    Vector<DisjointBoxLayout> grids_data;

    Vector<RealVect> vectDx;

    Vector<ProblemDomain> vectDomain;

    Real compute_sum(Vector<LevelData<FArrayBox> *> vars,
                     const Interval &a_interval)
    {
        return computeSum(vars, m_grid_params.refRatio,
                          m_grid_params.coarsestDx, a_interval);
    }

    Real compute_norm(Vector<LevelData<FArrayBox> *> vars,
                      const Interval &a_interval)
    {
        return computeNorm(vars, m_grid_params.refRatio,
                           m_grid_params.coarsestDx, a_interval);
    }

    Real compute_max(Vector<LevelData<FArrayBox> *> vars,
                     const Interval &a_interval)
    {
        return computeMax(vars, m_grid_params.refRatio, a_interval);
    }

  private:
    bool readin_matter_data;

    void set_domains_and_dx(Vector<ProblemDomain> &vectDomain,
                            Vector<RealVect> &vectDx);

    void set_tag_cells(Vector<LevelData<FArrayBox> *> &vect_tagging_criterion,
                       Vector<IntVectSet> &tagVect, Vector<RealVect> &vectDx,
                       Vector<ProblemDomain> &vectDomain, const int tags_grow,
                       const int baseLevel, int numLevels_tag);
};

#endif /* GRIDS_HPP_ */