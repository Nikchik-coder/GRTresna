/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#ifndef RHSTAGGING_HPP_
#define RHSTAGGING_HPP_

#include "DerivativeOperators.hpp"
#include "REAL.H"
#include "TaggingCriterion.hpp"
#include <cmath>

template <typename method_t, typename matter_t>
class RHSTagging : public TaggingCriterion
{
  public:
    RHSTagging(method_t *a_method, matter_t *a_matter, Real a_G_Newton)
        : TaggingCriterion(), method(a_method), matter(a_matter),
          G_Newton(a_G_Newton)
    {
    }

    ~RHSTagging() {}

    void set_regrid_condition(LevelData<FArrayBox> &a_condition,
                              LevelData<FArrayBox> &a_multigrid_vars,
                              const RealVect &a_dx,
                              const std::array<double, SpaceDim> center,
                              Real regrid_radius);

  private:
    method_t const *method;
    matter_t const *matter;

    const Real G_Newton;
};

template <typename method_t, typename matter_t>
void RHSTagging<method_t, matter_t>::set_regrid_condition(
    LevelData<FArrayBox> &a_condition, LevelData<FArrayBox> &a_multigrid_vars,
    const RealVect &a_dx, const std::array<double, SpaceDim> center,
    Real regrid_radius)
{
    DerivativeOperators derivs(a_dx);
    CH_assert(a_multigrid_vars.nComp() == NUM_MULTIGRID_VARS);

    method->initialise_method_vars(a_multigrid_vars, a_dx);
    matter->initialise_matter_vars(a_multigrid_vars, a_dx);

    DataIterator dit = a_condition.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &multigrid_vars_box = a_multigrid_vars[dit()];
        FArrayBox &condition_box = a_condition[dit()];
        condition_box.setVal(0.0, 0);

        Box unghosted_box = condition_box.box();
        BoxIterator bit(unghosted_box);
        for (bit.begin(); bit.ok(); ++bit)
        {
            // work out location on the grid
            IntVect iv = bit();
            RealVect loc;
            Grids::get_loc(loc, iv, a_dx, center);

            // Calculate the actual value of psi including BH part
            Real psi_reg = multigrid_vars_box(iv, c_psi_reg);
            Real psi_bh =
                method->psi_and_Aij_functions->compute_bowenyork_psi(loc);
            Real psi_0 = psi_reg + psi_bh;
            Real laplacian_psi_reg;
            derivs.scalar_Laplacian(laplacian_psi_reg, iv, multigrid_vars_box,
                                    c_psi_reg);

            // Get values of Aij
            Tensor<2, Real> Aij_reg;
            method->psi_and_Aij_functions->compute_ctt_Aij(
                Aij_reg, multigrid_vars_box, iv, a_dx, loc);
            Tensor<2, Real> Aij_bh;
            method->psi_and_Aij_functions->compute_bowenyork_Aij(Aij_bh, loc);
            // This is \bar  A_ij \bar A^ij
            Real A2_0 = 0.0;
            FOR2(i, j)
            {
                A2_0 += (Aij_reg[i][j] + Aij_bh[i][j]) *
                        (Aij_reg[i][j] + Aij_bh[i][j]);
            }

            // Compute emtensor components
            const auto emtensor =
                matter->compute_emtensor(iv, a_dx, multigrid_vars_box);

            if (regrid_radius > 0)
            {
                Real rr = sqrt(D_TERM(loc[0] * loc[0], +loc[1] * loc[1],
                                      +loc[2] * loc[2]));
                if (rr < regrid_radius)
                {
                    condition_box(iv, 0) = 1.0;
                }
            }
            else
            {
                // Tag on source magnitude. Exotic matter has rho < 0, so using
                // signed rho can cancel the very regions that need AMR.
                const Real safe_psi = (std::isfinite(psi_0) && psi_0 > 1.0e-12)
                                          ? psi_0
                                          : 1.0e-12;
                condition_box(iv, 0) =
                    2.0 * M_PI * G_Newton * abs(emtensor.rho) +
                    abs(0.125 * A2_0) + abs(log(safe_psi)) +
                    abs(laplacian_psi_reg) +
                    8.0 * M_PI * G_Newton *
                        (abs(emtensor.Si[0]) + abs(emtensor.Si[1]) +
                         abs(emtensor.Si[2]));
            }
        }
    }
}

#endif /* RHSTAGGING_HPP_*/