// Only compiled for examples that select the complex (boson-star) matter
// model; otherwise this translation unit is intentionally empty so it does not
// clash with the canonical ScalarField in shared builds.
#ifdef USE_COMPLEX_SCALAR_MATTER

#include "ComplexScalarField.hpp"
#include "DerivativeOperators.hpp"
#include "Grids.hpp"
#include "MultigridVariables.hpp"
#include <cmath>

namespace
{
Real potential_value(Real mod2, Real mass, Real lam)
{
    const Real mphi = mass * std::sqrt(mod2);
    return 0.5 * mphi * mphi - 0.25 * lam * mod2 * mod2;
}
} // namespace

void ComplexScalarField::initialise_matter_vars(
    LevelData<FArrayBox> &a_multigrid_vars, const RealVect &a_dx) const
{
    CH_assert(a_multigrid_vars.nComp() == NUM_MULTIGRID_VARS);

    const Real omega =
        (m_params.omega > 0.0) ? m_params.omega : m_params.scalar_mass;

    DataIterator dit = a_multigrid_vars.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &box = a_multigrid_vars[dit()];
        BoxIterator bit(box.box());
        for (bit.begin(); bit.ok(); ++bit)
        {
            IntVect iv = bit();
            RealVect loc;
            Grids::get_loc(loc, iv, a_dx, center);

            const Real r = std::sqrt(loc[0] * loc[0] + loc[1] * loc[1] +
                                     loc[2] * loc[2]);
            const Real phi1 =
                BosonStarParams::phi0_profile(r, m_params.phi_c,
                                                m_params.profile_width);

            box(iv, c_phi_re) = phi1;
            box(iv, c_phi_im) = 0.0;
            box(iv, c_Pi_re)  = 0.0;
            // alpha = 1 during initial paint; post-solve correction in wrapper.
            box(iv, c_Pi_im)  = -omega * phi1;
        }
    }
}

emtensor_t ComplexScalarField::compute_emtensor(
    const IntVect a_iv, const RealVect &a_dx,
    FArrayBox &a_multigrid_vars_box) const
{
    emtensor_t out;
    RealVect loc;
    Grids::get_loc(loc, a_iv, a_dx, center);

    Real psi_reg = a_multigrid_vars_box(a_iv, c_psi_reg);
    Real psi_bh  = psi_and_Aij_functions->compute_bowenyork_psi(loc);
    Real psi_0   = psi_reg + psi_bh;
    const Real chi = std::pow(psi_0, -4.0);

    const Real omega =
        (m_params.omega > 0.0) ? m_params.omega : m_params.scalar_mass;
    const Real r = std::sqrt(loc[0] * loc[0] + loc[1] * loc[1] +
                             loc[2] * loc[2]);

    const Real phi1 =
        BosonStarParams::phi0_profile(r, m_params.phi_c, m_params.profile_width);
    const Real phi2 = 0.0;
    const Real pi1  = 0.0;
    const Real pi2  = -omega * phi1;

    Real dphi1[SpaceDim];
    if (r > 1.0e-12)
    {
        const Real dphidr =
            BosonStarParams::dphi0_dr(r, m_params.phi_c, m_params.profile_width);
        for (int i = 0; i < SpaceDim; ++i)
            dphi1[i] = dphidr * loc[i] / r;
    }
    else
    {
        for (int i = 0; i < SpaceDim; ++i)
            dphi1[i] = 0.0;
    }
    Real dphi2[SpaceDim] = {0.0, 0.0, 0.0};

    Real grad1_sq = 0.0;
    Real grad2_sq = 0.0;
    for (int i = 0; i < SpaceDim; ++i)
    {
        grad1_sq += dphi1[i] * dphi1[i];
        grad2_sq += dphi2[i] * dphi2[i];
    }

    const Real mod2 = phi1 * phi1 + phi2 * phi2;
    const Real V    = potential_value(mod2, m_params.scalar_mass,
                                      m_params.scalar_lambda);

    // Phantom sign: flip entire T_ab for sign == -1 while preserving U(1).
    const Real sign = m_params.sign;
    out.rho = sign * (0.5 * pi1 * pi1 + 0.5 * chi * grad1_sq +
                      0.5 * pi2 * pi2 + 0.5 * chi * grad2_sq + V);

    for (int i = 0; i < SpaceDim; ++i)
        out.Si[i] = sign * (-pi1 * dphi1[i] - pi2 * dphi2[i]);

    return out;
}

#endif /* USE_COMPLEX_SCALAR_MATTER */
