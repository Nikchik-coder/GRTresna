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

void paint_boson_fields(const BosonStarParams::params_t &params,
                        const RealVect &loc, Real &phi1, Real &phi2,
                        Real &pi1, Real &pi2)
{
    phi1 = BosonStarParams::total_phi1(loc, params);
    phi2 = 0.0;
    pi1  = 0.0;
    const Real omega =
        (params.omega > 0.0) ? params.omega : params.scalar_mass;
    // alpha = 1 during initial paint; post-solve correction in wrapper.
    pi2 = -omega * phi1;
}
} // namespace

void ComplexScalarField::initialise_matter_vars(
    LevelData<FArrayBox> &a_multigrid_vars, const RealVect &a_dx) const
{
    CH_assert(a_multigrid_vars.nComp() == NUM_MULTIGRID_VARS);

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

            Real phi1, phi2, pi1, pi2;
            paint_boson_fields(m_params, loc, phi1, phi2, pi1, pi2);

            box(iv, c_phi_re) = phi1;
            box(iv, c_phi_im) = phi2;
            box(iv, c_Pi_re)  = pi1;
            box(iv, c_Pi_im)  = pi2;
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

    Real phi1, phi2, pi1, pi2;
    paint_boson_fields(m_params, loc, phi1, phi2, pi1, pi2);

    std::array<Real, 3> dphi1;
    BosonStarParams::total_grad_phi1(loc, m_params, dphi1);
    Real dphi2[3] = {0.0, 0.0, 0.0};

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

    const Real sign = m_params.sign;
    out.rho = sign * (0.5 * pi1 * pi1 + 0.5 * chi * grad1_sq +
                      0.5 * pi2 * pi2 + 0.5 * chi * grad2_sq + V);

    for (int i = 0; i < SpaceDim; ++i)
        out.Si[i] = sign * (-pi1 * dphi1[i] - pi2 * dphi2[i]);

    return out;
}

#endif /* USE_COMPLEX_SCALAR_MATTER */
