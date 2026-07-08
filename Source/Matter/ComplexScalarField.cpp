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
// V = 1/2 m^2 |Phi|^2 - 1/4 lambda |Phi|^4 + 1/6 mu |Phi|^6.  The sextic
// stabiliser (mu > 0) is REQUIRED for a genuine 3D Q-ball: a pure attractive
// quartic is the critical case and either collapses or disperses.  Must match
// GRTeclyn's ComplexScalarPotential.hpp so the constraint solve and the
// evolution share the same T_ab (otherwise the solved metric backs a different
// soliton than the one that is evolved, and the lump relaxes/disperses at t=0).
Real potential_value(Real mod2, Real mass, Real lam, Real mu)
{
    const Real mphi = mass * std::sqrt(mod2);
    return 0.5 * mphi * mphi - 0.25 * lam * mod2 * mod2 +
           (1.0 / 6.0) * mu * mod2 * mod2 * mod2;
}

void paint_boson_fields(const BosonStarParams::params_t &params,
                        const RealVect &loc, Real &phi1, Real &phi2,
                        Real &pi1, Real &pi2)
{
    phi1 = BosonStarParams::total_phi1(loc, params);
    phi2 = 0.0;
    pi1  = BosonStarParams::total_pi1(loc, params);
    const Real omega =
        (params.omega > 0.0) ? params.omega : params.scalar_mass;
    // Global U(1) phase velocity: Pi_im = -(omega/alpha(r)) phi1.  For a
    // self-gravitating star alpha(r) < 1 in the core, so this is the correct
    // stationary momentum; for a flat-space table alpha == 1 and it reduces to
    // the old -omega*phi1.
    pi2 = BosonStarParams::total_pi2(loc, params, omega);
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

            Real phi1 = 0.0, phi2 = 0.0, pi1 = 0.0, pi2 = 0.0;
            if (!m_params.lumps.empty())
            {
                // Per-lump painting so winding (rotating) lumps get a genuine
                // phi2 = f sin(m phi_az) and the correct two-channel U(1)
                // momentum; non-winding lumps keep the legacy boson-star ansatz.
                const Real omega = (m_params.omega > 0.0) ? m_params.omega
                                                          : m_params.scalar_mass;
                for (const auto &L : m_params.lumps)
                {
                    if (L.amp == 0.0)
                        continue;
                    const Real alpha_k = BosonStarParams::lump_alpha(loc, L);
                    if (L.winding != 0)
                    {
                        Real f1, f2;
                        BosonStarParams::lump_phi_winding(loc, L, f1, f2);
                        phi1 += f1;
                        phi2 += f2;
                        // Per-lump phase velocity (the wormhole rotation rate);
                        // fall back to the global omega only if unset.
                        const Real omega_k = (L.omega != 0.0) ? L.omega : omega;
                        // Pi = d_t Phi / alpha with Phi ~ e^{-i omega t}:
                        //   Pi1 = +(omega/alpha) phi2,  Pi2 = -(omega/alpha) phi1.
                        pi1 += (omega_k / alpha_k) * f2;
                        pi2 += -(omega_k / alpha_k) * f1;
                    }
                    else
                    {
                        const Real f = BosonStarParams::lump_phi1(loc, L);
                        phi1 += f;
                        pi1 += BosonStarParams::lump_pi1(loc, L);
                        pi2 += -(omega / alpha_k) * f;
                    }
                }
            }
            else
            {
                paint_boson_fields(m_params, loc, phi1, phi2, pi1, pi2);
            }

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

    // Independent-field matter model with PER-LUMP signs.  Each lump is a
    // separate complex scalar Phi_k = phi1_k + i phi2_k (phi2_k = 0 initially,
    // pi2_k = -omega phi1_k from the global U(1) phase velocity), and a lump
    // flagged EXOTIC contributes its whole stress-energy (kinetic + gradient +
    // potential + momentum) with a flipped sign -- a genuine ghost/phantom
    // source of negative energy.  Cross terms between lumps are dropped, exactly
    // as the real ScalarField does, so a config can mix normal and exotic lumps
    // (the warp/wormhole FTL geometry needs this).
    const Real omega =
        (m_params.omega > 0.0) ? m_params.omega : m_params.scalar_mass;

    Real rho_total = 0.0;
    Tensor<1, Real, SpaceDim> Si;
    for (int i = 0; i < SpaceDim; ++i)
        Si[i] = 0.0;

    if (!m_params.lumps.empty())
    {
        for (const auto &L : m_params.lumps)
        {
            if (L.amp == 0.0)
                continue;

            const Real sign = (L.exotic != 0) ? -1.0 : 1.0;
            const Real alpha_k = BosonStarParams::lump_alpha(loc, L);

            Real phi1_k, phi2_k, pi1_k, pi2_k;
            std::array<Real, 3> dphi1_k, dphi2_k;
            if (L.winding != 0)
            {
                // Genuine phase winding Phi = f e^{i(m phi_az)}: axisymmetric
                // |Phi|^2, and a real azimuthal momentum density that sources a
                // clean J_z via the momentum constraint (constraint-clean spin).
                BosonStarParams::lump_phi_winding(loc, L, phi1_k, phi2_k);
                const Real omega_k = (L.omega != 0.0) ? L.omega : omega;
                pi1_k = (omega_k / alpha_k) * phi2_k;
                pi2_k = -(omega_k / alpha_k) * phi1_k;
                BosonStarParams::lump_grad_phi_winding(loc, L, dphi1_k, dphi2_k);
            }
            else
            {
                phi1_k = BosonStarParams::lump_phi1(loc, L);
                phi2_k = 0.0;
                pi1_k  = BosonStarParams::lump_pi1(loc, L);
                // Stationary U(1) momentum uses the star's own lapse alpha(r):
                // Pi_im = -(omega/alpha) phi1 (alpha == 1 for a flat-space table).
                pi2_k = -(omega / alpha_k) * phi1_k;
                BosonStarParams::lump_grad_phi1(loc, L, dphi1_k);
                dphi2_k = {0.0, 0.0, 0.0}; // phi2_k = 0 => grad phi2_k = 0
            }

            Real grad1_sq = 0.0, grad2_sq = 0.0;
            for (int i = 0; i < SpaceDim; ++i)
            {
                grad1_sq += dphi1_k[i] * dphi1_k[i];
                grad2_sq += dphi2_k[i] * dphi2_k[i];
            }

            const Real mod2_k = phi1_k * phi1_k + phi2_k * phi2_k;
            const Real V_k    = potential_value(mod2_k, m_params.scalar_mass,
                                                m_params.scalar_lambda,
                                                m_params.scalar_mu);

            rho_total += sign * (0.5 * pi1_k * pi1_k + 0.5 * pi2_k * pi2_k +
                                 0.5 * chi * (grad1_sq + grad2_sq) + V_k);
            // Two-channel momentum density S_i = -(Pi1 d_i phi1 + Pi2 d_i phi2).
            for (int i = 0; i < SpaceDim; ++i)
                Si[i] += sign * (-(pi1_k * dphi1_k[i] + pi2_k * dphi2_k[i]));
        }
    }
    else
    {
        // Legacy single centered boson star (no lumps): use the global sign.
        Real phi1, phi2, pi1, pi2;
        paint_boson_fields(m_params, loc, phi1, phi2, pi1, pi2);

        std::array<Real, 3> dphi1;
        BosonStarParams::total_grad_phi1(loc, m_params, dphi1);

        Real grad1_sq = 0.0;
        for (int i = 0; i < SpaceDim; ++i)
            grad1_sq += dphi1[i] * dphi1[i];

        const Real mod2 = phi1 * phi1 + phi2 * phi2;
        const Real V    = potential_value(mod2, m_params.scalar_mass,
                                          m_params.scalar_lambda,
                                          m_params.scalar_mu);

        const Real sign = m_params.sign;
        rho_total = sign * (0.5 * pi1 * pi1 + 0.5 * chi * grad1_sq +
                            0.5 * pi2 * pi2 + V);
        for (int i = 0; i < SpaceDim; ++i)
            Si[i] = sign * (-pi1 * dphi1[i]);
    }

    out.rho = rho_total;
    for (int i = 0; i < SpaceDim; ++i)
        out.Si[i] = Si[i];

    return out;
}

#endif /* USE_COMPLEX_SCALAR_MATTER */
