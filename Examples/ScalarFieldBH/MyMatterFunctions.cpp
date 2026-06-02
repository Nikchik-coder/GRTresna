/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#include "ScalarField.hpp"
#include <array>
#include <cmath>

namespace
{
// Smooth (Cartesian) azimuthal modulation of the lump. m = 0 is axisymmetric
// (no angular momentum possible); m >= 1 breaks axisymmetry so that a rigidly
// rotating pattern carries net L_z. Using Cartesian harmonics avoids the
// coordinate singularity of cos(m*phi) on the z-axis.
inline Real angular_factor(int m, Real dx, Real dy, Real w)
{
    switch (m)
    {
    case 1:
        return dx / w; // dipole (m = 1)
    case 2:
        return (dx * dx - dy * dy) / (w * w); // quadrupole (m = 2)
    default:
        return 1.0; // axisymmetric
    }
}

// The momentum-carrying lump contribution to phi (without the constant
// offset). Returns 0 when the lump is disabled (amplitude == 0).
inline Real lump_phi_value(const RealVect &loc,
                           const MatterParams::params_t &p)
{
    if (p.lump_amp == 0.0)
    {
        return 0.0;
    }
    const Real dx = loc[0] - p.lump_center[0];
    const Real dy = loc[1] - p.lump_center[1];
    const Real dz = loc[2] - p.lump_center[2];
    const Real r2 = dx * dx + dy * dy + dz * dz;
    const Real w  = p.lump_width;
    const Real env = exp(-r2 / (2.0 * w * w));
    return p.lump_amp * angular_factor(p.lump_mode, dx, dy, w) * env;
}
} // namespace

Real ScalarField::my_potential_function(const Real &phi_here) const
{
    return 0.5 * pow(m_matter_params.scalar_mass * phi_here, 2.0);
}

Real ScalarField::my_phi_function(const RealVect &loc) const
{
    Real rr = sqrt(loc[0] * loc[0] + loc[1] * loc[1] + loc[2] * loc[2]);
    Real phi_legacy = m_matter_params.phi_0 +
                      m_matter_params.dphi *
                          exp(-rr / m_matter_params.dphi_length);
    return phi_legacy + lump_phi_value(loc, m_matter_params);
}

Real ScalarField::my_Pi_function(const RealVect &loc) const
{
    Real rr = sqrt(loc[0] * loc[0] + loc[1] * loc[1] + loc[2] * loc[2]);
    Real pi_legacy = m_matter_params.pi_0 +
                     m_matter_params.dpi * exp(-rr / m_matter_params.dpi_length);

    const auto &p = m_matter_params;

    const bool has_kinematics =
        (p.lump_velocity[0] != 0.0 || p.lump_velocity[1] != 0.0 ||
         p.lump_velocity[2] != 0.0 || p.lump_omega != 0.0);

    if (p.lump_amp == 0.0 || !has_kinematics)
    {
        return pi_legacy;
    }

    // Pi = (d/dt) phi for a pattern that is boosted with velocity v and
    // rigidly rotated about z with rate omega:
    //   Pi = -( v . grad phi ) - omega * d phi / d(azimuth).
    // This makes the momentum density S_i = -Pi d_i phi carry net linear
    // momentum P_i ~ v_i and net angular momentum L_z ~ omega (the latter
    // requires lump_mode >= 1 so that d phi / d(azimuth) != 0).
    const Real eps = 1.0e-3 * p.lump_width;
    auto lp = [&](Real x, Real y, Real z)
    {
        RealVect l;
        l[0] = x;
        l[1] = y;
        l[2] = z;
        return lump_phi_value(l, p);
    };

    const Real gx =
        (lp(loc[0] + eps, loc[1], loc[2]) - lp(loc[0] - eps, loc[1], loc[2])) /
        (2.0 * eps);
    const Real gy =
        (lp(loc[0], loc[1] + eps, loc[2]) - lp(loc[0], loc[1] - eps, loc[2])) /
        (2.0 * eps);
    const Real gz =
        (lp(loc[0], loc[1], loc[2] + eps) - lp(loc[0], loc[1], loc[2] - eps)) /
        (2.0 * eps);

    // Linear-momentum (boost) part: -v . grad phi
    const Real boost = -(p.lump_velocity[0] * gx + p.lump_velocity[1] * gy +
                         p.lump_velocity[2] * gz);

    // Angular-momentum (rigid rotation about z through lump centre) part:
    // -omega * d phi / d(azimuth) = -omega * (lx d_y phi - ly d_x phi)
    const Real lx  = loc[0] - p.lump_center[0];
    const Real ly  = loc[1] - p.lump_center[1];
    const Real rot = -p.lump_omega * (lx * gy - ly * gx);

    return pi_legacy + boost + rot;
}
