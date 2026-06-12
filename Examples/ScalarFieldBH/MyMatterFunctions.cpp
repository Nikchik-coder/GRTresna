/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#include "ScalarField.hpp"
#include <array>
#include <cmath>

// The per-lump phi / Pi profile helpers live in MatterParams.hpp so the matter
// painter (here) and the constraint-solve energy/momentum assembly
// (ScalarField::compute_emtensor) share exactly the same analytic cloud.

Real ScalarField::my_potential_function(const Real &phi_here) const
{
    const Real mphi = m_matter_params.scalar_mass * phi_here;
    const Real phi4 = pow(phi_here, 4.0);
    return 0.5 * mphi * mphi -
           0.25 * m_matter_params.scalar_lambda * phi4;
}

Real ScalarField::my_phi_function(const RealVect &loc) const
{
    Real rr = sqrt(loc[0] * loc[0] + loc[1] * loc[1] + loc[2] * loc[2]);
    Real phi = m_matter_params.phi_0 +
               m_matter_params.dphi * exp(-rr / m_matter_params.dphi_length);
    for (const auto &L : m_matter_params.lumps)
    {
        phi += MatterParams::lump_phi(loc, L);
    }
    return phi;
}

Real ScalarField::my_Pi_function(const RealVect &loc) const
{
    Real rr = sqrt(loc[0] * loc[0] + loc[1] * loc[1] + loc[2] * loc[2]);
    Real Pi = m_matter_params.pi_0 +
              m_matter_params.dpi * exp(-rr / m_matter_params.dpi_length);
    for (const auto &L : m_matter_params.lumps)
    {
        Pi += MatterParams::lump_pi(loc, L);
    }
    return Pi;
}
