/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#ifndef MATTERPARAMS_HPP_
#define MATTERPARAMS_HPP_

#include "GRParmParse.hpp"
#include "REAL.H"
#include <array>

namespace MatterParams
{

struct params_t
{
    // Legacy spherically-symmetric scalar profile (carries no net momentum)
    Real phi_0;
    Real dphi;
    Real dphi_length;
    Real pi_0;
    Real dpi;
    Real dpi_length;
    Real scalar_mass;

    // Momentum-carrying scalar "cloud" (a localised lump). All default to
    // zero/off so existing parameter files reproduce the legacy behaviour.
    // The lump is added on top of the legacy profile and its conjugate
    // momentum Pi is constructed so that the configuration carries net linear
    // and/or angular momentum (S_i = -Pi d_i phi is then non-zero), which the
    // momentum constraint solver responds to.
    Real lump_amp;                      // amplitude (0 => lump disabled)
    Real lump_width;                    // Gaussian width of the lump
    std::array<Real, 3> lump_center;    // lump centre, relative to grid centre
    std::array<Real, 3> lump_velocity;  // boost velocity v => linear momentum
    Real lump_omega;                    // rigid rotation rate about z => L_z
    int lump_mode;                      // azimuthal modulation m (>=1 for L_z)
};

inline void read_params(GRParmParse &pp, params_t &matter_params)
{
    pp.get("phi_0", matter_params.phi_0);
    pp.get("dphi", matter_params.dphi);
    pp.get("dphi_length", matter_params.dphi_length);
    pp.get("pi_0", matter_params.pi_0);
    pp.get("dpi", matter_params.dpi);
    pp.get("dpi_length", matter_params.dpi_length);
    pp.get("scalar_mass", matter_params.scalar_mass);

    // Momentum-carrying lump (optional, defaults reproduce legacy behaviour)
    pp.load("lump_amp", matter_params.lump_amp, 0.0);
    pp.load("lump_width", matter_params.lump_width, 1.0);
    pp.load("lump_omega", matter_params.lump_omega, 0.0);
    pp.load("lump_mode", matter_params.lump_mode, 0);
    if (pp.contains("lump_center"))
        pp.load("lump_center", matter_params.lump_center);
    else
        matter_params.lump_center = {0.0, 0.0, 0.0};
    if (pp.contains("lump_velocity"))
        pp.load("lump_velocity", matter_params.lump_velocity);
    else
        matter_params.lump_velocity = {0.0, 0.0, 0.0};
}

}; // namespace MatterParams

#endif
