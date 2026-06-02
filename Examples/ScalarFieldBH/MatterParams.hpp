/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#ifndef MATTERPARAMS_HPP_
#define MATTERPARAMS_HPP_

#include "GRParmParse.hpp"
#include "REAL.H"
#include "RealVect.H"
#include <array>
#include <cmath>
#include <string>
#include <vector>

namespace MatterParams
{

// One momentum-carrying scalar "cloud". Its conjugate momentum is built so the
// cloud carries net linear momentum (velocity) and/or angular momentum (omega
// with mode >= 1); a superposition of lumps paints an arbitrary energy- and
// momentum-density distribution, which GRTresna turns into constraint-
// satisfying initial data.
struct lump_t
{
    Real amp                   = 0.0; // amplitude (0 => disabled)
    Real width                 = 5.0; // Gaussian width
    std::array<Real, 3> center = {0.0, 0.0, 0.0};   // centre rel. grid centre
    std::array<Real, 3> velocity = {0.0, 0.0, 0.0}; // boost => linear momentum
    Real omega                 = 0.0; // rigid rotation about z => L_z
    int mode                   = 0;   // azimuthal modulation m (>=1 for L_z)
    // Phantom/ghost flag. 0 => canonical (positive-energy) scalar lump; != 0 =>
    // EXOTIC: this lump's kinetic energy density and momentum density enter the
    // constraint solve with a flipped sign, i.e. it is a ghost field that
    // sources genuinely negative energy (NEC violation). Treated as an
    // independent field, so a config can mix normal and exotic lumps.
    int exotic                 = 0;
};

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

    // Momentum-carrying scalar basis. Empty => pure legacy spherical data.
    std::vector<lump_t> lumps;
};

inline void read_lump(GRParmParse &pp, const std::string &prefix, lump_t &L)
{
    pp.load((prefix + "amp").c_str(), L.amp, 0.0);
    pp.load((prefix + "width").c_str(), L.width, 5.0);
    pp.load((prefix + "omega").c_str(), L.omega, 0.0);
    pp.load((prefix + "mode").c_str(), L.mode, 0);
    pp.load((prefix + "exotic").c_str(), L.exotic, 0);
    if (pp.contains((prefix + "center").c_str()))
        pp.load((prefix + "center").c_str(), L.center);
    else
        L.center = {0.0, 0.0, 0.0};
    if (pp.contains((prefix + "velocity").c_str()))
        pp.load((prefix + "velocity").c_str(), L.velocity);
    else
        L.velocity = {0.0, 0.0, 0.0};
}

// --- Shared lump field helpers ------------------------------------------------
//
// These describe a single scalar "cloud" analytically and are used in two
// places: when painting the initial (phi, Pi) onto the grid (MyMatterFunctions)
// and when assembling the per-lump (optionally exotic) energy/momentum density
// for the constraint solve (ScalarField::compute_emtensor). Keeping them here,
// inline, guarantees both call sites use exactly the same profile.

// Effective-amplitude scale applied to EXOTIC (phantom) lumps only.
//
// HISTORICAL: this used to be 0.08 to damp exotic lumps into the only regime
// the old CTTKHybrid K = sign*sqrt(24 pi G rho) ansatz could solve -- that
// ansatz takes the square root of a NEGATIVE number for rho < 0 and produced
// NaN for any non-trivial exotic strength. The solver now uses a maximal-slicing
// (K = 0) York/Lichnerowicz path (params: maximal_slicing) that sources the
// matter energy elliptically and handles rho < 0 directly, with under-relaxation
// + a psi-positivity floor for the indefinite operators that strong exotic
// matter produces.
//
// With the new solver the limit on exotic strength is now PHYSICAL, not a code
// bug: a single exotic lump converges to a clean Hamiltonian residual at
// effective amplitude 0.05 (~3%) and 0.10 (~7%), but at 0.15 the residual jumps
// to ~98% -- the Lichnerowicz/York existence boundary, beyond which no
// asymptotically-flat constraint-satisfying data exists for this configuration.
// The search amplitude lives in [0, 0.35]; this scale maps an exotic lump's
// effective amplitude into [0, ~0.09] so every exotic candidate stays inside
// the convergent regime (<~6% Ham) with margin, giving the optimizer a smooth
// gradient on exotic content. Canonical lumps are unaffected (full strength).
inline constexpr Real EXOTIC_AMP_SCALE = 0.25;

// Effective amplitude actually used for a lump: canonical lumps use their raw
// search amplitude; exotic lumps are damped into the convergent regime.
inline Real effective_amp(const lump_t &L)
{
    return (L.exotic != 0) ? (EXOTIC_AMP_SCALE * L.amp) : L.amp;
}

// Smooth (Cartesian) azimuthal modulation of a lump. m = 0 is axisymmetric (no
// angular momentum possible); m >= 1 breaks axisymmetry so a rigidly rotating
// pattern carries net L_z. Cartesian harmonics avoid the coordinate singularity
// of cos(m*phi) on the z-axis.
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

// Contribution of a single lump to phi (without the constant background).
inline Real lump_phi(const RealVect &loc, const lump_t &L)
{
    if (L.amp == 0.0)
    {
        return 0.0;
    }
    const Real dx  = loc[0] - L.center[0];
    const Real dy  = loc[1] - L.center[1];
    const Real dz  = loc[2] - L.center[2];
    const Real r2  = dx * dx + dy * dy + dz * dz;
    const Real w   = L.width;
    const Real env = exp(-r2 / (2.0 * w * w));
    return effective_amp(L) * angular_factor(L.mode, dx, dy, w) * env;
}

// Finite-difference spatial gradient of a single lump's phi.
inline void lump_grad_phi(const RealVect &loc, const lump_t &L,
                          std::array<Real, 3> &grad)
{
    if (L.amp == 0.0)
    {
        grad = {0.0, 0.0, 0.0};
        return;
    }
    const Real eps = 1.0e-3 * L.width;
    for (int i = 0; i < 3; ++i)
    {
        RealVect lp = loc;
        RealVect lm = loc;
        lp[i] += eps;
        lm[i] -= eps;
        grad[i] = (lump_phi(lp, L) - lump_phi(lm, L)) / (2.0 * eps);
    }
}

// Conjugate momentum contributed by a single lump: the convective derivative of
// a pattern boosted with velocity v and rigidly rotated about z with rate omega
//   Pi = -( v . grad phi ) - omega * d phi / d(azimuth).
// The momentum density S_i = -Pi d_i phi then carries net P_i ~ v_i and
// L_z ~ omega (the latter requires mode >= 1).
inline Real lump_pi(const RealVect &loc, const lump_t &L)
{
    const bool has_kinematics =
        (L.velocity[0] != 0.0 || L.velocity[1] != 0.0 ||
         L.velocity[2] != 0.0 || L.omega != 0.0);
    if (L.amp == 0.0 || !has_kinematics)
    {
        return 0.0;
    }

    std::array<Real, 3> g;
    lump_grad_phi(loc, L, g);

    const Real boost =
        -(L.velocity[0] * g[0] + L.velocity[1] * g[1] + L.velocity[2] * g[2]);

    const Real lx  = loc[0] - L.center[0];
    const Real ly  = loc[1] - L.center[1];
    const Real rot = -L.omega * (lx * g[1] - ly * g[0]);

    return boost + rot;
}

inline void read_params(GRParmParse &pp, params_t &matter_params)
{
    pp.get("phi_0", matter_params.phi_0);
    pp.get("dphi", matter_params.dphi);
    pp.get("dphi_length", matter_params.dphi_length);
    pp.get("pi_0", matter_params.pi_0);
    pp.get("dpi", matter_params.dpi);
    pp.get("dpi_length", matter_params.dpi_length);
    pp.get("scalar_mass", matter_params.scalar_mass);

    matter_params.lumps.clear();

    int num_lumps = 0;
    pp.load("num_lumps", num_lumps, 0);
    if (num_lumps > 0)
    {
        // Multi-lump basis: lump0_*, lump1_*, ...
        for (int k = 0; k < num_lumps; ++k)
        {
            lump_t L;
            read_lump(pp, "lump" + std::to_string(k) + "_", L);
            matter_params.lumps.push_back(L);
        }
    }
    else
    {
        // Backward-compatible single lump (un-indexed lump_* keys). Only added
        // when an amplitude is actually set, so legacy zero-momentum parameter
        // files are unchanged.
        lump_t L;
        pp.load("lump_amp", L.amp, 0.0);
        if (L.amp != 0.0)
        {
            read_lump(pp, "lump_", L);
            matter_params.lumps.push_back(L);
        }
    }
}

}; // namespace MatterParams

#endif
