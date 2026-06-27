#ifndef BOSONSTARPARAMS_HPP_
#define BOSONSTARPARAMS_HPP_

#include "GRParmParse.hpp"
#include "REAL.H"
#include "RealVect.H"
#include <array>
#include <cmath>
#include <string>
#include <vector>

namespace BosonStarParams
{

// One site in a bosonic shell: real part of Phi plus optional shell kinematics.
struct boson_lump_t
{
    Real amp                     = 0.0;
    Real width                   = 5.0;
    std::array<Real, 3> center   = {0.0, 0.0, 0.0};
    std::array<Real, 3> velocity = {0.0, 0.0, 0.0};
    Real omega                   = 0.0;
    int mode                     = 0;
    int profile                  = 0; // 0 = Gaussian, 1 = smoothed top-hat
    // Phantom/ghost flag. 0 => canonical (positive-energy) complex lump; != 0
    // => EXOTIC: this lump's kinetic + gradient + potential energy density and
    // its momentum density enter the constraint solve with a flipped sign, a
    // genuine ghost source of negative energy (NEC violation).  Treated as an
    // independent field, so a config can mix normal and exotic lumps -- this is
    // what the warp/wormhole FTL geometry needs (mirrors ScalarFieldBH).
    int exotic                   = 0;
};

struct params_t
{
    Real scalar_mass   = 0.1;
    Real scalar_lambda = 0.0;
    Real phi_c         = 0.08;
    Real profile_width = 8.0;
    Real omega         = 0.0; // global U(1) phase velocity; <=0 => scalar_mass
    Real sign          = 1.0;

    // Multi-site bosonic shell. Empty => single centered Gaussian (legacy).
    std::vector<boson_lump_t> lumps;
};

inline Real phi0_profile(Real r, Real phi_c, Real width)
{
    const Real rr = r * r;
    const Real w2 = width * width;
    return phi_c * std::exp(-0.5 * rr / w2);
}

inline Real dphi0_dr(Real r, Real phi_c, Real width)
{
    if (r <= 1.0e-12)
        return 0.0;
    return -phi0_profile(r, phi_c, width) * r / (width * width);
}

inline Real angular_factor(int m, Real dx, Real dy, Real w)
{
    switch (m)
    {
    case 1:
        return dx / w;
    case 2:
        return (dx * dx - dy * dy) / (w * w);
    default:
        return 1.0;
    }
}

inline Real lump_envelope(Real r2, Real r, Real w, int profile)
{
    if (profile == 1)
    {
        const Real soft = 0.25 * w;
        return 0.5 * (1.0 - std::tanh((r - w) / soft));
    }
    if (profile == 2)
    {
        // Bound boson lump: sech(r / w) with w = 1/sqrt(m^2 - omega^2).  Regular
        // at r=0 and decays as 2 e^{-r/w} -- the correct exponential tail of a
        // massive bound scalar.  Must match the Python painter and the GRTeclyn
        // pump controller (PROFILE_SECH_BOUND) so the solve, the initial data
        // and the trap target are consistent.
        return 1.0 / std::cosh(r / w);
    }
    return std::exp(-0.5 * r2 / (w * w));
}

inline Real lump_phi1(const RealVect &loc, const boson_lump_t &L)
{
    if (L.amp == 0.0)
        return 0.0;
    const Real dx = loc[0] - L.center[0];
    const Real dy = loc[1] - L.center[1];
    const Real dz = loc[2] - L.center[2];
    const Real r2 = dx * dx + dy * dy + dz * dz;
    const Real r  = std::sqrt(r2);
    const Real w  = L.width;
    return L.amp * angular_factor(L.mode, dx, dy, w) *
           lump_envelope(r2, r, w, L.profile);
}

inline void lump_grad_phi1(const RealVect &loc, const boson_lump_t &L,
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
        grad[i] = (lump_phi1(lp, L) - lump_phi1(lm, L)) / (2.0 * eps);
    }
}

// Real-channel Pi from shell boost / rigid rotation (same as ScalarFieldBH).
inline Real lump_pi1(const RealVect &loc, const boson_lump_t &L)
{
    const bool has_kinematics =
        (L.velocity[0] != 0.0 || L.velocity[1] != 0.0 ||
         L.velocity[2] != 0.0 || L.omega != 0.0);
    if (L.amp == 0.0 || !has_kinematics)
        return 0.0;

    std::array<Real, 3> g;
    lump_grad_phi1(loc, L, g);

    const Real boost =
        -(L.velocity[0] * g[0] + L.velocity[1] * g[1] + L.velocity[2] * g[2]);

    const Real lx  = loc[0] - L.center[0];
    const Real ly  = loc[1] - L.center[1];
    const Real rot = -L.omega * (lx * g[1] - ly * g[0]);

    return boost + rot;
}

inline Real total_phi1(const RealVect &loc, const params_t &p)
{
    if (!p.lumps.empty())
    {
        Real sum = 0.0;
        for (const auto &L : p.lumps)
            sum += lump_phi1(loc, L);
        return sum;
    }
    const Real r = std::sqrt(loc[0] * loc[0] + loc[1] * loc[1] +
                             loc[2] * loc[2]);
    return phi0_profile(r, p.phi_c, p.profile_width);
}

inline Real total_pi1(const RealVect &loc, const params_t &p)
{
    if (p.lumps.empty())
        return 0.0;
    Real sum = 0.0;
    for (const auto &L : p.lumps)
        sum += lump_pi1(loc, L);
    return sum;
}

inline void total_grad_phi1(const RealVect &loc, const params_t &p,
                            std::array<Real, 3> &grad)
{
    if (!p.lumps.empty())
    {
        grad = {0.0, 0.0, 0.0};
        for (const auto &L : p.lumps)
        {
            std::array<Real, 3> g;
            lump_grad_phi1(loc, L, g);
            for (int i = 0; i < 3; ++i)
                grad[i] += g[i];
        }
        return;
    }
    const Real r = std::sqrt(loc[0] * loc[0] + loc[1] * loc[1] +
                             loc[2] * loc[2]);
    if (r > 1.0e-12)
    {
        const Real dphidr = dphi0_dr(r, p.phi_c, p.profile_width);
        for (int i = 0; i < 3; ++i)
            grad[i] = dphidr * loc[i] / r;
    }
    else
    {
        grad = {0.0, 0.0, 0.0};
    }
}

inline void read_boson_lump(GRParmParse &pp, const std::string &prefix,
                            boson_lump_t &L)
{
    pp.load((prefix + "amp").c_str(), L.amp, 0.0);
    pp.load((prefix + "width").c_str(), L.width, 5.0);
    pp.load((prefix + "omega").c_str(), L.omega, 0.0);
    pp.load((prefix + "mode").c_str(), L.mode, 0);
    pp.load((prefix + "profile").c_str(), L.profile, 0);
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

inline void read_params(GRParmParse &pp, params_t &p)
{
    pp.load("scalar_mass", p.scalar_mass, 0.1);
    pp.load("scalar_lambda", p.scalar_lambda, 0.0);
    pp.load("scalar_sign", p.sign, 1.0);
    pp.load("bs_phi_c", p.phi_c, 0.08);
    pp.load("bs_profile_width", p.profile_width, 8.0);
    pp.load("bs_omega", p.omega, 0.0);

    p.lumps.clear();
    int num_lumps = 0;
    pp.load("num_lumps", num_lumps, 0);
    if (num_lumps > 0)
    {
        for (int k = 0; k < num_lumps; ++k)
        {
            boson_lump_t L;
            read_boson_lump(pp, "lump" + std::to_string(k) + "_", L);
            if (L.amp != 0.0)
                p.lumps.push_back(L);
        }
    }
}

} // namespace BosonStarParams

#endif /* BOSONSTARPARAMS_HPP_ */
