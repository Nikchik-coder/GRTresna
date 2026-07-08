#ifndef BOSONSTARPARAMS_HPP_
#define BOSONSTARPARAMS_HPP_

#include "GRParmParse.hpp"
#include "REAL.H"
#include "RealVect.H"
#include <array>
#include <cmath>
#include <fstream>
#include <mutex>
#include <sstream>
#include <string>
#include <unordered_map>
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
    // 0 = Gaussian, 1 = tanh shell, 2 = sech bound lump, 3 = tabulated Q-ball
    // ODE profile phi0(r) loaded from ``profile_path`` (see qball_load_table).
    int profile                  = 0;
    // Path to the tabulated flat-space Q-ball radial profile (two whitespace-
    // separated columns: r  phi0) used when profile == 3.  Written by the Python
    // qball_radial_ode solver so the GRTresna solve and the gridinit repaint use
    // the SAME phi0(r); empty unless profile == 3.
    std::string profile_path     = "";
    // Phantom/ghost flag. 0 => canonical (positive-energy) complex lump; != 0
    // => EXOTIC: this lump's kinetic + gradient + potential energy density and
    // its momentum density enter the constraint solve with a flipped sign, a
    // genuine ghost source of negative energy (NEC violation).  Treated as an
    // independent field, so a config can mix normal and exotic lumps -- this is
    // what the warp/wormhole FTL geometry needs (mirrors ScalarFieldBH).
    int exotic                   = 0;
    // Phase-winding flag for a rotating (Kleihaus-Kunz / Teo-class) wormhole.
    // 0 => legacy boson-star ansatz (phi2 = 0, real angular_factor modulation of
    // phi1, U(1) phase only in Pi_im).  != 0 => genuine phase winding
    //   Phi = f(r,theta) e^{i(m phi_az)},  f = amp * env(r) * (sin theta)^m,
    //   phi1 = f cos(m phi_az),  phi2 = f sin(m phi_az),
    //   Pi1 = (omega/alpha) phi2,  Pi2 = -(omega/alpha) phi1
    // (the e^{-i omega t} stationary rotation, matching the GRTeclyn complex_scalar
    // evolution convention).  ``mode`` is reused as the azimuthal winding number m.
    // This is what makes |Phi|^2 axisymmetric (no four-lobe dispersal, lesson L3)
    // and lets the two-channel momentum density source a clean J_z.
    int winding                  = 0;
};

struct params_t
{
    Real scalar_mass   = 0.1;
    Real scalar_lambda = 0.0;
    // Sextic self-interaction coupling.  V = 1/2 m^2|Phi|^2 - 1/4 lambda|Phi|^4
    // + 1/6 mu |Phi|^6.  Required (mu > 0 with lambda > 0) for a stable 3D
    // Q-ball; threaded into ComplexScalarField::potential_value so the solved
    // metric matches the GRTeclyn evolution potential.
    Real scalar_mu     = 0.0;
    Real phi_c         = 0.08;
    Real profile_width = 8.0;
    Real omega         = 0.0; // global U(1) phase velocity; <=0 => scalar_mass
    Real sign          = 1.0;
    // Shared tabulated Q-ball profile path for profile == 3 lumps (copied onto
    // each lump that does not set its own ``lump{k}_profile_path``).
    std::string qball_profile_path = "";

    // Multi-site bosonic shell. Empty => single centered Gaussian (legacy).
    std::vector<boson_lump_t> lumps;
};

// ---------------------------------------------------------------------------
// Tabulated flat-space Q-ball radial profile phi0(r) (lump profile == 3).
//
// The Python qball_radial_ode solver shoots the stationary ODE
//   phi0'' + (2/r) phi0' = (m^2 - omega^2) phi0 - lambda phi0^3 + mu phi0^5
// and writes the result to a two-column ``.dat`` file.  Loading that exact
// table here (instead of re-solving in C++) guarantees the constraint solve and
// the gridinit repaint paint the SAME phi0(r) -- the previous sech-in-solve /
// ODE-in-paint mismatch produced inconsistent constraint data and garbage
// fields.  The base table is stored at its own core amplitude phi_c = phi0(0);
// per-lump amplitude scaling (amp / phi_c) is applied in lump_phi1, matching the
// Python painter's profile_for_lump scaling.
// ---------------------------------------------------------------------------
struct qball_table_t
{
    std::vector<Real> r;
    std::vector<Real> phi0;
    // Optional lapse column alpha(r).  Present for the self-gravitating boson
    // star (3-column table); empty for the flat-space Q-ball (2-column table),
    // in which case alpha is treated as 1 everywhere (flat space).
    std::vector<Real> alpha;
    Real phi_c = 0.0; // phi0[0]
};

inline const qball_table_t &qball_load_table(const std::string &path)
{
    static std::unordered_map<std::string, qball_table_t> cache;
    static std::mutex cache_mutex;
    std::lock_guard<std::mutex> lock(cache_mutex);

    auto found = cache.find(path);
    if (found != cache.end())
        return found->second;

    qball_table_t table;
    if (!path.empty())
    {
        std::ifstream file(path);
        std::string line;
        // Parse line-by-line so that '#' comment / header lines are skipped
        // (a bare ``file >> r >> phi`` aborts on the first '#' and silently
        // yields an EMPTY table -> sech fallback -- the bug that made the
        // tabulated ODE profile never load).  Each data line is
        // ``r  phi0`` or ``r  phi0  alpha``.
        while (std::getline(file, line))
        {
            // Strip leading whitespace to detect comments/blank lines.
            std::size_t first = line.find_first_not_of(" \t\r\n");
            if (first == std::string::npos || line[first] == '#')
                continue;
            std::istringstream iss(line);
            Real rr = 0.0, pp = 0.0, aa = 0.0;
            if (!(iss >> rr >> pp))
                continue;
            table.r.push_back(rr);
            table.phi0.push_back(pp);
            if (iss >> aa)
                table.alpha.push_back(aa);
        }
    }
    if (!table.phi0.empty())
        table.phi_c = table.phi0.front();
    // Only keep the lapse column if it covers every row (otherwise it is
    // inconsistent and we fall back to alpha = 1).
    if (table.alpha.size() != table.phi0.size())
        table.alpha.clear();

    auto inserted = cache.emplace(path, std::move(table));
    return inserted.first->second;
}

// Linear interpolation of the (optional) lapse column alpha(r).  Returns 1.0
// when the table has no lapse column (flat-space Q-ball) or is empty.
inline Real qball_alpha_table_interp(const qball_table_t &table, Real r)
{
    const std::size_t n = table.alpha.size();
    if (n == 0)
        return 1.0;
    if (r <= table.r.front())
        return table.alpha.front();
    if (r >= table.r.back())
        return table.alpha.back();

    std::size_t lo = 0, hi = n - 1;
    while (hi - lo > 1)
    {
        const std::size_t mid = (lo + hi) / 2;
        if (table.r[mid] <= r)
            lo = mid;
        else
            hi = mid;
    }
    const Real span = table.r[hi] - table.r[lo];
    const Real t    = (span > 0.0) ? (r - table.r[lo]) / span : 0.0;
    return table.alpha[lo] + t * (table.alpha[hi] - table.alpha[lo]);
}

inline Real qball_table_interp(const qball_table_t &table, Real r)
{
    const std::size_t n = table.r.size();
    if (n == 0)
        return 0.0;
    if (r <= table.r.front())
        return table.phi0.front();
    if (r >= table.r.back())
        return table.phi0.back();

    std::size_t lo = 0, hi = n - 1;
    while (hi - lo > 1)
    {
        const std::size_t mid = (lo + hi) / 2;
        if (table.r[mid] <= r)
            lo = mid;
        else
            hi = mid;
    }
    const Real span = table.r[hi] - table.r[lo];
    const Real t    = (span > 0.0) ? (r - table.r[lo]) / span : 0.0;
    return table.phi0[lo] + t * (table.phi0[hi] - table.phi0[lo]);
}

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

// Normalised tabulated Q-ball envelope: phi0(r) / phi0(0) so envelope(0) = 1,
// matching the other lump_envelope shapes.  lump_phi1 then multiplies by L.amp,
// reproducing the Python painter's (amp / phi_c) * phi0(r) scaling exactly.
// Falls back to a sech of width L.width if the table is missing/empty so a
// misconfigured run degrades gracefully rather than painting zero field.
inline Real qball_envelope(Real r, const boson_lump_t &L)
{
    const qball_table_t &table = qball_load_table(L.profile_path);
    if (table.phi0.empty() || table.phi_c <= 0.0)
    {
        const Real w = (L.width > 0.0) ? L.width : 1.0;
        return 1.0 / std::cosh(r / w);
    }
    return qball_table_interp(table, r) / table.phi_c;
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
    const Real env =
        (L.profile == 3) ? qball_envelope(r, L) : lump_envelope(r2, r, w, L.profile);
    return L.amp * angular_factor(L.mode, dx, dy, w) * env;
}

// Stationary boson-star lapse alpha(r) at this location, from the lump's
// tabulated ODE profile (profile == 3 with a 3rd column).  Returns 1 for any
// other profile / a table without a lapse column (flat space).  A boson star is
// stationary, NOT static: alpha dips well below 1 in the core, and the U(1)
// momentum Pi_im = -(omega/alpha) phi0 must use it or the seed is not an
// equilibrium and radiates.
inline Real lump_alpha(const RealVect &loc, const boson_lump_t &L)
{
    if (L.profile != 3)
        return 1.0;
    const Real dx = loc[0] - L.center[0];
    const Real dy = loc[1] - L.center[1];
    const Real dz = loc[2] - L.center[2];
    const Real r  = std::sqrt(dx * dx + dy * dy + dz * dz);
    return qball_alpha_table_interp(qball_load_table(L.profile_path), r);
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

// ---------------------------------------------------------------------------
// Phase-winding (rotating wormhole) ansatz helpers.  Active per-lump when
// ``L.winding != 0``; ``L.mode`` is the azimuthal winding number m.
// ---------------------------------------------------------------------------

// (sin theta)^m polar factor, sin theta = sqrt(dx^2+dy^2)/r.  Regular on the
// z-axis (returns 0 for m>=1 there, so the toroidal field vanishes smoothly).
inline Real polar_factor(int m, Real dx, Real dy, Real dz)
{
    if (m <= 0)
        return 1.0;
    const Real rho2 = dx * dx + dy * dy;
    const Real r    = std::sqrt(rho2 + dz * dz);
    if (r <= 1.0e-12)
        return 0.0;
    const Real sinth = std::sqrt(rho2) / r;
    Real f           = 1.0;
    for (int k = 0; k < m; ++k)
        f *= sinth;
    return f;
}

// Axisymmetric modulus f(r,theta) = amp * env(r) * (sin theta)^m.
inline Real lump_winding_modulus(const RealVect &loc, const boson_lump_t &L)
{
    if (L.amp == 0.0)
        return 0.0;
    const Real dx = loc[0] - L.center[0];
    const Real dy = loc[1] - L.center[1];
    const Real dz = loc[2] - L.center[2];
    const Real r2 = dx * dx + dy * dy + dz * dz;
    const Real r  = std::sqrt(r2);
    const Real w  = L.width;
    const Real env =
        (L.profile == 3) ? qball_envelope(r, L) : lump_envelope(r2, r, w, L.profile);
    return L.amp * env * polar_factor(L.mode, dx, dy, dz);
}

// phi1 = f cos(m phi_az), phi2 = f sin(m phi_az).
inline void lump_phi_winding(const RealVect &loc, const boson_lump_t &L,
                             Real &phi1, Real &phi2)
{
    const Real f  = lump_winding_modulus(loc, L);
    const Real dx = loc[0] - L.center[0];
    const Real dy = loc[1] - L.center[1];
    const Real az = std::atan2(dy, dx);
    const Real m  = static_cast<Real>(L.mode);
    phi1          = f * std::cos(m * az);
    phi2          = f * std::sin(m * az);
}

// Finite-difference gradients of the winding fields phi1, phi2.
inline void lump_grad_phi_winding(const RealVect &loc, const boson_lump_t &L,
                                  std::array<Real, 3> &g1,
                                  std::array<Real, 3> &g2)
{
    if (L.amp == 0.0)
    {
        g1 = {0.0, 0.0, 0.0};
        g2 = {0.0, 0.0, 0.0};
        return;
    }
    const Real eps = 1.0e-3 * L.width;
    for (int i = 0; i < 3; ++i)
    {
        RealVect lp = loc;
        RealVect lm = loc;
        lp[i] += eps;
        lm[i] -= eps;
        Real p1p, p2p, p1m, p2m;
        lump_phi_winding(lp, L, p1p, p2p);
        lump_phi_winding(lm, L, p1m, p2m);
        g1[i] = (p1p - p1m) / (2.0 * eps);
        g2[i] = (p2p - p2m) / (2.0 * eps);
    }
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

// Imaginary-channel (U(1) phase) momentum Pi_im for a stationary boson star:
//   Pi_im = -(omega / alpha(r)) * phi1
// with alpha(r) the star's lapse (== 1 for a flat-space Q-ball table, which
// recovers the previous Pi_im = -omega phi1).  omega is the global U(1) phase
// velocity (bs_omega; falls back to the scalar mass when unset).
inline Real total_pi2(const RealVect &loc, const params_t &p, Real omega)
{
    if (!p.lumps.empty())
    {
        Real sum = 0.0;
        for (const auto &L : p.lumps)
        {
            const Real a = lump_alpha(loc, L);
            sum += -(omega / a) * lump_phi1(loc, L);
        }
        return sum;
    }
    // Legacy single centered star (no lumps): alpha = 1 (Gaussian ansatz).
    return -omega * total_phi1(loc, p);
}

// Stationary-lapse superposition for the initial gauge:
//   alpha(loc) = 1 + sum_k (alpha_k(loc) - 1)
// which is exactly alpha(r) for a single isolated star and -> 1 far away.
// Setting the output lapse to this (instead of the flat alpha = 1) starts the
// evolution in the boson star's own gauge and avoids a large 1+log transient.
inline Real total_alpha(const RealVect &loc, const params_t &p)
{
    if (p.lumps.empty())
        return 1.0;
    Real alpha = 1.0;
    for (const auto &L : p.lumps)
        alpha += (lump_alpha(loc, L) - 1.0);
    return (alpha > 1.0e-3) ? alpha : 1.0e-3;
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
    pp.load((prefix + "winding").c_str(), L.winding, 0);
    // Optional per-lump tabulated-profile path (profile == 3).  Usually left
    // unset and inherited from the global ``qball_profile_path`` in read_params.
    pp.load((prefix + "profile_path").c_str(), L.profile_path, std::string(""));
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
    pp.load("scalar_mu", p.scalar_mu, 0.0);
    pp.load("scalar_sign", p.sign, 1.0);
    pp.load("bs_phi_c", p.phi_c, 0.08);
    pp.load("bs_profile_width", p.profile_width, 8.0);
    pp.load("bs_omega", p.omega, 0.0);
    pp.load("qball_profile_path", p.qball_profile_path, std::string(""));

    p.lumps.clear();
    int num_lumps = 0;
    pp.load("num_lumps", num_lumps, 0);
    if (num_lumps > 0)
    {
        for (int k = 0; k < num_lumps; ++k)
        {
            boson_lump_t L;
            read_boson_lump(pp, "lump" + std::to_string(k) + "_", L);
            // Inherit the shared tabulated profile path when a profile==3 lump
            // does not specify its own, so all Q-ball lumps read the same phi0(r).
            if (L.profile == 3 && L.profile_path.empty())
                L.profile_path = p.qball_profile_path;
            if (L.amp != 0.0)
                p.lumps.push_back(L);
        }
    }
}

} // namespace BosonStarParams

#endif /* BOSONSTARPARAMS_HPP_ */
