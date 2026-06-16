#ifndef BOSONSTARPARAMS_HPP_
#define BOSONSTARPARAMS_HPP_

#include "GRParmParse.hpp"
#include "REAL.H"
#include <cmath>

namespace BosonStarParams
{

struct params_t
{
    Real scalar_mass   = 0.1;
    Real scalar_lambda = 0.0;
    Real phi_c         = 0.08;
    Real profile_width = 8.0;
    Real omega         = 0.0; // if <= 0, defaults to scalar_mass at paint time
    Real sign          = 1.0; // +1 canonical, -1 phantom (flips T_ab)
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

inline void read_params(GRParmParse &pp, params_t &p)
{
    pp.load("scalar_mass", p.scalar_mass, 0.1);
    pp.load("scalar_lambda", p.scalar_lambda, 0.0);
    pp.load("scalar_sign", p.sign, 1.0);
    pp.load("bs_phi_c", p.phi_c, 0.08);
    pp.load("bs_profile_width", p.profile_width, 8.0);
    pp.load("bs_omega", p.omega, 0.0);
}

} // namespace BosonStarParams

#endif /* BOSONSTARPARAMS_HPP_ */
