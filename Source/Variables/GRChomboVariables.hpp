/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#ifndef GRCHOMBOVARIABLES_HPP
#define GRCHOMBOVARIABLES_HPP

#include "ArrayTools.hpp"
#include "ParityDefinitions.hpp"

// assign an enum to each variable
enum
{
    c_chi,

    c_h11,
    c_h12,
    c_h13,
    c_h22,
    c_h23,
    c_h33,

    c_K,

    c_A11,
    c_A12,
    c_A13,
    c_A22,
    c_A23,
    c_A33,

    c_Theta,

    c_Gamma1,
    c_Gamma2,
    c_Gamma3,

    c_lapse,

    c_shift1,
    c_shift2,
    c_shift3,

    c_B1,
    c_B2,
    c_B3,

    c_phi, // matter field added
    c_Pi,  //(minus) conjugate momentum

#ifdef USE_COMPLEX_SCALAR_MATTER
    c_phi2, // imaginary part of complex scalar
    c_Pi2,  // (minus) conjugate momentum of imaginary part
#endif

    NUM_GRCHOMBO_VARS
};

namespace GRChomboVariables
{
static constexpr char const *variable_names[NUM_GRCHOMBO_VARS] = {
    "chi",

    "h11",    "h12",    "h13",    "h22", "h23", "h33",

    "K",

    "A11",    "A12",    "A13",    "A22", "A23", "A33",

    "Theta",

    "Gamma1", "Gamma2", "Gamma3",

    "lapse",

    "shift1", "shift2", "shift3",

    "B1",     "B2",     "B3",

    "phi",    "Pi"
#ifdef USE_COMPLEX_SCALAR_MATTER
    ,         "phi2",   "Pi2"
#endif
};

static constexpr std::array<int, NUM_GRCHOMBO_VARS> const vars_parity = {
    EVEN,   EVEN,   ODD_XY, ODD_XZ, EVEN,  ODD_YZ, EVEN,  EVEN,  EVEN,
    ODD_XY, ODD_XZ, EVEN,   ODD_YZ, EVEN,  EVEN,   ODD_X, ODD_Y, ODD_Z,
    EVEN,   ODD_X,  ODD_Y,  ODD_Z,  ODD_X, ODD_Y,  ODD_Z, EVEN,  EVEN
#ifdef USE_COMPLEX_SCALAR_MATTER
    ,       EVEN,   EVEN
#endif
};

} // namespace GRChomboVariables

#endif /* GRCHOMBOVARIABLES_HPP */