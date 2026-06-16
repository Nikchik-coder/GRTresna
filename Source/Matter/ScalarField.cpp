/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

// The boson-star example selects the complex matter model and provides its own
// MatterParams/MultigridVariables, so the canonical real ScalarField source is
// compiled out of that build to avoid duplicate/incompatible symbols.
#ifndef USE_COMPLEX_SCALAR_MATTER

#include "ScalarField.hpp"
#include "DerivativeOperators.hpp"
#include "EMTensor.hpp"
#include "FArrayBox.H"
#include "GRParmParse.hpp"
#include "Grids.hpp"
#include "IntVect.H"
#include "LevelData.H"
#include "MultigridVariables.hpp"
#include "PsiAndAijFunctions.hpp"
#include "REAL.H"
#include "RealVect.H"
#include "Tensor.hpp"
#include <array>
#include <cmath>

void ScalarField::initialise_matter_vars(LevelData<FArrayBox> &a_multigrid_vars,
                                         const RealVect &a_dx) const
{
    CH_assert(a_multigrid_vars.nComp() == NUM_MULTIGRID_VARS);

    DataIterator dit = a_multigrid_vars.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        // These contain the vars in the boxes, set them all to zero
        FArrayBox &multigrid_vars_box = a_multigrid_vars[dit()];

        // Iterate over the box and set non zero comps
        Box ghosted_box = multigrid_vars_box.box();
        BoxIterator bit(ghosted_box);
        for (bit.begin(); bit.ok(); ++bit)
        {

            // work out location on the grid
            IntVect iv = bit();
            RealVect loc;
            Grids::get_loc(loc, iv, a_dx, center);

            multigrid_vars_box(iv, c_phi_0) = my_phi_function(loc);
            multigrid_vars_box(iv, c_Pi_0) = my_Pi_function(loc);
        }
    }
}

// template <class data_t>
emtensor_t ScalarField::compute_emtensor(const IntVect a_iv,
                                         const RealVect &a_dx,
                                         FArrayBox &a_multigrid_vars_box) const
{
    emtensor_t out;

    RealVect loc;
    Grids::get_loc(loc, a_iv, a_dx, center);

    Real psi_reg = a_multigrid_vars_box(a_iv, c_psi_reg);
    Real psi_bh  = psi_and_Aij_functions->compute_bowenyork_psi(loc);
    Real psi_0   = psi_reg + psi_bh;
    Real phi_0   = a_multigrid_vars_box(a_iv, c_phi_0);
    const Real chi = pow(psi_0, -4.0);

    // Independent-field matter model. The background radial profile and each
    // lump are assembled as separate scalar fields (cross terms dropped), so a
    // lump flagged EXOTIC contributes its kinetic-energy and momentum density
    // with a flipped sign -- a genuine ghost/phantom source of negative energy
    // (NEC violation) that the constraint solve bakes into the initial
    // geometry. With no exotic lumps and well-separated clouds this matches the
    // canonical single-field source up to the (small) inter-lump cross terms.
    Real rho_kinetic = 0.0;
    Tensor<1, Real, SpaceDim> Si;
    FOR1(i) { Si[i] = 0.0; }

    // --- Background (always canonical): spherical phi_0 + dphi e^{-r/L} cloud,
    //     carrying no net momentum.
    {
        Real rr = sqrt(loc[0] * loc[0] + loc[1] * loc[1] + loc[2] * loc[2]);
        Real env = m_matter_params.dphi *
                   exp(-rr / m_matter_params.dphi_length);
        Real Pi_bg = m_matter_params.pi_0 +
                     m_matter_params.dpi *
                         exp(-rr / m_matter_params.dpi_length);
        Tensor<1, Real, SpaceDim> d1_bg;
        Real radial = (rr > 1.0e-12)
                          ? (-env / (m_matter_params.dphi_length * rr))
                          : 0.0;
        Real d1_bg_squared = 0.0;
        FOR1(i)
        {
            d1_bg[i] = radial * loc[i];
            d1_bg_squared += d1_bg[i] * d1_bg[i];
        }
        rho_kinetic += 0.5 * chi * d1_bg_squared + 0.5 * Pi_bg * Pi_bg;
        FOR1(i) { Si[i] += -Pi_bg * d1_bg[i]; }
    }

    // --- Lumps (each canonical or exotic).
    for (const auto &L : m_matter_params.lumps)
    {
        if (L.amp == 0.0)
        {
            continue;
        }
        const Real sign = (L.exotic != 0) ? -1.0 : 1.0;
        std::array<Real, 3> g;
        MatterParams::lump_grad_phi(loc, L, g);
        Real Pi_k          = MatterParams::lump_pi(loc, L);
        Real d1_k_squared = 0.0;
        FOR1(i) { d1_k_squared += g[i] * g[i]; }
        rho_kinetic += sign * (0.5 * chi * d1_k_squared + 0.5 * Pi_k * Pi_k);
        FOR1(i) { Si[i] += sign * (-Pi_k * g[i]); }
    }

    Real V_of_phi = my_potential_function(phi_0);

    out.rho = rho_kinetic + V_of_phi;
    FOR1(i) { out.Si[i] = Si[i]; }

    return out;
}

#endif /* USE_COMPLEX_SCALAR_MATTER */
