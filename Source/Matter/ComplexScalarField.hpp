#ifndef COMPLEXSCALARFIELD_HPP_
#define COMPLEXSCALARFIELD_HPP_

#include "BosonStarParams.hpp"
#include "EMTensor.hpp"
#include "FArrayBox.H"
#include "GRParmParse.hpp"
#include "IntVect.H"
#include "LevelData.H"
#include "PsiAndAijFunctions.hpp"
#include "REAL.H"
#include "RealVect.H"
#include "Tensor.hpp"

class ComplexScalarField
{
  public:
    using params_t = BosonStarParams::params_t;

    ComplexScalarField(params_t a_params,
                       PsiAndAijFunctions *a_psi_and_Aij_functions,
                       const std::array<double, SpaceDim> a_center,
                       RealVect a_domainLength)
        : m_params(a_params),
          psi_and_Aij_functions(a_psi_and_Aij_functions), center(a_center),
          domainLength(a_domainLength)
    {
    }

    emtensor_t compute_emtensor(const IntVect a_iv, const RealVect &a_dx,
                                FArrayBox &a_multigrid_vars_box) const;

    void initialise_matter_vars(LevelData<FArrayBox> &a_multigrid_vars,
                                const RealVect &a_dx) const;

    static void read_params(GRParmParse &pp, params_t &matter_params)
    {
        BosonStarParams::read_params(pp, matter_params);
    }

    params_t m_params;
    PsiAndAijFunctions *psi_and_Aij_functions;

  private:
    const std::array<double, SpaceDim> center;
    RealVect domainLength;
};

#endif /* COMPLEXSCALARFIELD_HPP_ */
