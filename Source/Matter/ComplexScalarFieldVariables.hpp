#ifndef COMPLEXSCALARFIELDVARIABLES_HPP_
#define COMPLEXSCALARFIELDVARIABLES_HPP_

#include "MetricVariables.hpp"
#include "ParityDefinitions.hpp"

enum
{
    c_phi_re = NUM_METRIC_VARS,
    c_phi_im,
    c_Pi_re,
    c_Pi_im,
    NUM_MULTIGRID_VARS
};

namespace MatterVariables
{

static const std::array<std::string, NUM_MULTIGRID_VARS - NUM_METRIC_VARS>
    variable_names = {"phi_re", "phi_im", "Pi_re", "Pi_im"};

static constexpr std::array<int, NUM_MULTIGRID_VARS - NUM_METRIC_VARS> const
    vars_parity = {EVEN, EVEN, EVEN, EVEN};

} // namespace MatterVariables

#endif // COMPLEXSCALARFIELDVARIABLES_HPP_
