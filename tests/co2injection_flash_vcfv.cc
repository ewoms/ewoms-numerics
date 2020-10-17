// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the eWoms project.

  eWoms is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  eWoms is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with eWoms.  If not, see <http://www.gnu.org/licenses/>.
*/
/*!
 * \file
 *
 * \brief Test for the isothermal compositional model based on flash
 *        calculations.
 */
#include "config.h"

#if HAVE_QUAD
#include <ewoms/common/quad.hh>
#endif

#include <ewoms/numerics/utils/start.hh>
#include <ewoms/numerics/models/flash/flashmodel.hh>
#include <ewoms/numerics/discretizations/vcfv/vcfvdiscretization.hh>
#include "problems/co2injectionflash.hh"
#include "problems/co2injectionproblem.hh"

namespace Ewoms {
namespace Co2Injection {
#include <ewoms/material/components/co2tables.inc.cc>
}}

BEGIN_PROPERTIES

NEW_TYPE_TAG(Co2InjectionFlashVcfvProblem, INHERITS_FROM(FlashModel, Co2InjectionBaseProblem));
SET_TAG_PROP(Co2InjectionFlashVcfvProblem, SpatialDiscretizationSplice, VcfvDiscretization);

// use the flash solver adapted to the CO2 injection problem
SET_TYPE_PROP(
    Co2InjectionFlashVcfvProblem, FlashSolver,
    Ewoms::Co2InjectionFlash<GET_PROP_TYPE(TypeTag, Scalar),
                           GET_PROP_TYPE(TypeTag, FluidSystem)>);

// the flash model has serious problems with the numerical
// precision. if quadruple precision math is available, we use it,
// else we increase the tolerance of the Newton solver
#if HAVE_QUAD
SET_TYPE_PROP(Co2InjectionFlashVcfvProblem, Scalar, quad);

// the default linear solver used for this problem (-> AMG) cannot be used with quadruple
// precision scalars... (this seems to only apply to Dune >= 2.4)
SET_TAG_PROP(Co2InjectionFlashVcfvProblem, LinearSolverSplice, ParallelBiCGStabLinearSolver);
#else
SET_SCALAR_PROP(Co2InjectionFlashVcfvProblem, NewtonTolerance, 1e-5);
#endif

END_PROPERTIES

int main(int argc, char **argv)
{
    using VcfvProblemTypeTag = TTAG(Co2InjectionFlashVcfvProblem);
    return Ewoms::start<VcfvProblemTypeTag>(argc, argv);
}
