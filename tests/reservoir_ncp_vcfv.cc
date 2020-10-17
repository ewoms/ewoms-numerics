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
 * \brief Test for the reservoir problem using the NCP model, the VCFV discretization and
 *        finite differences.
 */
#include "config.h"

#include <ewoms/numerics/utils/start.hh>
#include <ewoms/numerics/models/ncp/ncpmodel.hh>
#include <ewoms/numerics/discretizations/vcfv/vcfvdiscretization.hh>
#include "problems/reservoirproblem.hh"


namespace Ewoms {
namespace CO2DefaultTables {
#include <ewoms/material/components/co2tables.inc.cc>
}}

BEGIN_PROPERTIES

NEW_TYPE_TAG(ReservoirNcpVcfvProblem, INHERITS_FROM(NcpModel, ReservoirBaseProblem));

// Select the vertex centered finite volume method as spatial discretization
SET_TAG_PROP(ReservoirNcpVcfvProblem, SpatialDiscretizationSplice, VcfvDiscretization);

// enable the storage cache for this problem so that the storage cache receives wider
// testing
SET_BOOL_PROP(ReservoirNcpVcfvProblem, EnableStorageCache, true);

// reduce the base epsilon for the finite difference method to 10^-11. for some reason
// the simulator converges better with this. (TODO: use automatic differentiation?)
SET_SCALAR_PROP(ReservoirNcpVcfvProblem, BaseEpsilon, 1e-11);

END_PROPERTIES

int main(int argc, char **argv)
{
    using ProblemTypeTag = TTAG(ReservoirNcpVcfvProblem);
    return Ewoms::start<ProblemTypeTag>(argc, argv);
}
