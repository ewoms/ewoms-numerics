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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \brief Test for the black-oil model using the VCFV discretization.
 */
#include "config.h"

#include <ewoms/numerics/utils/start.hh>
#include <ewoms/numerics/models/blackoil/blackoilmodel.hh>
#include <ewoms/numerics/discretizations/vcfv/vcfvdiscretization.hh>
#include "problems/reservoirproblem.hh"

namespace Ewoms {
namespace CO2DefaultTables {
#include <ewoms/material/components/co2tables.inc.cc>
}}

BEGIN_PROPERTIES

NEW_TYPE_TAG(ReservoirBlackOilVcfvProblem, INHERITS_FROM(BlackOilModel, ReservoirBaseProblem));

// Select the vertex centered finite volume method as spatial discretization
SET_TAG_PROP(ReservoirBlackOilVcfvProblem, SpatialDiscretizationSplice, VcfvDiscretization);

END_PROPERTIES

int main(int argc, char **argv)
{
    using ProblemTypeTag = TTAG(ReservoirBlackOilVcfvProblem);
    return Ewoms::start<ProblemTypeTag>(argc, argv);
}
