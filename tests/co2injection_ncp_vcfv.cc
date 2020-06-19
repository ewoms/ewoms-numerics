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
 * \brief Test for the isothermal VCVF discretization based on non-linear
 *        complementarity problems.
 */
#include "config.h"

#include <ewoms/numerics/utils/start.hh>
#include <ewoms/numerics/models/ncp/ncpmodel.hh>
#include <ewoms/numerics/discretizations/vcfv/vcfvdiscretization.hh>

#include "problems/co2injectionproblem.hh"

BEGIN_PROPERTIES

NEW_TYPE_TAG(Co2InjectionNcpVcfvProblem, INHERITS_FROM(NcpModel, Co2InjectionBaseProblem));
SET_TAG_PROP(Co2InjectionNcpVcfvProblem, SpatialDiscretizationSplice, VcfvDiscretization);

END_PROPERTIES

int main(int argc, char **argv)
{
    using VcfvProblemTypeTag = TTAG(Co2InjectionNcpVcfvProblem);
    return Ewoms::start<VcfvProblemTypeTag>(argc, argv);
}
