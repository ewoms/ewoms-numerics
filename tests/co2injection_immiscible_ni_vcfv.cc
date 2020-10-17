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
 * \brief Simulation of the injection problem using the VCVF discretization
 *        assuming immisicibility and with energy enabled.
 */
#include "config.h"

#include <ewoms/numerics/utils/start.hh>
#include <ewoms/numerics/models/immiscible/immisciblemodel.hh>
#include <ewoms/numerics/discretizations/vcfv/vcfvdiscretization.hh>
#include "problems/co2injectionproblem.hh"

namespace Ewoms {
namespace Co2Injection {
#include <ewoms/material/components/co2tables.inc.cc>
}}

BEGIN_PROPERTIES

NEW_TYPE_TAG(Co2InjectionImmiscibleNiVcfvProblem, INHERITS_FROM(ImmiscibleModel,
                                                                Co2InjectionBaseProblem));
SET_TAG_PROP(Co2InjectionImmiscibleNiVcfvProblem, SpatialDiscretizationSplice, VcfvDiscretization);

SET_BOOL_PROP(Co2InjectionImmiscibleNiVcfvProblem, EnableEnergy, true);

END_PROPERTIES

////////////////////////
// the main function
////////////////////////
int main(int argc, char **argv)
{
    using VcfvProblemTypeTag = TTAG(Co2InjectionImmiscibleNiVcfvProblem);
    return Ewoms::start<VcfvProblemTypeTag>(argc, argv);
}
