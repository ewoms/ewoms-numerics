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
 * \brief Test for the Richards model using the VCFV discretization.
 */
#include "config.h"

#include <ewoms/numerics/utils/start.hh>
#include <ewoms/numerics/discretizations/vcfv/vcfvdiscretization.hh>

#include "problems/richardslensproblem.hh"

BEGIN_PROPERTIES

NEW_TYPE_TAG(RichardsLensVcfvProblem, INHERITS_FROM(RichardsLensProblem));
SET_TAG_PROP(RichardsLensVcfvProblem, SpatialDiscretizationSplice, VcfvDiscretization);

END_PROPERTIES

int main(int argc, char **argv)
{
    using ProblemTypeTag = TTAG(RichardsLensVcfvProblem);
    return Ewoms::start<ProblemTypeTag>(argc, argv);
}
