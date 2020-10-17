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
 * \brief Two-phase test for the immiscible model which uses the
 *        vertex-centered finite volume discretization
 */
#include "config.h"

#include <ewoms/numerics/utils/start.hh>
#include <ewoms/numerics/models/immiscible/immisciblemodel.hh>
#include "problems/lensproblem.hh"

BEGIN_PROPERTIES

NEW_TYPE_TAG(LensProblemVcfvAd, INHERITS_FROM(ImmiscibleTwoPhaseModel, LensBaseProblem));

// use automatic differentiation for this simulator
SET_TAG_PROP(LensProblemVcfvAd, LocalLinearizerSplice, AutoDiffLocalLinearizer);

// use linear finite element gradients if dune-localfunctions is available
#if HAVE_DUNE_LOCALFUNCTIONS
SET_BOOL_PROP(LensProblemVcfvAd, UseP1FiniteElementGradients, true);
#endif

END_PROPERTIES

int main(int argc, char **argv)
{
    using ProblemTypeTag = TTAG(LensProblemVcfvAd);
    return Ewoms::start<ProblemTypeTag>(argc, argv);
}
