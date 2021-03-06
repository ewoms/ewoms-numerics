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
 * \brief Two-phase test for the immiscible model which uses the element-centered finite
 *        volume discretization in conjunction with automatic differentiation
 */
#ifndef EWOMS_LENS_IMMISCIBLE_ECFV_AD_HH
#define EWOMS_LENS_IMMISCIBLE_ECFV_AD_HH

#include <ewoms/numerics/models/immiscible/immisciblemodel.hh>
#include <ewoms/numerics/discretizations/ecfv/ecfvdiscretization.hh>
#include "problems/lensproblem.hh"

BEGIN_PROPERTIES

NEW_TYPE_TAG(LensProblemEcfvAd, INHERITS_FROM(ImmiscibleTwoPhaseModel, LensBaseProblem));

// use the element centered finite volume spatial discretization
SET_TAG_PROP(LensProblemEcfvAd, SpatialDiscretizationSplice, EcfvDiscretization);

// use automatic differentiation for this simulator
SET_TAG_PROP(LensProblemEcfvAd, LocalLinearizerSplice, AutoDiffLocalLinearizer);

// this problem works fine if the linear solver uses single precision scalars
SET_TYPE_PROP(LensProblemEcfvAd, LinearSolverScalar, float);

END_PROPERTIES

#endif // EWOMS_LENS_IMMISCIBLE_ECFV_AD_HH
