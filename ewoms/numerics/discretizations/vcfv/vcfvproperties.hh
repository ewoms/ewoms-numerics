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
 * \brief Declares the basic properties used by the common infrastructure of
 *        the vertex-centered finite volume discretization.
 */
#ifndef EWOMS_VCFV_PROPERTIES_HH
#define EWOMS_VCFV_PROPERTIES_HH

#include <ewoms/numerics/discretizations/common/fvbaseproperties.hh>

#include <ewoms/common/propertysystem.hh>

BEGIN_PROPERTIES

//! The type tag for models based on the VCFV-scheme
NEW_TYPE_TAG(VcfvDiscretization, INHERITS_FROM(FvBaseDiscretization));

//! Use P1 finite-elements gradients instead of two-point gradients. Note that setting
//! this property to true requires the dune-localfunctions module to be available.
NEW_PROP_TAG(UseP1FiniteElementGradients);

END_PROPERTIES

#endif
