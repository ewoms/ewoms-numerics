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
 * \brief Declare the basic properties used by the common infrastructure of
 *        the element-centered finite volume discretization.
 */
#ifndef EWOMS_ECFV_PROPERTIES_HH
#define EWOMS_ECFV_PROPERTIES_HH

#include <ewoms/numerics/discretizations/common/fvbaseproperties.hh>

#include <ewoms/common/propertysystem.hh>

BEGIN_PROPERTIES

//! The type tag for models based on the ECFV-scheme
NEW_TYPE_TAG(EcfvDiscretization, INHERITS_FROM(FvBaseDiscretization));

END_PROPERTIES

#endif
