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
 * \ingroup NcpModel
 *
 * \brief Declares the properties required for the NCP compositional
 *        multi-phase model.
 */
#ifndef EWOMS_NCP_PROPERTIES_HH
#define EWOMS_NCP_PROPERTIES_HH

#include <ewoms/numerics/common/multiphasebaseproperties.hh>
#include <ewoms/numerics/io/vtkcompositionmodule.hh>
#include <ewoms/numerics/io/vtkenergymodule.hh>
#include <ewoms/numerics/io/vtkdiffusionmodule.hh>

BEGIN_PROPERTIES

//! Enable the energy equation?
NEW_PROP_TAG(EnableEnergy);

//! Enable diffusive fluxes?
NEW_PROP_TAG(EnableDiffusion);

//! The unmodified weight for the pressure primary variable
NEW_PROP_TAG(NcpPressureBaseWeight);
//! The weight for the saturation primary variables
NEW_PROP_TAG(NcpSaturationsBaseWeight);
//! The unmodified weight for the fugacity primary variables
NEW_PROP_TAG(NcpFugacitiesBaseWeight);

//! The themodynamic constraint solver which calculates the
//! composition of any phase given all component fugacities.
NEW_PROP_TAG(NcpCompositionFromFugacitiesSolver);

END_PROPERTIES

#endif
