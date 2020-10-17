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
 * \ingroup PvsModel
 *
 * \brief Declares the properties required for the compositional
 *        multi-phase primary variable switching model.
 */
#ifndef EWOMS_PVS_PROPERTIES_HH
#define EWOMS_PVS_PROPERTIES_HH

#include <ewoms/numerics/common/multiphasebaseproperties.hh>
#include <ewoms/numerics/common/diffusionmodule.hh>
#include <ewoms/numerics/common/energymodule.hh>
#include <ewoms/numerics/io/vtkcompositionmodule.hh>
#include <ewoms/numerics/io/vtkphasepresencemodule.hh>
#include <ewoms/numerics/io/vtkdiffusionmodule.hh>
#include <ewoms/numerics/io/vtkenergymodule.hh>

BEGIN_PROPERTIES

//! Specifies whether energy is considered as a conservation quantity or not
NEW_PROP_TAG(EnableEnergy);
//! Enable diffusive fluxes?
NEW_PROP_TAG(EnableDiffusion);

//! The verbosity of the model (0 -> do not print anything, 2 -> spam stdout a lot)
NEW_PROP_TAG(PvsVerbosity);
//! The basis value for the weight of the pressure primary variable
NEW_PROP_TAG(PvsPressureBaseWeight);
//! The basis value for the weight of the saturation primary variables
NEW_PROP_TAG(PvsSaturationsBaseWeight);
//! The basis value for the weight of the mole fraction primary variables
NEW_PROP_TAG(PvsMoleFractionsBaseWeight);

END_PROPERTIES

#endif
