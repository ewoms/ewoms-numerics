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
 * \ingroup FlashModel
 *
 * \brief Declares the properties required by the compositional
 *        multi-phase model based on flash calculations.
 */
#ifndef EWOMS_FLASH_PROPERTIES_HH
#define EWOMS_FLASH_PROPERTIES_HH

#include <ewoms/numerics/common/multiphasebaseproperties.hh>
#include <ewoms/numerics/io/vtkcompositionmodule.hh>
#include <ewoms/numerics/io/vtkenergymodule.hh>
#include <ewoms/numerics/io/vtkdiffusionmodule.hh>

BEGIN_PROPERTIES

//! Provides the thermodynamic relations
NEW_PROP_TAG(FluidSystem);
//! The type of the flash constraint solver
NEW_PROP_TAG(FlashSolver);
//! The maximum accepted error of the flash solver
NEW_PROP_TAG(FlashTolerance);

//! The thermal conduction law which ought to be used
NEW_PROP_TAG(ThermalConductionLaw);
//! The parameters of the thermal conduction law
NEW_PROP_TAG(ThermalConductionLawParams);

//! Specifies whether energy should be considered as a conservation quantity or not
NEW_PROP_TAG(EnableEnergy);
//! Enable diffusive fluxes?
NEW_PROP_TAG(EnableDiffusion);

END_PROPERTIES

#endif
