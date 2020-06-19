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
 * \copydoc Ewoms::ImmiscibleIndices
 */
#ifndef EWOMS_IMMISCIBLE_INDICES_HH
#define EWOMS_IMMISCIBLE_INDICES_HH

#include "immiscibleproperties.hh"
#include <ewoms/numerics/common/energymodule.hh>

namespace Ewoms {

/*!
 * \ingroup ImmiscibleModel
 *
 * \brief The indices for the isothermal multi-phase model.
 */
template <class TypeTag, int PVOffset>
struct ImmiscibleIndices
    : public EnergyIndices<PVOffset + GET_PROP_VALUE(TypeTag, NumPhases),
                           GET_PROP_VALUE(TypeTag, EnableEnergy)>
{
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    using EnergyIndices = Ewoms::EnergyIndices<PVOffset + numPhases, enableEnergy>;

public:
    // number of equations/primary variables
    static const int numEq = numPhases + EnergyIndices::numEq_;

    // Primary variable indices

    //! Index for wetting/non-wetting phase pressure
    //! (depending on formulation) in a solution vector
    static const int pressure0Idx = PVOffset + 0;
    //! Index of the saturation of the non-wetting/wetting phase
    static const int saturation0Idx = PVOffset + 1;

    // indices of the equations

    //! Index of the continuity equation of the first phase
    static const int conti0EqIdx = PVOffset + 0;
};
} // namespace Ewoms

#endif
