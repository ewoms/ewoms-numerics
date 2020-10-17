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
 * \copydoc Ewoms::RichardsExtensiveQuantities
 */
#ifndef EWOMS_RICHARDS_EXTENSIVE_QUANTITIES_HH
#define EWOMS_RICHARDS_EXTENSIVE_QUANTITIES_HH

#include "richardsproperties.hh"

#include <ewoms/numerics/common/multiphasebaseextensivequantities.hh>

namespace Ewoms {

/*!
 * \ingroup RichardsModel
 * \ingroup ExtensiveQuantities
 *
 * \brief Calculates and stores the data which is required to
 *        calculate the flux of fluid over a face of a finite volume.
 */
template <class TypeTag>
class RichardsExtensiveQuantities
    : public MultiPhaseBaseExtensiveQuantities<TypeTag>
{
};

} // namespace Ewoms

#endif
