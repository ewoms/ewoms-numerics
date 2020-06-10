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
 * \copydoc Ewoms::FlashExtensiveQuantities
 */
#ifndef EWOMS_FLASH_EXTENSIVE_QUANTITIES_HH
#define EWOMS_FLASH_EXTENSIVE_QUANTITIES_HH

#include "flashproperties.hh"

#include <ewoms/numerics/common/multiphasebaseextensivequantities.hh>
#include <ewoms/numerics/common/energymodule.hh>
#include <ewoms/numerics/common/diffusionmodule.hh>

namespace Ewoms {

/*!
 * \ingroup FlashModel
 * \ingroup ExtensiveQuantities
 *
 * \brief This template class contains the data which is required to
 *        calculate all fluxes of components over a face of a finite
 *        volume for the compositional multi-phase model based on
 *        flash calculations.
 *
 * This means pressure and concentration gradients, phase densities at
 * the integration point, etc.
 */
template <class TypeTag>
class FlashExtensiveQuantities
    : public MultiPhaseBaseExtensiveQuantities<TypeTag>
    , public EnergyExtensiveQuantities<TypeTag, GET_PROP_VALUE(TypeTag, EnableEnergy)>
    , public DiffusionExtensiveQuantities<TypeTag, GET_PROP_VALUE(TypeTag, EnableDiffusion)>
{
    typedef MultiPhaseBaseExtensiveQuantities<TypeTag> ParentType;

    typedef GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    enum { enableDiffusion = GET_PROP_VALUE(TypeTag, EnableDiffusion) };
    typedef Ewoms::DiffusionExtensiveQuantities<TypeTag, enableDiffusion> DiffusionExtensiveQuantities;

    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    typedef Ewoms::EnergyExtensiveQuantities<TypeTag, enableEnergy> EnergyExtensiveQuantities;

public:
    /*!
     * \copydoc MultiPhaseBaseExtensiveQuantities::update
     */
    void update(const ElementContext& elemCtx, unsigned scvfIdx, unsigned timeIdx)
    {
        ParentType::update(elemCtx, scvfIdx, timeIdx);
        DiffusionExtensiveQuantities::update_(elemCtx, scvfIdx, timeIdx);
        EnergyExtensiveQuantities::update_(elemCtx, scvfIdx, timeIdx);
    }

    /*!
     * \copydoc MultiPhaseBaseExtensiveQuantities::updateBoundary
     */
    template <class Context, class FluidState>
    void updateBoundary(const Context& context,
                        unsigned bfIdx,
                        unsigned timeIdx,
                        const FluidState& fluidState)
    {
        ParentType::updateBoundary(context, bfIdx, timeIdx, fluidState);
        DiffusionExtensiveQuantities::updateBoundary_(context, bfIdx, timeIdx, fluidState);
        EnergyExtensiveQuantities::updateBoundary_(context, bfIdx, timeIdx, fluidState);
    }
};

} // namespace Ewoms

#endif
