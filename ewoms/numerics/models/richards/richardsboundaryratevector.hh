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
 * \copydoc Ewoms::RichardsBoundaryRateVector
 */
#ifndef EWOMS_RICHARDS_BOUNDARY_RATE_VECTOR_HH
#define EWOMS_RICHARDS_BOUNDARY_RATE_VECTOR_HH

#include <ewoms/common/valgrind.hh>
#include <ewoms/material/constraintsolvers/ncpflash.hh>

#include "richardsintensivequantities.hh"

namespace Ewoms {

/*!
 * \ingroup RichardsModel
 *
 * \brief Implements a boundary vector for the fully implicit Richards model.
 */
template <class TypeTag>
class RichardsBoundaryRateVector : public GET_PROP_TYPE(TypeTag, RateVector)
{
    using ParentType = GET_PROP_TYPE(TypeTag, RateVector);
    using ExtensiveQuantities = GET_PROP_TYPE(TypeTag, ExtensiveQuantities);
    using Evaluation = GET_PROP_TYPE(TypeTag, Evaluation);
    using Indices = GET_PROP_TYPE(TypeTag, Indices);

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { contiEqIdx = Indices::contiEqIdx };
    enum { liquidPhaseIdx = GET_PROP_VALUE(TypeTag, LiquidPhaseIndex) };

    using Toolbox = Ewoms::MathToolbox<Evaluation>;

public:
    RichardsBoundaryRateVector() : ParentType()
    {}

    /*!
     * \copydoc
     * ImmiscibleBoundaryRateVector::ImmiscibleBoundaryRateVector(Scalar)
     */
    RichardsBoundaryRateVector(const Evaluation& value)
        : ParentType(value)
    {}

    /*!
     * \copydoc ImmiscibleBoundaryRateVector::ImmiscibleBoundaryRateVector(const
     * ImmiscibleBoundaryRateVector& )
     */
    RichardsBoundaryRateVector(const RichardsBoundaryRateVector& value) = default;
    RichardsBoundaryRateVector& operator=(const RichardsBoundaryRateVector& value) = default;

    /*!
     * \copydoc ImmiscibleBoundaryRateVector::setFreeFlow
     */
    template <class Context, class FluidState>
    void setFreeFlow(const Context& context, unsigned bfIdx, unsigned timeIdx, const FluidState& fluidState)
    {
        ExtensiveQuantities extQuants;
        extQuants.updateBoundary(context, bfIdx, timeIdx, fluidState);
        const auto& insideIntQuants = context.intensiveQuantities(bfIdx, timeIdx);
        unsigned focusDofIdx = context.focusDofIndex();
        unsigned interiorDofIdx = context.interiorScvIndex(bfIdx, timeIdx);

        ////////
        // advective fluxes of all components in all phases
        ////////
        (*this) = Evaluation(0.0);

        unsigned phaseIdx = liquidPhaseIdx;
        Evaluation density;
        if (fluidState.pressure(phaseIdx) > insideIntQuants.fluidState().pressure(phaseIdx)) {
            if (focusDofIdx == interiorDofIdx)
                density = fluidState.density(phaseIdx);
            else
                density = Ewoms::getValue(fluidState.density(phaseIdx));
        }
        else if (focusDofIdx == interiorDofIdx)
            density = insideIntQuants.fluidState().density(phaseIdx);
        else
            density = Ewoms::getValue(insideIntQuants.fluidState().density(phaseIdx));

        // add advective flux of current component in current
        // phase
        (*this)[contiEqIdx] += extQuants.volumeFlux(phaseIdx) * density;

#ifndef NDEBUG
        for (unsigned i = 0; i < numEq; ++i) {
            Ewoms::Valgrind::CheckDefined((*this)[i]);
        }
        Ewoms::Valgrind::CheckDefined(*this);
#endif
    }

    /*!
     * \copydoc ImmiscibleBoundaryRateVector::setInFlow
     */
    template <class Context, class FluidState>
    void setInFlow(const Context& context,
                   unsigned bfIdx,
                   unsigned timeIdx,
                   const FluidState& fluidState)
    {
        this->setFreeFlow(context, bfIdx, timeIdx, fluidState);

        // we only allow fluxes in the direction opposite to the outer unit normal
        for (unsigned eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            Evaluation& val = this->operator[](eqIdx);
            val = Toolbox::min(0.0, val);
        }
    }

    /*!
     * \copydoc ImmiscibleBoundaryRateVector::setOutFlow
     */
    template <class Context, class FluidState>
    void setOutFlow(const Context& context,
                    unsigned bfIdx,
                    unsigned timeIdx,
                    const FluidState& fluidState)
    {
        this->setFreeFlow(context, bfIdx, timeIdx, fluidState);

        // we only allow fluxes in the same direction as the outer unit normal
        for (unsigned eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            Evaluation& val = this->operator[](eqIdx);
            val = Toolbox::max(0.0, val);
        }
    }

    /*!
     * \copydoc ImmiscibleBoundaryRateVector::setNoFlow
     */
    void setNoFlow()
    { (*this) = Evaluation(0.0); }
};

} // namespace Ewoms

#endif
