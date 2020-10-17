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
 * \copydoc Ewoms::FlashBoundaryRateVector
 */
#ifndef EWOMS_FLASH_BOUNDARY_RATE_VECTOR_HH
#define EWOMS_FLASH_BOUNDARY_RATE_VECTOR_HH

#include "flashproperties.hh"

#include <ewoms/numerics/common/energymodule.hh>
#include <ewoms/common/valgrind.hh>

namespace Ewoms {

/*!
 * \ingroup FlashModel
 *
 * \brief Implements a boundary vector for the fully implicit
 *        compositional multi-phase model which is based on flash
 *        calculations.
 */
template <class TypeTag>
class FlashBoundaryRateVector : public GET_PROP_TYPE(TypeTag, RateVector)
{
    using ParentType = GET_PROP_TYPE(TypeTag, RateVector);
    using ExtensiveQuantities = GET_PROP_TYPE(TypeTag, ExtensiveQuantities);
    using Evaluation = GET_PROP_TYPE(TypeTag, Evaluation);
    using Indices = GET_PROP_TYPE(TypeTag, Indices);

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };

    using EnergyModule = Ewoms::EnergyModule<TypeTag, enableEnergy>;
    using Toolbox = Ewoms::MathToolbox<Evaluation>;

public:
    FlashBoundaryRateVector() : ParentType()
    {}

    /*!
     * \copydoc
     * ImmiscibleBoundaryRateVector::ImmiscibleBoundaryRateVector(Scalar)
     */
    FlashBoundaryRateVector(const Evaluation& value) : ParentType(value)
    {}

    /*!
     * \copydoc ImmiscibleBoundaryRateVector::ImmiscibleBoundaryRateVector(const
     * ImmiscibleBoundaryRateVector& )
     */
    FlashBoundaryRateVector(const FlashBoundaryRateVector& value) = default;
    FlashBoundaryRateVector& operator=(const FlashBoundaryRateVector& value) = default;

    /*!
     * \copydoc ImmiscibleBoundaryRateVector::setFreeFlow
     */
    template <class Context, class FluidState>
    void setFreeFlow(const Context& context,
                     unsigned bfIdx,
                     unsigned timeIdx,
                     const FluidState& fluidState)
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
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
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

            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                Evaluation molarity;
                if (fluidState.pressure(phaseIdx) > insideIntQuants.fluidState().pressure(phaseIdx)) {
                    if (focusDofIdx == interiorDofIdx)
                        molarity = fluidState.molarity(phaseIdx, compIdx);
                    else
                        molarity = Ewoms::getValue(fluidState.molarity(phaseIdx, compIdx));
                }
                else if (focusDofIdx == interiorDofIdx)
                    molarity = insideIntQuants.fluidState().molarity(phaseIdx, compIdx);
                else
                    molarity = Ewoms::getValue(insideIntQuants.fluidState().molarity(phaseIdx, compIdx));

                // add advective flux of current component in current
                // phase
                (*this)[conti0EqIdx + compIdx] += extQuants.volumeFlux(phaseIdx)*molarity;
            }

            if (enableEnergy) {
                Evaluation specificEnthalpy;
                if (fluidState.pressure(phaseIdx) > insideIntQuants.fluidState().pressure(phaseIdx)) {
                    if (focusDofIdx == interiorDofIdx)
                        specificEnthalpy = fluidState.enthalpy(phaseIdx);
                    else
                        specificEnthalpy = Ewoms::getValue(fluidState.enthalpy(phaseIdx));
                }
                else if (focusDofIdx == interiorDofIdx)
                    specificEnthalpy = insideIntQuants.fluidState().enthalpy(phaseIdx);
                else
                    specificEnthalpy = Ewoms::getValue(insideIntQuants.fluidState().enthalpy(phaseIdx));

                Evaluation enthalpyRate = density*extQuants.volumeFlux(phaseIdx)*specificEnthalpy;
                EnergyModule::addToEnthalpyRate(*this, enthalpyRate);
            }
        }

        // thermal conduction
        EnergyModule::addToEnthalpyRate(*this, EnergyModule::thermalConductionRate(extQuants));

#ifndef NDEBUG
        for (unsigned i = 0; i < numEq; ++i) {
            Ewoms::Valgrind::CheckDefined((*this)[i]);
        }
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
